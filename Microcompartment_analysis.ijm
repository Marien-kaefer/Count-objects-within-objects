/*


												- Written by Marie Held [mheldb@liverpool.ac.uk] January 2024
												  Liverpool CCI (https://cci.liverpool.ac.uk/)
________________________________________________________________________________________________________________________

BSD 2-Clause License

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*
*/

//get input parameters
#@ String(value="Please specify the files and parameters to be applied for the processing.", visibility="MESSAGE") message
#@ File (label = "File for segmentation:", style = "open") segmentation_inputFile
#@ String(label = "Is the file for segmentation fluorescence light?", choices={"No", "Yes"}, style="radioButtonHorizontal") segmentation_file_choice
#@ Integer(label="If so, which channel should be segmented?" , value = 1) fluorescence_channel_number_to_segment
#@ File (label = "Fluorescence file:", style = "open") fluorescence_inputFile
#@ Boolean(label="Microcompartment Counting") analysis_choice_microcomp_count
#@ Boolean(label="Heat Map Generation") analysis_choice_heatmap
#@ Boolean(label="Measure all ROIs (all parameters)") analysis_choice_measure_all
//#@ Boolean(label="First Microcompartment Genesis Site") analysis_choice_first_genesis_site
#@ File(label = "Labkit classifier file (if using): ", style = "file") classifier_file
#@ Integer(label="Ch1 Prominence: " , value = 1000) Ch1_prominence
#@ Integer(label="Ch2_prominence: " , value = 1000) Ch2_prominence
#@ Double(label="Bacteria size minimum: ", value = 0.4) bacteria_size_minimum
#@ Double(label="Bacteria size maximum: ", value = 3.0) bacteria_size_maximum
#@ Double(label="Bacteria circularity minimum: ", value = 0.0) bacteria_circularity_minimum
#@ Double(label="Bacteria circularity maximum: ", value = 0.9) bacteria_circularity_maximum
#@ Double(label="Intensity histogram kurtosis cutoff value: ", value = 3.0) kurtosis_cutoff

start = getTime(); 
print(TimeStamp() + " The analysis is running");

pre_clean_up();
//variable definitions and calculations
file_path = File.getDirectory(segmentation_inputFile);
segmentation_title = File.getName(segmentation_inputFile);

fluorescence_title = File.getName(fluorescence_inputFile); 
segmentation_title_no_ext = file_name_remove_extension(segmentation_title); 
//print("Segmentation title: " + segmentation_title_no_ext); 
fluorescence_title_no_ext = file_name_remove_extension(fluorescence_title); 
//print("Fluorescence title: " + fluorescence_title_no_ext); 
file_base_title = substring(segmentation_title_no_ext, 0, (lengthOf(segmentation_title_no_ext)-2));
ROI_set_title = file_base_title + "_bacteriaROIs.zip"; 
ROI_set_combined_title = file_base_title + "_bacteriaROIs-combined.roi";

run("Bio-Formats Windowless Importer", "open=[" + file_path + File.separator + fluorescence_title + "]");
getDimensions(width, height, channels, slices, frames);
SR_channel_delta = channels / 2 ;  
close(); 


if (segmentation_file_choice == "No") {
	TL_scaled_title = scaling_TL_image(file_path, segmentation_title, fluorescence_title);
	TL_scaled_mask = segment_TL_image(TL_scaled_title, classifier_file, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum); 
	
	ROI_set_title = determine_and_apply_xy_offset(TL_scaled_mask, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum, file_path, segmentation_title_no_ext, ROI_set_title, ROI_set_combined_title);
}
else if (segmentation_file_choice == "Yes") {
	ROI_set_title = segment_FL_image(file_path, fluorescence_title, fluorescence_channel_number_to_segment, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum, ROI_set_title, ROI_set_combined_title);
}

if (ROI_set_title != "NaN"){
	//setBatchMode("hide");
	if (analysis_choice_microcomp_count == true){
		print("Counting spots.");
		spot_counting(fluorescence_title, file_path, file_base_title, Ch1_prominence, Ch2_prominence, SR_channel_delta); 
		print("Counting spots finished."); 
	}
	
	if (analysis_choice_heatmap == true){
		print("Creating spot heat map.");
		spot_heatmap(fluorescence_title, file_path, ROI_set_combined_title, file_base_title, Ch1_prominence, Ch2_prominence, ROI_set_title, kurtosis_cutoff, SR_channel_delta);
		print("Heatmap generation finished."); 
	}
	
		if (analysis_choice_measure_all == true){
		print("Measuring all parameters for all ROIs.");
		measure_all(fluorescence_title, file_path, ROI_set_combined_title, file_base_title, Ch1_prominence, Ch2_prominence, ROI_set_title, kurtosis_cutoff);
		print("Measuring all parameters finished."); 
	}
	
	
	//setBatchMode("show");
	
}

//let user know the process has finished and how long it took
stop = getTime(); 
duration = stop - start;
duration_String = duration_conversion(duration);
print("The analysis took " + duration_conversion(duration));
print(TimeStamp() + ": Processing of [" + segmentation_title + "] complete.");
beep();


function pre_clean_up(){
	setForegroundColor(255, 255, 255);
	setBackgroundColor(0, 0, 0);
	run("Options...", "iterations=1 count=3 black");
	close("*");
	roiManager("reset");
	run("Clear Results");
}

function file_name_remove_extension(file_title){
	dotIndex = lastIndexOf(file_title, "." ); 
	file_name_without_extension = substring(file_title, 0, dotIndex );
	//print( "Name without extension: " + file_name_without_extension );
	return file_name_without_extension;
}

function scaling_TL_image(file_path, segmentation_title, fluorescence_title){
	run("Bio-Formats Windowless Importer", "open=[" + file_path + File.separator + segmentation_title + "]");
	run("Bio-Formats Windowless Importer", "open=[" + file_path + File.separator + fluorescence_title + "]");
	
	selectWindow(segmentation_title); 
	getDimensions(width, height, channels, slices, frames);
	TL_width = width; 
	TL_height = height; 
	getPixelSize(unit, pixelWidth, pixelHeight);
	TL_pixelWidth = pixelWidth;
	TL_pixelHeight = pixelHeight; 
	
	selectWindow(fluorescence_title); 
	getDimensions(width, height, channels, slices, frames);
	SIM_width = width; 
	SIM_height = height; 
	getPixelSize(unit, pixelWidth, pixelHeight);
	SIM_pixelWidth = pixelWidth;
	SIM_pixelHeight = pixelHeight; 
		
	SIM_sampling_rate_factor = SIM_width / TL_width; 
	//print("SIM sampling rate factor: " + SIM_sampling_rate_factor);
	optovar_factor = (TL_pixelWidth / SIM_pixelWidth) / SIM_sampling_rate_factor; 
	//print("Optovar factor: " + optovar_factor);
	
	if (optovar_factor == 1) {
		selectWindow(segmentation_title); 
		run("Duplicate...", " ");
		rename("temp_title"); 
} 
	else if(optovar_factor > 1) {  //for optovar > 1, i.e., SIM image optovar > TL optovar 
		rectangle_width = floor(TL_width / optovar_factor);
		//print("Rectangle width " + rectangle_width); 
		rectangle_height = floor(TL_height / optovar_factor);
		//print("Rectangle height " + rectangle_height); 
		rectangle_x = (TL_width - rectangle_width)/2;
		//print("Rectangle x " + rectangle_x); 
		rectangle_y = (TL_height - rectangle_height)/2;
		//print("Rectangle x " + rectangle_x);  
		selectWindow(segmentation_title); 
		makeRectangle(rectangle_x, rectangle_y, rectangle_width , rectangle_height);
		run("Duplicate...", " ");
		rename("temp_title"); 
	}
	else  //for optovar < 1, i.e., SIM image optovar < TL optovar 
	{
		rectangle_width = floor(SIM_width / optovar_factor);
		//print("Rectangle width " + rectangle_width); 
		rectangle_height = floor(SIM_height / optovar_factor);
		//print("Rectangle height " + rectangle_height); 
		rectangle_x = (SIM_width - rectangle_width)/2;
		//print("Rectangle x " + rectangle_x); 
		rectangle_y = (SIM_height - rectangle_height)/2;
		//print("Rectangle x " + rectangle_x);  
		selectWindow(fluorescence_title); 
		makeRectangle(rectangle_x, rectangle_y, rectangle_width , rectangle_height);run("Duplicate...", " ");
		rename("temp_title"); 
	}
	selectWindow("temp_title"); 
	run("Size...", "width=" + SIM_width + " height=" + SIM_height + " depth=1 constrain average interpolation=None");
	saveAs("Tiff", file_path + File.separator + segmentation_title_no_ext + "_scaled.tif");
	TL_scaled_title = getTitle; 
	return TL_scaled_title; 
}

function segment_TL_image(TL_scaled_title, classifier_file, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum){
	selectWindow(fluorescence_title); 
	getDimensions(width, height, channels, slices, frames);
	SIM_width = width; 
	SIM_height = height; 
	
	satifaction_score = "No";
    do {
    	selectWindow(segmentation_title); 
		run("Duplicate...", " ");
		TL_duplicate_ID = getImageID();
		
		segmentation_choice_array = newArray("Labkit (Mengru)", "Minimum (Ping)", "Variance (Kuo)", "Median & Minimum (Jinlu)"); 
		Dialog.create("Choose segmentation method");
		Dialog.addMessage("Please select a segmentation method.");	
		Dialog.addChoice("Segmentation approach: ", segmentation_choice_array);
		Dialog.show();
	
		user_segmentation_choice = Dialog.getChoice();

		if (user_segmentation_choice == "Labkit (Mengru)") {
			//Mengru
			run("Segment Image With Labkit", "segmenter_file=" + classifier_file + " use_gpu=false");
			setAutoThreshold("Default dark");
			//run("Threshold...");
			run("Convert to Mask"); 	
			TL_mask_ID = getImageID();
		}
		
		else if (user_segmentation_choice == "Minimum (Ping)") {
			//PING
			run("Enhance Contrast", "saturated=0.35");
			run("Duplicate...", " ");
			raw_duplicate_ID = getImageID();
			raw_duplicate_title = getTitle(); 
			run("Duplicate...", " ");
			raw_Gauss_ID = getImageID();
			raw_Gauss_title = getTitle();
			run("Gaussian Blur...", "sigma=20");
			imageCalculator("Subtract create 32-bit", raw_duplicate_title, raw_Gauss_title);
			BG_subtracted_ID = getImageID();
			selectImage(raw_duplicate_ID);
			close(); 
			selectImage(raw_Gauss_ID);
			close(); 

			selectImage(BG_subtracted_ID);
			run("Median...", "radius=3");
			run("Minimum...", "radius=0.5");
			//run("Threshold...");
			setAutoThreshold("Otsu");
			setOption("BlackBackground", true);
			run("Convert to Mask");			
			TL_mask_ID = getImageID();			
		}
		
		else if (user_segmentation_choice == "Variance (Kuo)") { 
			run("Enhance Contrast", "saturated=0.35");
			run("Duplicate...", " ");
			run("Subtract Background...", "rolling=20 light");
			run("Median...", "radius=2");
			run("Variance...", "radius=2");
			setAutoThreshold("Huang");
			setOption("BlackBackground", true);
			run("Convert to Mask");			
			TL_mask_ID = getImageID();
			run("Dilate");
			run("Dilate");
			run("Dilate");
		}
		else if (user_segmentation_choice == "Median & Minimum (Jinlu)"){
			run("Enhance Contrast", "saturated=0.35");
			run("Duplicate...", " ");
			run("Median...", "radius=2");
			run("Minimum...", "radius=2");
			//run("Threshold...");
			setAutoThreshold("Default");
			setOption("BlackBackground", true);
			run("Convert to Mask");	
			run("Fill Holes");
			TL_mask_ID = getImageID();
		}
		
		run("Analyze Particles...", "size=" + bacteria_size_minimum + "-" + bacteria_size_maximum + " circularity=" + bacteria_circularity_minimum + "-" + bacteria_circularity_maximum + " exclude clear add");
		selectImage(TL_duplicate_ID);
		roiManager("Show All without labels");	

		satisfaction_score_array = newArray("Yes", "No"); 
		Dialog.create("Segmentation Satisfaction");
		Dialog.addMessage("Are you satisfied with the segmentation result? If ''No'' is selected, you will presented with the segmentation choice again.");
		Dialog.addChoice("Satisfaction score: ", satisfaction_score_array);
		Dialog.show();

		satifaction_score = Dialog.getChoice();
		if (satifaction_score == "No"){
			selectImage(TL_mask_ID);
			close();
			selectImage(TL_duplicate_ID);
			close();
			roiManager("reset");
			run("Select None");
			run("Hide Overlay");
		}
		
   } while (satifaction_score != "Yes");

	roiManager("reset");
	run("Select None");
	run("Hide Overlay");
	selectImage(TL_mask_ID);
	run("Size...", "width=" + SIM_width + " height=" + SIM_height + " depth=1 constrain average interpolation=None");
	saveAs("Tiff", file_path + File.separator + segmentation_title_no_ext + "_scaled_mask.tif");
	TL_scaled_mask = getTitle();
	
	return TL_scaled_mask; 
	
	}

function segment_FL_image(file_path, fluorescence_title, fluorescence_channel_number_to_segment, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum, ROI_set_title, ROI_set_combined_title){
	
	run("Bio-Formats Windowless Importer", "open=[" + file_path + File.separator + fluorescence_title + "]");
	getDimensions(width, height, channels, slices, frames);
	if ((fluorescence_channel_number_to_segment<0 || fluorescence_channel_number_to_segment>channels)){ 
		waitForUser("Something went wrong.", "The fluorescence image does not have the specified channel. Please start the analysis again.");
		ROI_set_title = "NaN"; 
	}
	else {
	 
	 	run("Duplicate...", "duplicate channels=" + fluorescence_channel_number_to_segment);
		run("Median...", "radius=2");
		setAutoThreshold("Default dark");
		//run("Threshold...");
		run("Convert to Mask");
		run("Analyze Particles...", "size=" + bacteria_size_minimum + "-" + bacteria_size_maximum + " circularity=" + bacteria_circularity_minimum + "-" + bacteria_circularity_maximum + " exclude clear add");
		fluorescence_segmentation_mask_title = segmentation_title_no_ext + "_fluorescence_mask_Ch" +  fluorescence_channel_number_to_segment; 
		roiManager("show none");
		saveAs("Tiff", file_path + File.separator + fluorescence_segmentation_mask_title + ".tif");
		fluorescence_segmentation_mask_title = getTitle(); 
		
		roiManager("deselect");
		roiManager("save", file_path + File.separator + ROI_set_title);
		roiManager("combine");
		roiManager("add");
		ROI_count = roiManager("count") - 1;
		roiManager("select", ROI_count);
		roiManager("save", file_path + File.separator + ROI_set_combined_title);	
	}
	return ROI_set_title;
}

function determine_and_apply_xy_offset(TL_scaled_mask, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum, file_path, segmentation_title_no_ext, ROI_set_title, ROI_set_combined_title){
	selectWindow(TL_scaled_mask); 
	setAutoThreshold("Default dark");
	//run("Threshold...");
	run("Convert to Mask"); 
	run("Analyze Particles...", "size=" + bacteria_size_minimum + "-" + bacteria_size_maximum + " circularity=" + bacteria_circularity_minimum + "-" + bacteria_circularity_maximum + " exclude clear add");

	amend_ROI_set(TL_scaled_mask);	
	
	for (i = 0; i < roiManager("count"); i++) {
		roiManager("select", i);
		roiManager("Rename", IJ.pad((i), 3));
	}

	selectWindow(fluorescence_title); 
	run("Enhance Contrast", "saturated=0.35");
	roiManager("Show All");

	//ROI_to_move = Math.ceil(roiManager("count")/2); 
	//roiManager("Select", ROI_to_move);
	waitForUser("Select an ROI by clicking on its number label but do not move the ROI yet. Then confirm with 'OK'.");
	ROI_to_move = roiManager("index");
	getSelectionBounds(x, y, w, h);
	original_position_x = x; 
	original_position_y = y;
	waitForUser("Move the selected ROI to the correct position then click 'OK'.");
	
	roiManager("Update");
	getSelectionBounds(x, y, w, h);
	adjusted_position_x = x; 
	adjusted_position_y = y;
	
	roiManager("Select", ROI_to_move);
	setSelectionLocation(original_position_x, original_position_y);
	
	ROI_movement_x = adjusted_position_x - original_position_x; 
	ROI_movement_y = adjusted_position_y - original_position_y; 
	//print("Movement in x " + ROI_movement_x + " and y " + ROI_movement_y); 
	
	//move all ROIs to the correct position
	number_of_ROI = roiManager("count");
	if (number_of_ROI==0)
	  exit("The ROI Manager is empty");
	for (i=0; i<number_of_ROI; i++) {
	  roiManager('select', i);
	  getSelectionBounds(x, y, w, h);
	  setSelectionLocation(x+ROI_movement_x, y+ROI_movement_y);
	  roiManager('update');
	}
	roiManager("deselect");
	roiManager("save", file_path + File.separator + ROI_set_title);
	roiManager("combine");
	roiManager("add");
	ROI_count = roiManager("count") - 1;
	roiManager("select", ROI_count);
	roiManager("save", file_path + File.separator + ROI_set_combined_title);
	
	selectWindow(TL_scaled_title); 
	run("Translate...", "x=" + ROI_movement_x +" y=" + ROI_movement_y +" interpolation=None");
	saveAs("Tiff", file_path + File.separator + segmentation_title_no_ext + "_scaled.tif");
	close(); 
	selectWindow(TL_scaled_mask); 
	close(); 
	roiManager("reset");
	run("Select None");	
	return ROI_set_title; 
}

function amend_ROI_set(TL_scaled_mask){
	selectWindow(TL_scaled_mask); 
	roiManager("deselect");
	waitForUser("In case any ROIs should be rejected for analysis, select the ROI by clicking on the numerical label in the image and then delete. \n Click ''OK'' when done."); 
}

// ##################### SPOT COUNTING ##################### 
function spot_counting(fluorescence_title, file_path, file_base_title, Ch1_prominence, Ch2_prominence, SR_channel_delta){
	selectWindow(fluorescence_title); 
	run("Set Measurements...", "area mean shape kurtosis redirect=None decimal=3");
	roiManager("Open", file_path + File.separator + ROI_set_title);
	roiManager("Deselect");
	roiManager("multi-measure measure_all");
	
	number_of_bacteria = roiManager("count") ; 
	Ch1_kurtosis = newArray(number_of_bacteria);
	Ch1_mean = newArray(number_of_bacteria);
	Ch1_class = newArray(number_of_bacteria); 
	Ch2_kurtosis = newArray(number_of_bacteria);
	Ch2_mean = newArray(number_of_bacteria);
	Ch2_class = newArray(number_of_bacteria); 
	
	for (i = 0; i < number_of_bacteria; i++) {
	
		Ch1_kurtosis[i] = getResult("Kurt", i);
		Ch2_kurtosis[i] = getResult("Kurt", i + (SR_channel_delta * number_of_bacteria));
		if (Ch1_kurtosis[i] > kurtosis_cutoff) {
			Ch1_class[i] = "spots";
		} 
		else{
			Ch1_class[i] = "homogeneous";
		}
		Ch1_mean[i] = getResult("Mean", i);
		Ch2_mean[i] = getResult("Mean", i + (SR_channel_delta * number_of_bacteria));
		if (Ch2_kurtosis[i] > kurtosis_cutoff) {
			Ch2_class[i] = "spots";
		} 
		else{
			Ch2_class[i] = "homogeneous";
		}
	}
	
	selectWindow("Results"); 
	run("Close");			
			
	Ch1_spot_count = newArray(number_of_bacteria);
	Ch2_spot_count = newArray(number_of_bacteria);
	//loopStop = 20; 
	loopStop = roiManager("count");
	
	// loop to count the number of spots per bacterium 
	
	selectWindow(fluorescence_title); 
	//run("Brightness/Contrast...");
	//run("Enhance Contrast", "saturated=0.01");
	for(i=0; i<loopStop; i++) {
		roiManager("select", i);
		 Stack.setChannel(1);
		//roiManager("rename", "Object " + i + 1);
		run("Find Maxima...", "noise="+Ch1_prominence+" output=[Count]");
		Ch1_spot_count[i] = getResult("Count", i);
		run("Find Maxima...", "noise="+Ch1_prominence+" output=[Point Selection]");
		run("Add Selection...");
	}

	selectWindow(fluorescence_title); 
	roiManager("Show None");
	run("Select None");
	run("Duplicate...", "ignore"); 
	saveAs("Tiff", file_path + File.separator + file_base_title + "_Ch1-w-spots.tif");
	close();
	selectWindow("Results"); 
	run("Close");
	
	run("Remove Overlay");
	selectWindow(fluorescence_title); 
	//run("Brightness/Contrast...");
	run("Enhance Contrast", "saturated=0.01");
	for(i=0; i<loopStop; i++) {
		roiManager("select", i);
		Stack.setChannel(SR_channel_delta + 1);
		run("Find Maxima...", "noise="+Ch2_prominence+" output=[Count]");
		Ch2_spot_count[i] = getResult("Count", i);
		run("Find Maxima...", "noise="+Ch2_prominence+" output=[Point Selection]");
		run("Add Selection...");
	}
	roiManager("Show None");
	run("Select None");
	run("Remove Overlay");
	run("Duplicate...", "ignore"); 
	saveAs("Tiff", file_path + File.separator + file_base_title + "_Ch2-w-spots.tif");
	close(); 
	run("Clear Results");
		
	//clean up
	if (isOpen(segmentation_title)){
		selectWindow(segmentation_title);
		run("Close");
	}
	roiManager("reset");
	assemble_and_save_results(number_of_bacteria, Ch1_spot_count, Ch1_mean, Ch1_kurtosis, Ch1_class, Ch2_spot_count, Ch2_mean, Ch2_kurtosis, Ch2_class); 

}

function assemble_and_save_results(number_of_bacteria, Ch1_spot_count, Ch1_mean, Ch1_kurtosis, Ch1_class, Ch2_spot_count, Ch2_mean, Ch2_kurtosis, Ch2_class){
	run("Clear Results");
	//write results into a new results window 
	for (i = 0; i < number_of_bacteria; i++) {
		setResult("Bacterium ID", i, IJ.pad(i,3));
		setResult("Ch1 spot count", i, Ch1_spot_count[i]);
		setResult("Ch1 mean", i, Ch1_mean[i]); 
		setResult("Ch1 Kurtosis", i, Ch1_kurtosis[i]);
		setResult("Ch1 Class", i , Ch1_class[i]); 
		setResult("Ch2 spot count", i, Ch2_spot_count[i]);  
		setResult("Ch2 mean", i, Ch2_mean[i]);
		setResult("Ch2 Kurtosis", i, Ch2_kurtosis[i]);
		setResult("Ch2 Class", i , Ch2_class[i]); 
	}
	saveAs("Results", file_path + File.separator + file_base_title + "_spot_count_results.csv");
	close("Results");
}

// ##################### HEATMAP GENERATION #####################  

function spot_heatmap(fluorescence_title, file_path, ROI_set_combined_title, file_base_title, Ch1_prominence, Ch2_prominence, ROI_set_title, kurtosis_cutoff, SR_channel_delta){
	//generate mask of bacteria bodies
	roiManager("reset");
	run("Select None");
	selectWindow(fluorescence_title); 
	getDimensions(width, height, channels, slices, frames);
	FL_channels = channels; 
	FL_width = width; 
	FL_height = height;
	newImage("Untitled", "8-bit black", FL_width, FL_height, 1);
	roiManager("Open", file_path + File.separator + ROI_set_combined_title);
	roiManager("select", 0);
	run("Add...", "value=255");
	bacteria_mask_file_title = file_base_title + "_objects";
	rename(bacteria_mask_file_title); 
	saveAs("TIFF", file_path + File.separator + bacteria_mask_file_title + ".tif");
	rename(bacteria_mask_file_title); 
	run("Select None");
	roiManager("reset");
	run("Set Measurements...", "fit kurtosis display redirect=None decimal=3");
	
	Dialog.create("Heatmap options");
	Dialog.addMessage("Please select the type(s) of heatmap to be created.");
	Dialog.addCheckbox("Spot peak locations", true);
	Dialog.addCheckbox("Spot sizes", false);
	Dialog.addCheckbox("Intensity sum (spots only)", false);
	Dialog.addCheckbox("Intensity sum (all bacteria)", false);
	Dialog.addNumber("Heatmap aspect ratio", 2.6, 1, 10, " ");
	Dialog.show();
	
	spot_peak_location_heatmap_choice = Dialog.getCheckbox();
	spot_size_heatmap_choice = Dialog.getCheckbox();
	intensity_sum_heatmap_choice = Dialog.getCheckbox();
	intensity_sum_all_heatmap_choice = Dialog.getCheckbox();
	aspect_ratio = Dialog.getNumber();
	spots_image_title = "spots-locations";
	spots_file_title = file_base_title + "_spots_locations";
	newImage(spots_file_title, "8-bit grayscale-mode", FL_width, FL_height, FL_channels, 1, 1);
	
	// ##################### HEATMAP PRE-PROCESSING ##################### 
	
	//set parameters for heatmap and create empty first heatmap image

	roiManager("open", file_path + File.separator + ROI_set_title);
	ROI_initial_count = roiManager("count");
	k = 1; 
	
	for (i = 0; i < channels; i+=SR_channel_delta) {
		channel_prominence_variable_name = "Ch" + k + "_prominence";
		
		//print(channel_prominence_variable_name); 		
		
		spots_file_title = file_base_title + "_spots_locations";
		
		if (spot_peak_location_heatmap_choice == true){ 	
			roiManager("reset");
			run("Select None"); 
			spots_file_title = spots_file_title + "_Ch" + k; 
			newImage(spots_file_title, "8-bit black", FL_width, FL_height, 1);
			spots_file_ID = getImageID();
			roiManager("reset");			
			print("Channel " + k + " - Spot location heatmap generation");
			points_ROI_file_title = file_base_title + "_Ch" + k + "_points.roi";
			selectWindow(fluorescence_title);
			Stack.setChannel(i+1);
			run("Select None");			
			run("Hide Overlay");		
			run("Find Maxima...", "prominence=" + channel_prominence_variable_name + " output=[Point Selection]");
			roiManager("Add"); 
			roiManager("select", 0);
			roiManager("save", file_path + File.separator + points_ROI_file_title);
			selectImage(spots_file_ID);
			Stack.setChannel(i+1);
			roiManager("select", 0);
			run("Draw", "slice");
			run("Subtract...", "value=254");
			setMinAndMax(0, 1);
			selectWindow(fluorescence_title);
			roiManager("reset");
			run("Select None");
			selectImage(spots_file_ID);
			//run("Hide Overlay");
			saveAs("Tiff", file_path + File.separator + spots_file_title + ".tif");
			rename(spots_file_title); 
			roiManager("open", file_path + File.separator + ROI_set_title);
			all_bacteria_choice = "no"; 
			heatmap_stack_title = "Spot_peak_location";
			heatmap_creation_loop(ROI_initial_count, spots_file_title, bacteria_mask_file_title, bacteria_size_minimum, file_path, i, heatmap_stack_title, all_bacteria_choice, aspect_ratio);
			run("Select None");
		}
		if (spot_size_heatmap_choice == true){
			print("Channel " + k + " - Spot size heatmap generation");
			maxima_in_tolerance_file_title = "Channel " + k + "_maxima_in_tolerance";
			selectWindow(fluorescence_title);
			roiManager("deselect");
			run("Select None");
			Stack.setChannel(i+1);	
			//print("Duplicating channel: " + (i+1)); 
			run("Duplicate...", "duplicate channels="+(i+1));
			run("Find Maxima...", "prominence=" + channel_prominence_variable_name + " output=[Maxima Within Tolerance]"); 
			run("Subtract...", "value=254");
			setMinAndMax(0, 1);
			rename(maxima_in_tolerance_file_title);		
			all_bacteria_choice = "no"; 
			heatmap_stack_title = "Maxima_in_tolerance";
			heatmap_creation_loop(ROI_initial_count,maxima_in_tolerance_file_title, bacteria_mask_file_title, bacteria_size_minimum, file_path, i, heatmap_stack_title, all_bacteria_choice, aspect_ratio);
			roiManager("deselect");
		}
		
		if (intensity_sum_heatmap_choice == true){
			print("Channel " + k + " - Intensity heatmap generation - spots only");
			
			selectWindow(fluorescence_title);
			Stack.setChannel(i+1);
			all_bacteria_choice = "no"; 
			heatmap_stack_title = "Intensity";
			heatmap_creation_loop(ROI_initial_count, fluorescence_title, bacteria_mask_file_title, bacteria_size_minimum, file_path, i, heatmap_stack_title, all_bacteria_choice, aspect_ratio);
			roiManager("deselect");
			//rescale heatmap to between 0 and 1
			//getMinAndMax(min, max);
			//print(min); 
			//print(max); 		
			//run("Subtract...", "value=" + min);
			//run("Divide...", "value=" + (max - min));
			//resetMinAndMax;
			//saveAs("TIFF", file_path + File.separator + new_heatmap + ".tif");
		}		
		if (intensity_sum_all_heatmap_choice == true){
			print("Channel " + k + " - Intensity heatmap generation all bacteria");
			
			selectWindow(fluorescence_title);
			Stack.setChannel(i+1);
			all_bacteria_choice = "yes"; 
			heatmap_stack_title = "Intensity_spots_and_homogeneous";
			heatmap_creation_loop(ROI_initial_count, fluorescence_title, bacteria_mask_file_title, bacteria_size_minimum, file_path, i, heatmap_stack_title, all_bacteria_choice, aspect_ratio);
			roiManager("deselect");
			//rescale heatmap to between 0 and 1
			//getMinAndMax(min, max);
			//print(min); 
			//print(max); 		
			//run("Subtract...", "value=" + min);
			//run("Divide...", "value=" + (max - min));
			//resetMinAndMax;
			//saveAs("TIFF", file_path + File.separator + new_heatmap + ".tif");
		}	
	k += 1; 
	}
}


function heatmap_creation_loop(ROI_initial_count, input_data_image, bacteria_mask_file_title, bacteria_size_minimum, file_path, i, heatmap_stack_title, all_bacteria_choice, aspect_ratio){
		counter = 0;
		
		//create an empty stack  
		heatmap_height = 26; 
		heatmap_width = Math.ceil(heatmap_height * aspect_ratio);
		newImage("Heatmap Stack", "16-bit black", heatmap_width, heatmap_height, ROI_initial_count);
		heatmap_stack_ID = getImageID();
		//loop over all ROIs
		//for (r = 0; r < 20; r++) {
		for (r = 0; r < ROI_initial_count; r++) {
			//print("Roi index: " + r); 
			open_images = getList("image.titles");
			//Array.print(open_images); 	
			
			ROI_of_choice = r;
			//print("Processing object " + (r+1) +"/" + ROI_initial_count + ", Channel " + (i+1)); 
			//print("input data image: " + input_data_image); 
			//print("old heat map: " + old_heatmap_file_name); 
		
			//duplicate signal of input channel 
			selectWindow(fluorescence_title);
			roiManager("select", ROI_of_choice);
			object_angle = getValue("Angle");
			object_kurtosis = getValue("Kurt");
			//print("Object angle: " + object_angle);
			//print("Object kurtosis: " + object_kurtosis);			

			selectWindow(input_data_image);
			roiManager("select", ROI_of_choice);	
			run("Duplicate...", " duplicate channels=" + (i+1) + " title=Spots_ROI_CH" + k + "_object_" + r);		
			run("Clear Outside");
			//rotate duplicated bacterial spots image by fit ellipse angle
			run("Rotate... ", "angle=" + object_angle + " grid=1 interpolation=Bicubic enlarge"); 
			rotated_object = getTitle();		
			
			selectWindow(bacteria_mask_file_title); 
			roiManager("select", ROI_of_choice);
			run("Duplicate...", "title=Mask_ROI_" + r);
			run("Clear Outside");
			run("Rotate... ", "angle=" + object_angle + " grid=1 interpolation=Bicubic enlarge"); //rotate duplicated bacterium image by fit ellipse angle
			rotated_single_object_mask = getTitle();
									 
			roiManager("deselect");
			//create new ROI from mask
			selectWindow(rotated_single_object_mask);
			setAutoThreshold("Default dark no-reset");
			//run("Threshold...");
			//setThreshold(255, 255);
			run("Convert to Mask");
			roiManager("deselect");
			run("Select None");
			run("Analyze Particles...", "size="+ bacteria_size_minimum + "-Infinity add");
		
			selectWindow(rotated_object); 
			ROI_new_count =  roiManager("count") - 1;
			roiManager("select", ROI_new_count);
			run("Duplicate...", "duplicate title=final_rotated_spots");
			final_rotated_object_title = getTitle();
		
			run("Size...", "width=" + heatmap_width + " height=" + heatmap_height + " depth=1 average interpolation=Bicubic");
			resized_rotated_spots_title = getTitle();
			resized_rotated_spots_ID = getImageID();
			
			if (all_bacteria_choice == "no"){
				if (object_kurtosis < kurtosis_cutoff) {
					run("Set...", "value=0");
					counter = counter + 1;
					print("Rejected bacteria due to kurtosis: " + counter);  
				}
				else if (object_kurtosis == "NaN") {
					print("Rejected bacteria due to ROI located outside of the image area");  
				}
			}
			
			selectImage(resized_rotated_spots_ID);
			run("Select None");
			run("Copy");
			selectImage(heatmap_stack_ID); 
			setSlice(r + 1);
			run("Paste");
			run("Enhance Contrast", "saturated=0.35");
			
			selectWindow(rotated_object);
			close();
			selectWindow(rotated_single_object_mask);
			close();
			selectImage(resized_rotated_spots_ID); 
			close(); 
			
			roiManager("deselect");
			run("Select None");
			
		}
		selectImage(heatmap_stack_ID); 
		resetMinAndMax; 
		saveAs("TIFF", file_path + File.separator + file_base_title + "_Ch" + k + "_" + heatmap_stack_title + "_Heatmap_stack.tif");
		selectImage(heatmap_stack_ID); 
		run("Z Project...", "projection=[Sum Slices]");
		run("mpl-viridis");
		run("Enhance Contrast", "saturated=0.35");
		saveAs("TIFF", file_path + File.separator + file_base_title + "_Ch" + k + "_" + heatmap_stack_title + "_Heatmap_sum.tif");
		close(); 
		print("Heatmap creation loop complete."); 
	}


function measure_all(fluorescence_title, file_path, ROI_set_combined_title, file_base_title, Ch1_prominence, Ch2_prominence, ROI_set_title, kurtosis_cutoff){
	selectWindow(fluorescence_title);
	roiManager("Open", file_path + File.separator + ROI_set_title);
	roiManager("Deselect");
	run("Clear Results");
	run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis display redirect=None decimal=3");
	roiManager("multi-measure measure_all");
	saveAs("Results", file_path + File.separator + file_base_title + "_all_bacteria_measurements.csv");
	roiManager("reset");
	}

// ##################### GENERIC ##################### 

// set up time string for print statements
function TimeStamp(){
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	month = month + 1;
	TimeString ="["+year+"-";
	if (month<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+month + "-";
	if (dayOfMonth<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+dayOfMonth + " --- ";
	if (hour<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+hour+":";
	if (minute<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+minute+":";
	if (second<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+second + "]";
	return TimeString;
}

//convert time from ms to more appropriate time unit
function duration_conversion(duration){
	if (duration < 1000){
		duration_String = duration + " ms";
	} 
	else if (duration <60000){
		duration = duration / 1000;
		duration_String = d2s(duration, 0) + " s";
	}
	else if (duration <3600000){
		duration = duration / 60000;
		duration_String = d2s(duration, 1) +  "min";
	}
	else if (duration <86400000){
		duration = duration / 3600000;
		duration_String = d2s(duration, 0) + " hr";
	}
	else if (duration <604800000){
		duration = duration / 86400000;
		duration_String = d2s(duration, 0) + " d";
	}
	//print("Duration string: " + duration_String);	
	return duration_String;
}

close("*");
run("Clear Results");


