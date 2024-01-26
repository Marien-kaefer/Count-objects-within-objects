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
#@ String(value="Please select the transmitted light file of the file set you wish to process.", visibility="MESSAGE") message
#@ File (label = "Transmitted light file:", style = "open") inputFile
#@ Boolean(label="Microcompartment Counting") analysis_choice_microcomp_count
#@ Boolean(label="Heat Map") analysis_choice_heatmap
#@ String(label = "Path to TL classifier file: ", description = "Z:/private/Marie/Image_Analysis/2023-06-14-LIU-Mengru-count-spots-in-bacteria/input/Labkit/bacteria_TL_segmentation_2560.classifier") classifier_file
#@ Integer(label="Ch1 Prominence: " , value = 1000) Ch1_prominence
#@ Integer(label="Ch2_prominence: " , value = 1000) Ch2_prominence
#@ Double(label="Bacteria size minimum: ", value = 0.4) bacteria_size_minimum
#@ Double(label="Bacteria size maximum: ", value = 3.0) bacteria_size_maximum
#@ Double(label="Bacteria circularity minimum: ", value = 0.0) bacteria_circularity_minimum
#@ Double(label="Bacteria circularity maximum: ", value = 0.9) bacteria_circularity_maximum
#@ Double(label="Intensity histogram kurtosis cutoff value: ", value = 3.0) kurtosis_cutoff

start = getTime(); 
print("The analysis is running");

pre_clean_up();
//variable definitions and calculations
file_path = File.getDirectory(inputFile);
TL_title = File.getName(inputFile);
TL_title_no_ext = file_name_remove_extension(TL_title); 
file_base_title = substring(TL_title_no_ext, 0, (lengthOf(TL_title_no_ext)-2));
FL_title_no_ext = file_base_title + "-1"; 
FL_title = FL_title_no_ext + ".czi";
ROI_set_title = file_base_title + "_bacteriaROIs.zip"; 
ROI_set_combined_title = file_base_title + "_bacteriaROIs-combined.roi";

TL_scaled_title = scaling_TL_image(file_path, TL_title, FL_title);
TL_scaled_mask = segment_TL_image(TL_scaled_title, classifier_file); 
ROI_set_title = determine_and_apply_xy_offset(TL_scaled_mask, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum, file_path, TL_title_no_ext, ROI_set_title, ROI_set_combined_title);

setBatchMode("hide");
if (analysis_choice_microcomp_count == true){
	print("Counting spots.");
	spot_counting(FL_title, file_path, file_base_title, Ch1_prominence, Ch2_prominence); 
}

if (analysis_choice_heatmap == true){
	print("Creating spot heat map.");
	spot_heatmap(FL_title, file_path, ROI_set_combined_title, file_base_title, Ch1_prominence, Ch2_prominence, ROI_set_title, kurtosis_cutoff);
}
setBatchMode("show");

//let user know the process has finished and how long it took
stop = getTime(); 
duration = stop - start;
duration_String = duration_conversion(duration);
print("The analysis took " + duration_conversion(duration));
print(TimeStamp() + ": Processing of [" + TL_title + "] complete.");
beep();


function pre_clean_up(){
	setForegroundColor(255, 255, 255);
	setBackgroundColor(0, 0, 0);
	close("*");
	roiManager("reset");
	run("Clear Results");
}

function file_name_remove_extension(TL_title){
	dotIndex = lastIndexOf(TL_title, "." ); 
	file_name_without_extension = substring(TL_title, 0, dotIndex );
	//print( "Name without extension: " + file_name_without_extension );
	return file_name_without_extension;
}

function scaling_TL_image(file_path, TL_title, FL_title){
	run("Bio-Formats Windowless Importer", "open=[" + file_path + File.separator + TL_title + "]");
	run("Bio-Formats Windowless Importer", "open=[" + file_path + File.separator + FL_title + "]");
	
	selectWindow(TL_title); 
	getDimensions(width, height, channels, slices, frames);
	TL_width = width; 
	TL_height = height; 
	getPixelSize(unit, pixelWidth, pixelHeight);
	TL_pixelWidth = pixelWidth;
	TL_pixelHeight = pixelHeight; 
	
	selectWindow(FL_title); 
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
		selectWindow(TL_title); 
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
		selectWindow(TL_title); 
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
		selectWindow(FL_title); 
		makeRectangle(rectangle_x, rectangle_y, rectangle_width , rectangle_height);run("Duplicate...", " ");
		rename("temp_title"); 
	}
	selectWindow("temp_title"); 
	run("Size...", "width=" + SIM_width + " height=" + SIM_height + " depth=1 constrain average interpolation=None");
	saveAs("Tiff", file_path + File.separator + TL_title_no_ext + "_scaled.tif");
	TL_scaled_title = getTitle; 
	return TL_scaled_title; 
}


function segment_TL_image(TL_scaled_title, classifier_file){
	selectWindow(FL_title); 
	getDimensions(width, height, channels, slices, frames);
	SIM_width = width; 
	SIM_height = height; 
	
	selectWindow(TL_title); 
	//classifier_file = "Y:/private/Marie/Image_Analysis/2023-06-14-LIU-Mengru-count-spots-in-bacteria/input/Labkit/bacteria_TL_segmentation.classifier";
	run("Segment Image With Labkit", "segmenter_file=" + classifier_file + " use_gpu=false");
	saveAs("Tiff", file_path + File.separator + TL_title_no_ext + "_scaled_mask.tif");
	run("Size...", "width=" + SIM_width + " height=" + SIM_height + " depth=1 constrain average interpolation=None");
	TL_scaled_mask = getTitle();
	return TL_scaled_mask; 
	}

function determine_and_apply_xy_offset(TL_scaled_mask, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum, file_path, TL_title_no_ext, ROI_set_title, ROI_set_combined_title){
	selectWindow(TL_scaled_mask); 
	setAutoThreshold("Default dark");
	//run("Threshold...");
	run("Convert to Mask"); 
	run("Analyze Particles...", "size=" + bacteria_size_minimum + "-" + bacteria_size_maximum + " circularity=" + bacteria_circularity_minimum + "-" + bacteria_circularity_maximum + " exclude clear add");

	selectWindow(FL_title); 
	run("Enhance Contrast", "saturated=0.35");
	roiManager("Show All");

	ROI_to_move = Math.ceil(roiManager("count")/2); 
	roiManager("Select", ROI_to_move);
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
	saveAs("Tiff", file_path + File.separator + TL_title_no_ext + "_scaled.tif");
	close(); 
	selectWindow(TL_scaled_mask); 
	close(); 
	roiManager("reset");
	run("Select None");	
	return ROI_set_title; 
}


// ##################### SPOT COUNTING ##################### 
function spot_counting(FL_title, file_path, file_base_title, Ch1_prominence, Ch2_prominence){
	selectWindow(FL_title); 
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
		Ch2_kurtosis[i] = getResult("Kurt", i + number_of_bacteria);
		if (Ch1_kurtosis[i] > kurtosis_cutoff) {
			Ch1_class[i] = "spots";
		} 
		else{
			Ch1_class[i] = "homogeneous";
		}
		Ch1_mean[i] = getResult("Mean", i);
		Ch2_mean[i] = getResult("Mean", i + number_of_bacteria);
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
	//loopStop = 5; 
	loopStop = roiManager("count");
	
	// loop to count the number of spots per bacterium 
	
	selectWindow(FL_title); 
	setSlice(1); 
	//run("Brightness/Contrast...");
	//run("Enhance Contrast", "saturated=0.01");
	for(i=0; i<loopStop; i++) {
		roiManager("select", i);
		roiManager("rename", "Object " + i + 1);
		run("Find Maxima...", "noise="+Ch1_prominence+" output=[Count]");
		Ch1_spot_count[i] = getResult("Count", i);
		run("Find Maxima...", "noise="+Ch1_prominence+" output=[Point Selection]");
		run("Add Selection...");
	}

	selectWindow(FL_title); 
	roiManager("Show None");
	run("Select None");
	run("Duplicate...", "ignore"); 
	saveAs("Tiff", file_path + File.separator + file_base_title + "_Ch1-w-spots.tif");
	close();
	selectWindow("Results"); 
	run("Close");
	
	run("Remove Overlay");
	selectWindow(FL_title); 
	setSlice(2);  
	//run("Brightness/Contrast...");
	run("Enhance Contrast", "saturated=0.01");
	for(i=0; i<loopStop; i++) {
		roiManager("select", i);
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
	selectWindow(TL_title); 
	close(); 
	roiManager("reset");

	assemble_and_save_results(number_of_bacteria, Ch1_spot_count, Ch1_mean, Ch1_kurtosis, Ch1_class, Ch2_spot_count, Ch2_mean, Ch2_kurtosis, Ch2_class); 
}

function assemble_and_save_results(number_of_bacteria, Ch1_spot_count, Ch1_mean, Ch1_kurtosis, Ch1_class, Ch2_spot_count, Ch2_mean, Ch2_kurtosis, Ch2_class){
	run("Clear Results");
	//write results into a new results window 
	for (i = 0; i < number_of_bacteria; i++) {
		setResult("Bacterium ID", i, i+1);
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

function spot_heatmap(FL_title, file_path, ROI_set_combined_title, file_base_title, Ch1_prominence, Ch2_prominence, ROI_set_title, kurtosis_cutoff){
	//generate mask of bacteria bodies
	selectWindow(FL_title); 
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
	Dialog.addCheckbox("Intensity sum", false);
	Dialog.show();
	
	spot_peak_location_heatmap_choice = Dialog.getCheckbox();
	spot_size_heatmap_choice = Dialog.getCheckbox();
	intensity_sum_heatmap_choice = Dialog.getCheckbox();
	spots_image_title = "spots-locations";
	spots_file_title = file_base_title + "_spots_locations";
	newImage(spots_file_title, "8-bit grayscale-mode", FL_width, FL_height, FL_channels, 1, 1);
	
	// ##################### HEATMAP PRE-PROCESSING ##################### 
	
	//set parameters for heatmap and create empty first heatmap image
	aspect_ratio = 3.0;  
	heatmap_height = 26; 
	heatmap_width = Math.ceil(heatmap_height * aspect_ratio);
	roiManager("open", file_path + File.separator + ROI_set_title);
	ROI_initial_count = roiManager("count");
	
	for (i = 0; i < channels; i++) {
		channel_prominence_variable_name = "Ch" + (i+1) + "_prominence";
		//print(channel_prominence_variable_name); 		
		
		if (spot_peak_location_heatmap_choice == true){
			roiManager("reset");			
			print("Channel" + (i+1) + " - Spot location heatmap generation");
			points_ROI_file_title = "Channel" + (i+1) + "_points.roi";
			selectWindow(FL_title);
			Stack.setChannel(i+1);
			run("Select None");			
			run("Hide Overlay");		
			run("Find Maxima...", "prominence=" + channel_prominence_variable_name + " output=[Point Selection]");
			roiManager("Add");
			roiManager("select", 0);
			roiManager("save", file_path + File.separator + points_ROI_file_title);
			selectWindow(spots_file_title);
			Stack.setChannel(i+1);
			roiManager("select", 0);
			run("Draw", "slice");
			run("Subtract...", "value=254");
			setMinAndMax(0, 1);
			selectWindow(FL_title);
			roiManager("reset");
			run("Select None");
			selectWindow(spots_file_title);
			//run("Hide Overlay");
			saveAs("Tiff", file_path + File.separator + spots_file_title + ".tif");
			rename(spots_file_title); 
			roiManager("open", file_path + File.separator + ROI_set_title);
			old_heatmap_spot_peak_location = file_base_title + "_Ch_" + (i+1) + "_spot_peak_location_HeatMap_0";
			newImage(old_heatmap_spot_peak_location, "32-bit black", heatmap_width, heatmap_height, 1);
			run("Green Fire Blue");
			//rename(old_heatmap_spot_peak_location);
			old_heatmap_spot_peak_location = file_base_title + "_Ch_" + (i+1) + "_spot_peak_location_HeatMap_";
			heatmap_creation_loop(ROI_initial_count, spots_file_title, old_heatmap_spot_peak_location, bacteria_mask_file_title, bacteria_size_minimum, file_path);
			run("Select None");
		}
		if (spot_size_heatmap_choice == true){
			print("Channel" + (i+1) + " - Spot size heatmap generation");
			maxima_in_tolerance_file_title = "Channel" + (i+1) + "_maxima_in_tolerance";
			selectWindow(FL_title);
			roiManager("deselect");
			run("Select None");
			Stack.setChannel(i+1);			
			run("Find Maxima...", "prominence=" + channel_prominence_variable_name + " output=[Maxima Within Tolerance]");
			run("Subtract...", "value=254");
			setMinAndMax(0, 1);
			rename(maxima_in_tolerance_file_title);		
			old_heatmap_spot_size = file_base_title + "_Ch_" + (i+1) + "_spot_size_HeatMap_0";
			newImage(old_heatmap_spot_size, "32-bit black", heatmap_width, heatmap_height, 1);
			run("Green Fire Blue");
			//rename(old_heatmap_spot_size);
			old_heatmap_spot_size = file_base_title + "_Ch_" + (i+1) + "_spot_size_HeatMap_";
			heatmap_creation_loop(ROI_initial_count, maxima_in_tolerance_file_title, old_heatmap_spot_size, bacteria_mask_file_title, bacteria_size_minimum, file_path);
			roiManager("deselect");
		}
		
		if (intensity_sum_heatmap_choice == true){
			print("Channel" + (i+1) + " - Intensity heatmap generation");
			selectWindow(FL_title);
			Stack.setChannel(i+1);
			old_heatmap_intensity_sum = file_base_title + "_Ch_" + (i+1) + "_intensity_sum_location_HeatMap_0";
			newImage(old_heatmap_intensity_sum, "32-bit black", heatmap_width, heatmap_height, 1);
			run("Green Fire Blue");
			//rename(old_heatmap_intensity_sum);
			old_heatmap_intensity_sum = file_base_title + "_Ch_" + (i+1) + "_intensity_sum_location_HeatMap_";
			heatmap_creation_loop(ROI_initial_count, FL_title, old_heatmap_intensity_sum, bacteria_mask_file_title, bacteria_size_minimum, file_path);
			roiManager("deselect");
			//rescale heatmap to between 0 and 1
			getMinAndMax(min, max);
			//print(min); 
			//print(max); 		
			run("Subtract...", "value=" + min);
			run("Divide...", "value=" + (max - min));
			resetMinAndMax;
		}		
	}
}

function heatmap_creation_loop(ROI_initial_count, input_data_image, old_heatmap_file_name, bacteria_mask_file_title, bacteria_size_minimum, file_path){
		counter = 0; 
		//loop over all ROIs
		//for (k = 0; k < 3; k++) {
		for (k = 0; k < ROI_initial_count; k++) {
			
			ROI_of_choice = k;
			print("Processing object " + (k+1) +"/" + ROI_initial_count + ", Channel " + (i+1)); 
			//print("input data image: " + input_data_image); 
			//print("old heat map: " + old_heatmap_file_name); 
		
			//duplicate signal of input channel 
			selectWindow(FL_title);
			roiManager("select", ROI_of_choice);
			object_angle = getValue("Angle");
			object_kurtosis = getValue("Kurt");
			//print("Object angle: " + object_angle);
			//print("Object kurtosis: " + object_kurtosis);			

			selectWindow(input_data_image);
			roiManager("select", ROI_of_choice);		
			run("Duplicate...", " duplicate channels=" + (i+1) + " title=Spots_ROI_CH" + (i+1) + "_object_" + k);		
			run("Clear Outside");
			//rotate duplicated bacterial spots image by fit ellipse angle
			run("Rotate... ", "angle=" + object_angle + " grid=1 interpolation=Bicubic enlarge"); 
			rotated_object = getTitle();		
			
			selectWindow(bacteria_mask_file_title); 
			roiManager("select", ROI_of_choice);
			run("Duplicate...", "title=Mask_ROI_" + k);
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
		
			//scale image to certain aspect ratio --> check with Mengru which aspect ratio is desired
			run("Size...", "width=" + heatmap_width + " height=" + heatmap_height + " depth=1 average interpolation=Bicubic");
			resized_rotated_spots_title = getTitle();
			
			if (object_kurtosis < kurtosis_cutoff) {
				run("Set...", "value=0");
				counter = counter + 1;
				print("Rejected bacteria due to kurtosis: " + counter);  
			}
			
			old_heatmap = old_heatmap_file_name + k;
			//print("Old heat map: " + old_heatmap); 
			imageCalculator("Add 32-bit", old_heatmap, resized_rotated_spots_title);
			new_heatmap = old_heatmap_file_name + (k+1);
			//print("New heat map: " + new_heatmap); 
			rename(new_heatmap); 
			resetMinAndMax;
			
			//saveAs("TIFF", file_path + File.separator + new_heatmap + ".tif");
			//rename(new_heatmap); 
			
			selectWindow(old_heatmap); 
			close(); 			
			selectWindow(rotated_object);
			close();
			selectWindow(rotated_single_object_mask);
			close();
			selectWindow(resized_rotated_spots_title); 
			close(); 
			
			roiManager("deselect");
			
		}
		resetMinAndMax();
		saveAs("TIFF", file_path + File.separator + new_heatmap + ".tif");
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



//close("*");