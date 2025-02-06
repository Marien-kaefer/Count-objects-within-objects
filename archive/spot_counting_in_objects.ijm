/*
Macro to segment bacteria using the transmitted light channels and count the number of puncta per bacterium in two fluorescence channels.

												- Written by Marie Held [mheldb@liverpool.ac.uk] June 2023
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
#@ Integer(label="SIM sampling rate applied during processing:", value=2, min=2, max=4, style="slider") SIM_sampling
#@ String(value="Please choose the following processing parameters.", visibility="MESSAGE") message2
#@ Integer(label="TL background subtraction radius: ", value = 10) BG_subtr_radius
#@ Integer(label="TL median filter radius: ", value = 2) median_radius
#@ String (choices={"Otsu","Default", "Huang","Intermodes","IsoData","IJ_IsoData","Li","MaxEntropy","Mean","MinError","Minimum","Moments","Percentile","RenyiEntropy","Shanbhag","Triangle","Yen"}, style="listBox") threshold_algorithm
#@ Integer(label="TL prominence: " , value = 100) TL_prominence
#@ String(value="To determine the offset between the transmitted light and fluorescence channels, choose 'Yes'.", visibility="MESSAGE") message3
#@ String(choices={"Yes", "No"}, style="radioButtonHorizontal") offset_determination_choice
#@ String(value="If offset between the transmitted light and fluorescence channels, is known, choose 'No' above and enter values below.", visibility="MESSAGE") message4
#@ Integer(label="TL image offset in x: " , value = 0) offset_TL_Fluo_x
#@ Integer(label="TL image offset in y: " , value = 0) offset_TL_Fluo_y
#@ Integer(label="Ch1 Prominence: " , value = 1000) Ch1_prominence
#@ Integer(label="Ch2_prominence: " , value = 1000) Ch2_prominence
#@ Double(label="Bacteria size minimum: ", value = 0.8) bacteria_size_minimum
#@ Double(label="Bacteria size maximum: ", value = 3.0) bacteria_size_maximum
#@ Double(label="Bacteria circularity minimum: ", value = 0.0) bacteria_circularity_minimum
#@ Double(label="Bacteria circularity maximum: ", value = 0.8) bacteria_circularity_maximum
#@ Double(label="Intensity histogram kurtosis cutoff value: ", value = 3.0) kurtosis_cutoff

pre_clean_up();
//variable definitions and calculations
file_path = File.getDirectory(inputFile);
TL_title = File.getName(inputFile);
originalTitle = TL_title; 
TL_title = file_name_remove_extension(TL_title); 
file_base_title = substring(TL_title, 0, (lengthOf(TL_title)-2))
FL_title = file_base_title + "-1"; 
offset_TL_Fluo = newArray(2); 
offset_TL_Fluo_x = offset_TL_Fluo[0] / SIM_sampling; 
offset_TL_Fluo_y = offset_TL_Fluo[1] / SIM_sampling; 

//offset determination choice conditional
if (offset_determination_choice == "Yes") {
	get_channel_offset(inputFile, offset_TL_Fluo_x, offset_TL_Fluo_y, offset_TL_Fluo, BG_subtr_radius, originalTitle, threshold_algorithm, SIM_sampling, file_path, file_base_title, TL_prominence, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum); 
	close("*");	
	offset_TL_Fluo_x = offset_TL_Fluo[0] / SIM_sampling; 
	offset_TL_Fluo_y = offset_TL_Fluo[1] / SIM_sampling; 
}	

start = getTime(); 
TL_preprocessing(inputFile, offset_TL_Fluo_x, offset_TL_Fluo_y, BG_subtr_radius, originalTitle, threshold_algorithm, SIM_sampling, file_path, file_base_title, TL_prominence, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum);

// -- masurements
run("Bio-Formats Windowless Importer", "open=[" + file_path + File.separator + FL_title + ".czi]");
rename(FL_title); 

run("Set Measurements...", "area mean shape kurtosis redirect=None decimal=3");
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

// -- foci segmentation start

Ch1_spot_count = newArray(number_of_bacteria);
Ch2_spot_count = newArray(number_of_bacteria);
//loopStop = 5; 
loopStop = roiManager("count");

// loop to count the number of foci per bacterium 

selectWindow(FL_title); 
setSlice(1); 
//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.01");
for(i=0; i<loopStop; i++) {
	roiManager("select", i);
	roiManager("rename", "Bacterium " + i + 1);
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
run("Duplicate...", "ignore"); 
saveAs("Tiff", file_path + File.separator + file_base_title + "_Ch2-w-spots.tif");
close(); 
run("Clear Results");
// -- foci segmentation end

write_input_parameters_to_file(file_path, file_base_title, BG_subtr_radius, median_radius, threshold_algorithm, TL_prominence, offset_TL_Fluo_x, offset_TL_Fluo_y, Ch1_prominence, Ch2_prominence, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum, number_of_bacteria); 
assemble_and_save_results(number_of_bacteria, Ch1_spot_count, Ch1_mean, Ch1_kurtosis, Ch1_class, Ch2_spot_count, Ch2_mean, Ch2_kurtosis, Ch2_class); 
clean_up(); 

//let user know the process has finished and how long it took
stop = getTime(); 
duration = stop - start;
duration_String = duration_conversion(duration);
print("The processing took " + duration_conversion(duration));
print(TimeStamp() + ": Processing of [" + originalTitle + "] complete.");
beep();



//----------------------------- FUNCTIONS -----------------------------//

function pre_clean_up(){
	roiManager("reset");
	run("Clear Results");
}

function file_name_remove_extension(originalTitle){
	dotIndex = lastIndexOf(originalTitle, "." ); 
	file_name_without_extension = substring(originalTitle, 0, dotIndex );
	//print( "Name without extension: " + file_name_without_extension );
	return file_name_without_extension;
}

function TL_preprocessing(inputFile, offset_TL_Fluo_x, offset_TL_Fluo_y, BG_subtr_radius, originalTitle, threshold_algorithm, SIM_sampling, file_path, file_base_title, TL_prominence, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum){
	run("Bio-Formats Windowless Importer", "open=[" + inputFile + "]");
	
	getDimensions(width, height, channels, slices, frames);
	
	run("Translate...", "x=" + offset_TL_Fluo_x +" y=" + offset_TL_Fluo_y +" interpolation=None");
	//run("Brightness/Contrast...");
	run("Enhance Contrast", "saturated=0.35");
	
	run("Subtract Background...", "rolling=" + BG_subtr_radius + " light");
	selectWindow(originalTitle);
	run("Median...", "radius=" + median_radius);
	run("Duplicate...", " ");
	duplicateTitle = getTitle();
	
	//create merged channel image
	run("Duplicate...", " ");
	duplicate_for_merge_Title = getTitle();
	run("Size...", "width=" + (SIM_sampling * width) + " height=" + (SIM_sampling * width) + " depth=1 constrain average interpolation=None");
	//run("Brightness/Contrast...");
	resetMinAndMax();
	run("Bio-Formats Windowless Importer", "open=[" + file_path + File.separator + FL_title + ".czi]");
	temp_title = getTitle();
	run("Split Channels");
	Ch1_title = "C1-" + temp_title;
	Ch2_title = "C2-" + temp_title;
	run("Merge Channels...", "c1=" + Ch1_title + " c2=" + Ch2_title + " c3=" + duplicate_for_merge_Title + " create");
	setSlice(1);
	run("Green");
	setSlice(2);
	run("Magenta");
	setSlice(3);
	run("Grays");
	saveAs("Tiff", file_path + File.separator + file_base_title + "_merge.tif");
	close(); 
	
	selectWindow(originalTitle);
	//run("Threshold...");
	setAutoThreshold(threshold_algorithm + " no-reset");
	run("Convert to Mask");
	
	run("Size...", "width=" + (SIM_sampling * width) + " height=" + (SIM_sampling * width) + " depth=1 constrain average interpolation=None");

	saveAs("Tiff", file_path + File.separator + file_base_title + "_TL-mask.tif");
	rename(originalTitle); 
	
	selectWindow(duplicateTitle);
	run("Find Maxima...", "prominence=" + TL_prominence + " light output=[Segmented Particles]"); 
	segmentedParticlesTitle = getTitle();

	run("Size...", "width=" + (SIM_sampling * width) + " height=" + (SIM_sampling * width) + " depth=1 constrain average interpolation=None");
	
	run("Paste Control...");
	setPasteMode("AND");
	selectWindow(segmentedParticlesTitle); 
	run("Copy");
	selectWindow(originalTitle);
	run("Paste");
	run("Select None");
	
	run("Analyze Particles...", "size=" + bacteria_size_minimum + "-" + bacteria_size_maximum + " circularity=" + bacteria_circularity_minimum + "-" + bacteria_circularity_maximum + " exclude clear add");
	roiManager("deselect");
	roiManager("save", file_path + File.separator + file_base_title + "_bacteriaROIs.zip");
}

//enable user to determine the offset between the different imaging modalities 
function get_channel_offset(inputFile, offset_TL_Fluo_x, offset_TL_Fluo_y, offset_TL_Fluo, BG_subtr_radius, originalTitle, threshold_algorithm, SIM_sampling, file_path, file_base_title, TL_prominence, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum){
	TL_preprocessing(inputFile, offset_TL_Fluo_x, offset_TL_Fluo_y, BG_subtr_radius, originalTitle, threshold_algorithm, SIM_sampling, file_path, file_base_title, TL_prominence, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum);
	run("Bio-Formats Windowless Importer", "open=[" + file_path + File.separator + FL_title + ".czi]");
	//run("Brightness/Contrast...");
	run("Enhance Contrast", "saturated=0.35");
	roiManager("Show All");

	waitForUser("Determine the offset between the TL and FL images by making a note of the number of pixels the ROI should be shifted to the right (+x), left (-x), down (+y) or up (-y). ");

	title = "Please enter the offsets in x and/or y. ";
	Dialog.create("Offset entry");
	Dialog.addNumber("x:", offset_TL_Fluo_x);
	Dialog.addNumber("y:", offset_TL_Fluo_y);
	Dialog.show();
	offset_TL_Fluo_x = Dialog.getNumber();
	offset_TL_Fluo_y = Dialog.getNumber();
	offset_TL_Fluo[0] = offset_TL_Fluo_x; 
	offset_TL_Fluo[1] = offset_TL_Fluo_y;
	return offset_TL_Fluo;
}

//write all the used input parameters into a text file
function write_input_parameters_to_file(file_path, file_base_title, BG_sutr_radius, median_radius, threshold_algorithm, TL_prominence, offset_TL_Fluo_x, offset_TL_Fluo_y, Ch1_prominence, Ch2_prominence, bacteria_size_minimum, bacteria_size_maximum, bacteria_circularity_minimum, bacteria_circularity_maximum, number_of_bacteria){ 
	parameters_output_file = File.open(file_path + File.separator + file_base_title + "_analysis_parameters.txt"); 
	print(parameters_output_file, "TL background subtraction radius: " + BG_subtr_radius);
    print(parameters_output_file, "TL median filter radius: " + median_radius);
    print(parameters_output_file, "Threshold algorithm: " + threshold_algorithm);
    print(parameters_output_file, "TL prominence: " + TL_prominence);
    print(parameters_output_file, "TL image offset in x: " + (offset_TL_Fluo_x * SIM_sampling));
    print(parameters_output_file, "TL image offset in y: " + (offset_TL_Fluo_y * SIM_sampling));
    print(parameters_output_file, "Ch1 Prominence: " + Ch1_prominence);
    print(parameters_output_file, "Ch2 Prominence: " + Ch2_prominence);
    print(parameters_output_file, "Bacteria size minimum: " + bacteria_size_minimum);
    print(parameters_output_file, "Bacteria size maximum: " + bacteria_size_maximum);
    print(parameters_output_file, "Bacteria circularity minimum: " + bacteria_circularity_minimum);
    print(parameters_output_file, "Bacteria circularity maximum: " + bacteria_circularity_maximum);
	print(parameters_output_file, "Number of objects measured: " + number_of_bacteria);
	File.close(parameters_output_file)
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
	saveAs("Results", file_path + File.separator + file_base_title + "_results.csv");
	close("*");
}

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

//clean up: close results window, close all open image windows
function clean_up(){
	close("*");
	run("Clear Results");
}











