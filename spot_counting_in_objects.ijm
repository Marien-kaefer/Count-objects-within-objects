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

*/


#@ File (label = "Input directory", style = "open") inputFile
#@ Integer(label="SIM sampling rate applied", value=2, min=2, max=4, style="slider") SIM_sampling
#@ String(value="You will be required to choose the Transmitted light file of the file set you wish to process. Please choose the following parameters.", visibility="MESSAGE") message
#@ Integer(label="TL background subtraction radius: ", value = 10) BG_subtr_radius
#@ Integer(label="TL median filter radius: ", value = 2) median_radius
#@ String (choices={"Otsu","Default", "Huang","Intermodes","IsoData","IJ_IsoData","Li","MaxEntropy","Mean","MinError","Minimum","Moments","Percentile","RenyiEntropy","Shanbhag","Triangle","Yen"}, style="listBox") threshold_algorithm
#@ Integer(label="TL prominence: " , value = 100) TL_prominence
#@ Integer(label="TL image offset in x: " , value = -2) offset_TL_Fluo_x
#@ Integer(label="TL image offset in y: " , value = -1) offset_TL_Fluo_y
#@ Integer(label="Ch1 Prominence: " , value = 1000) Ch1_prominence
#@ Integer(label="Ch2_prominence: " , value = 1000) Ch2_prominence
#@ Double(label="Bacteria size minimum: ", value = 0.6) bacteria_size_minimum
#@ Double(label="Bacteria size maximum: ", value = 3.0) bacteria_size_maximum
#@ Double(label="Bacteria circularity minimum: ", value = 0.0) bacteria_circularity_minimum
#@ Double(label="Bacteria circularity maximum: ", value = 0.8) bacteria_circularity_maximum










run("Bio-Formats Windowless Importer", "open=[" + inputFile + "]");
file_path = File.getDirectory(inputFile);
TL_title = File.getName(inputFile);
originalTitle = TL_title; 

TL_title = file_name_remove_extension(TL_title); 
selectWindow(TL_title + ".czi"); 
rename(TL_title); 
file_base_title = substring(TL_title, 0, (lengthOf(TL_title)-1))
FL_title = file_base_title + "1"; 


//offset_TL_Fluo_x = -2; 
//offset_TL_Fluo_y = -1; 


offset_TL_Fluo_x = offset_TL_Fluo_x / SIM_sampling; 
offset_TL_Fluo_y = offset_TL_Fluo_y / SIM_sampling; 

roiManager("reset");
run("Clear Results");
run("Bio-Formats Windowless Importer", "open=[" + inputFile + "]");
originalTitle = getTitle();
getDimensions(width, height, channels, slices, frames);

run("Translate...", "x=" + offset_TL_Fluo_x +" y=" + offset_TL_Fluo_y +" interpolation=None");
//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");

run("Subtract Background...", "rolling=" + BG_subtr_radius + " light");
selectWindow(originalTitle);
run("Median...", "radius=" + median_radius);
run("Duplicate...", " ");
duplicateTitle = getTitle();

selectWindow(originalTitle);
//run("Threshold...");
setAutoThreshold(threshold_algorithm + " no-reset");
run("Convert to Mask");


run("Size...", "width=" + (SIM_sampling * width) + " height=" + (SIM_sampling * width) + " depth=1 constrain average interpolation=None");

saveAs("Tiff", file_path + File.separator + file_base_title + "TL.tif");
rename(originalTitle); 
run("Duplicate...", " ");


selectWindow(duplicateTitle);
run("Find Maxima...", "prominence=" + TL_prominence + " light output=[Segmented Particles]"); 
segmentedParticlesTitle = getTitle();

//run("Size...", "width=2560 height=2560 depth=1 constrain average interpolation=None");//soft code!!!!!!!!!!!!!!!
run("Size...", "width=" + (SIM_sampling * width) + " height=" + (SIM_sampling * width) + " depth=1 constrain average interpolation=None");

run("Paste Control...");
setPasteMode("AND");
selectWindow(segmentedParticlesTitle); 
run("Copy");
selectWindow(originalTitle);
run("Paste");
run("Select None");

//run("Analyze Particles...", "size=0.60-3.00 circularity=0.00-0.8 exclude clear add");
run("Analyze Particles...", "size=" + bacteria_size_minimum + "-" + bacteria_size_maximum + " circularity=" + bacteria_circularity_minimum + "-" + bacteria_circularity_minimum + " exclude clear add");
roiManager("deselect");
roiManager("save", file_path + File.separator + file_base_title + "_bacteriaROIs.zip");




// -- foci segmentation start

run("Bio-Formats Windowless Importer", "open=[" + file_path + File.separator + FL_title + ".czi]");
rename(FL_title); 
//run("Subtract Background...", "rolling=5 stack");

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
	if (Ch1_kurtosis[i] > 0) {
		Ch1_class[i] = "spots";
	} 
	else{
		Ch1_class[i] = "homogeneous";
	}
	Ch1_mean[i] = getResult("Mean", i);
	Ch2_mean[i] = getResult("Mean", i + number_of_bacteria);
	if (Ch2_kurtosis[i] > 0) {
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
loopStop = 5; 

// loop to count the number of foci per cell 

selectWindow(FL_title); 
setSlice(1); 
//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.15");
for(i=0; i<loopStop; i++) {
//for(i=0; i<roiManager("count"); i++) {
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
saveAs("Tiff", file_path + File.separator + file_base_title + "Ch1-w-spots.tif");
close();
selectWindow("Results"); 
run("Close");

run("Remove Overlay");
selectWindow(FL_title); 
setSlice(2);  
//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.15");
for(i=0; i<loopStop; i++) {
//for(i=0; i<roiManager("count"); i++) {
	roiManager("select", i);
	run("Find Maxima...", "noise="+Ch2_prominence+" output=[Count]");
	Ch2_spot_count[i] = getResult("Count", i);
	run("Find Maxima...", "noise="+Ch2_prominence+" output=[Point Selection]");
	run("Add Selection...");
}
roiManager("Show None");
run("Select None");
run("Duplicate...", "ignore"); 
saveAs("Tiff", file_path + File.separator + file_base_title + "Ch2-w-spots.tif");
close(); 
run("Clear Results");
// -- foci segmentation end


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

function file_name_remove_extension(originalTitle){
	dotIndex = lastIndexOf(originalTitle, "." ); 
	file_name_without_extension = substring(originalTitle, 0, dotIndex );
	//print( "Name without extension: " + file_name_without_extension );
	return file_name_without_extension;
}