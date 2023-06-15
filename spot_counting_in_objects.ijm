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








offset_TL_Fluo_x = -10; 
offset_TL_Fluo_y = 0; 
SIM_sampling = 2; 

offset_TL_Fluo_x = offset_TL_Fluo_x / SIM_sampling; 
offset_TL_Fluo_y = offset_TL_Fluo_y / SIM_sampling; 
























roiManager("reset");
open("Z:/private/Marie/Image Analysis/2023-06-14-LIU-Mengru-count-spots-in-bacteria/input/2-2.czi");
originalTitle = getTitle();

run("Translate...", "x=" + offset_TL_Fluo_x +" y=" + offset_TL_Fluo_y +" interpolation=None");
//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");



run("Subtract Background...", "rolling=10 light");
selectWindow(originalTitle);
run("Median...", "radius=2");
run("Duplicate...", " ");
duplicateTitle = getTitle();

selectWindow(originalTitle);
//run("Threshold...");
setAutoThreshold("Otsu no-reset");
run("Convert to Mask");

run("Size...", "width=2560 height=2560 depth=1 constrain average interpolation=None");
run("Duplicate...", " ");




selectWindow(duplicateTitle);
run("Find Maxima...", "prominence=100 light output=[Segmented Particles]");
segmentedParticlesTitle = getTitle();

run("Size...", "width=2560 height=2560 depth=1 constrain average interpolation=None");



run("Paste Control...");
setPasteMode("AND");
selectWindow(segmentedParticlesTitle); 
run("Copy");
selectWindow(originalTitle);
run("Paste");
run("Select None");

run("Analyze Particles...", "size=0.60-3.00 circularity=0.00-0.8 exclude clear add");





//Parameters to adjust if required
foci_maxima_segmentation_smoothing = 3;
background_subtration_rolling_ball_radius = 50;
prominence = 1000;



// -- foci segmentation start
open("Z:/private/Marie/Image Analysis/2023-06-14-LIU-Mengru-count-spots-in-bacteria/input/2-1.czi");
FL_title = getTitle();

//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");
number_of_bacteria = roiManager("count") + 1; 

Ch1_spot_count = newArray(number_of_bacteria);
Ch2_spot_count = newArray(number_of_bacteria);
loopStop = 20; 

// loop to count the number of foci per cell 

selectWindow(FL_title); 
setSlice(1); 
for(i=0; i<loopStop; i++) {
//for(i=0; i<roiManager("count"); i++) {
	roiManager("select", i);
	roiManager("rename", "Bacterium " + i + 1);
	run("Find Maxima...", "noise="+prominence+" output=[Count]");
	Ch1_spot_count[i] = getResult("Count", i);
	run("Find Maxima...", "noise="+prominence+" output=[Point Selection]");
	run("Add Selection...");
}
roiManager("Show None");
selectWindow("Results"); 
run("Close");
saveAs("Tiff", "Z:/private/Marie/Image Analysis/2023-06-14-LIU-Mengru-count-spots-in-bacteria/Playground/2-spots.tif");
results_image_title = getTitle();


run("Remove Overlay");
selectWindow(results_image_title); 
setSlice(2);  
for(i=0; i<loopStop; i++) {
//for(i=0; i<roiManager("count"); i++) {
	roiManager("select", i);
	run("Find Maxima...", "noise="+prominence+" output=[Count]");
	Ch2_spot_count[i] = getResult("Count", i);
	run("Find Maxima...", "noise="+prominence+" output=[Point Selection]");
	run("Add Selection...");
}
roiManager("Show None");
selectWindow("Results"); 
run("Close");
saveAs("Tiff", "Z:/private/Marie/Image Analysis/2023-06-14-LIU-Mengru-count-spots-in-bacteria/Playground/2-spots.tif");
// -- foci segmentation end


	
//write results into a new results window 
for (i = 0; i < number_of_bacteria; i++) {
	setResult("Bacterium ID", i, i+1);
	setResult("Ch1 spot count", i, Ch1_spot_count[i]);
	setResult("Ch2 spot count", i, Ch2_spot_count[i]);  
}

saveAs("Results", "Z:/private/Marie/Image Analysis/2023-06-14-LIU-Mengru-count-spots-in-bacteria/spot-counts.csv");
selectWindow("Results"); 
close(); 
close("*");