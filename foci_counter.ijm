/*


												- Written by Marie Held [mheldb@liverpool.ac.uk] March 2024
												  Liverpool CCI (https://cci.liverpool.ac.uk/)
________________________________________________________________________________________________________________________

BSD 2-Clause License

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*
*/


//Parameters to adjust if required
min_cell_size = 350;  								// default: 500; increase to reject more smaller cells, 
cell_segmentation_smoothing = 6; 					// default: 6; decrease if cell features are lost or artifacts appear
foci_maxima_segmentation_smoothing = 3; 			// default: 3; decrease if foci are lost
background_subtration_rolling_ball_radius = 10; 	// default: 10; increase if cell background is not removed sufficiently
prominence = 15000; 								// default 10000; decrease if not all foci are identified
line_width = 2; 									// default: 2; do not change

// ---- Automated script start ----
getPixelSize(unit, pixelWidth, pixelHeight);
getDimensions(width, height, channels, slices, frames);


//scaled parameters/
scaling_factor = pixelWidth / 0.207567;
scaled_image_width = Math.ceil(scaling_factor * width); 
scaled_image_height = Math.ceil(scaling_factor * height); 
run("Size...", "width=" + scaled_image_width + " height=" + scaled_image_height + " depth=" + channels + " constrain average interpolation=Bilinear");

// -- cell segmentation start --
title = getTitle();
selectWindow(title);
setSlice(1); 
run("Duplicate...", " ");
run("Median...", "radius="  + cell_segmentation_smoothing); // pre-processing - smoothing
setAutoThreshold("Huang dark"); // classifying pixels into two categories: background and object (cell) pixels
//run("Threshold...");
run("Convert to Mask");
run("Fill Holes");
//run("Watershed");
mask_title = getTitle();
setForegroundColor(0, 0, 0);
setBackgroundColor(255, 255, 255);
setTool("freeline");
run("Line Width...", "line="+line_width);
//ask user to check segmentation and add cell borders where necessary
waitForUser("Manually separate touching objects by drawing a freehand line and then pressing [d] to overlay the line. ");
run("Convert to Mask");
roiManager("Reset");
run("Analyze Particles...", "size=" + min_cell_size + "-Infinity clear add"); // assign a unique identifier to each object and create outline
// -- cell segmentation end --

// -- foci segmentation start --
selectWindow(title);
run("Duplicate...", " ");
foci_duplicate = getTitle(); 
run("Median...", "radius=foci_maxima_segmentation_smoothing");
run("Subtract Background...", "rolling=" + background_subtration_rolling_ball_radius);
//run("Median...", "radius=foci_maxima_segmentation_smoothing");
//run("Brightness/Contrast...");
resetMinAndMax();


// loop to count the number of foci per cell 
for(i=0; i<roiManager("count"); i++) {
	selectWindow(foci_duplicate);
	roiManager("select", i);
	roiManager("rename", "Cell " + i + 1);
	run("Find Maxima...", "noise="+prominence+" output=[Count]");
	run("Find Maxima...", "noise="+prominence+" output=[Point Selection]");
	run("Add Selection...");
}
// -- foci segmentation end --

run("Line Width...", "line=1");
setTool("line");

//notification that the script has finished
print(title + " done!")
beep();



TEST