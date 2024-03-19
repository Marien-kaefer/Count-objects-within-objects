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
scaled_image_width = scaling_factor * width; 
scaled_image_height = scaling_factor * height; 
run("Size...", "width=" + scaled_image_width + " height=" + scaled_image_height + " depth=" + channels + " constrain average interpolation=Bilinear");

// -- cell segmentation start --
title = getTitle();
selectWindow(title);
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