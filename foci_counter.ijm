//Parameters to adjust if required
min_cell_size = 500;
cell_segmentation_smoothing = 6; 
foci_maxima_segmentation_smoothing = 3;
background_subtration_rolling_ball_radius = 50;
prominence = 50;
line_width = 2; 

// ---- Automated script start

// -- cell segmentation start
title = getTitle();
selectWindow(title);
run("Duplicate...", " ");
run("Median...", "radius="  + cell_segmentation_smoothing); // pre-processing - smoothing
setAutoThreshold("Huang dark"); // classifying pixels into two categories: background and object
//run("Threshold...");
run("Convert to Mask");
//run("Watershed");
mask_title = getTitle();
setTool("freeline");
run("Line Width...", "line="+line_width);
waitForUser("Manually separate touching objects by drawing a freehand line and then pressing [d] to overlay the line. ");
run("Convert to Mask");

roiManager("Reset");
run("Analyze Particles...", "size=" + min_cell_size + "-Infinity clear add"); // assign a unique identifier to each object and create outline
// -- cell segmentation end

// -- foci segmentation start
selectWindow(title);
run("Duplicate...", " ");
run("Median...", "radius=foci_maxima_segmentation_smoothing");
run("Subtract Background...", "rolling=" + background_subtration_rolling_ball_radius);
//run("Brightness/Contrast...");
resetMinAndMax();
run("8-bit");

// loop to count the number of foci per cell 
for(i=0; i<roiManager("count"); i++) {
	roiManager("select", i);
	roiManager("rename", "Cell " + i + 1);
	run("Find Maxima...", "noise="+prominence+" output=[Count]");
	run("Find Maxima...", "noise="+prominence+" output=[Point Selection]");
	run("Add Selection...");
}
// -- foci segmentation end

//notification that the script has finished
print(title + " done!")
beep();
run("Line Width...", "line="+line_width);