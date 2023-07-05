FileList = getList("image.titles");
//Array.print(FileList);


title = "Files and Parameters";
Dialog.create("Calibration Dialog");
Dialog.addChoice("Transmitted light file", FileList);
Dialog.addChoice("Fluorescence file", FileList);
Dialog.addNumber("Transmitted light prominence", 500);
Dialog.addNumber("SIM prominence", 10000);
Dialog.addNumber("Bacteria size minimum: ", 0.3)
Dialog.addNumber("Bacteria size maximum: ", 2.0) 
Dialog.show();

phaseStack = Dialog.getChoice();
fluorescenceStack = Dialog.getChoice();
TL_prominence = Dialog.getNumber();
Ch1_prominence = Dialog.getNumber();
bacteria_size_minimum = Dialog.getNumber();
bacteria_size_maximum = Dialog.getNumber();

selectWindow(phaseStack); 
getDimensions(width, height, channels, slices, frames);
TL_width = width; 
TL_height = height; 
getPixelSize(unit, pixelWidth, pixelHeight);
TL_pixelWidth = pixelWidth;
TL_pixelHeight = pixelHeight; 

TL_original_title = getTitle();
file_path = getDirectory("image"); 

file_base_title = substring(TL_original_title, 0, indexOf(TL_original_title, "_"));


selectWindow(fluorescenceStack); 
FL_original_title = getTitle();
getDimensions(width, height, channels, slices, frames);
SIM_width = width; 
SIM_height = height; 
getPixelSize(unit, pixelWidth, pixelHeight);
SIM_pixelWidth = pixelWidth;
SIM_pixelHeight = pixelHeight; 

SIM_to_TL = SIM_width / TL_width; 
scale_factor = (TL_pixelWidth / SIM_pixelWidth) / SIM_to_TL ;
//print("Scale factor " + scale_factor); 
rectangle_width = floor(TL_width / scale_factor);
//print("Rectangle width " + rectangle_width); 
rectangle_height = floor(TL_height / scale_factor);
//print("Rectangle height " + rectangle_height); 
rectangle_x = (TL_width - rectangle_width)/2;
//print("Rectangle x " + rectangle_x); 
rectangle_y = (TL_height - rectangle_height)/2;
//print("Rectangle x " + rectangle_x); 

selectWindow(TL_original_title); 
makeRectangle(rectangle_x, rectangle_y, rectangle_width , rectangle_height);
run("Size...", "width=" + SIM_width + " height=" + SIM_height + " depth=1 constrain average interpolation=None");
run("Duplicate...", " ");
TL_crop_title = getTitle();





run("Median...", "radius=2");
run("Duplicate...", " ");
TL_threshold_title = getTitle();
selectWindow(TL_original_title);
run("Duplicate...", " ");
TL_segmented_particles_title = getTitle();


selectWindow(TL_threshold_title);
setAutoThreshold("Otsu");
//run("Threshold...");
run("Convert to Mask");
selectWindow(TL_segmented_particles_title);
run("Find Maxima...", "prominence=" + TL_prominence +" light output=[Segmented Particles]");
rename(TL_segmented_particles_title); 
run("Paste Control...");
setPasteMode("AND");
run("Copy");
selectWindow(TL_threshold_title);
run("Paste");
run("Select None");

run("Analyze Particles...", "size=" + bacteria_size_minimum + "-" + bacteria_size_maximum + " exclude clear add");



// determine required offset between TL and SIM images 
selectWindow(FL_original_title); 
roiManager("Show All");

ROI_to_move = roiManager("count")/2; 
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

roiManager("save", file_path + File.separator + file_base_title + "_bacteriaROIs.zip");
//measure parameters
run("Set Measurements...", "area mean standard modal min shape integrated median skewness kurtosis redirect=None decimal=3");
roiManager("Deselect");
roiManager("multi-measure");
saveAs("Results", file_path + File.separator + file_base_title + "_results.csv");
run("Clear Results");

prominence = 10000; 
Ch1_spot_count = newArray(number_of_ROI); 
//count spots in bacteria
for(i=0; i<number_of_ROI; i++) {
	roiManager("select", i);
	run("Find Maxima...", "noise="+Ch1_prominence+" output=[Count]");
	Ch1_spot_count[i] = getResult("Count", i);
	run("Find Maxima...", "noise="+Ch1_prominence+" output=[Point Selection]");
	run("Add Selection...");
}

saveAs("Results", file_path + File.separator + file_base_title + "_spot-counts.csv");
close("Results");

selectWindow(fluorescenceStack); 
saveAs("Tiff", file_path + File.separator + file_base_title + "_w-spots.tif");

roiManager("reset");
close("*");
