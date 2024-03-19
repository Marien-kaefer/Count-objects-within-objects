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



//get input parameters
#@ String(value="Please select the file you wish to process and set the relevant parameters.", visibility="MESSAGE") message
#@ Double(label="StarDist probability/score threshold: ", value = 0.35) probability_threshold
#@ Double(label="StarDist overlap threshold: ", value = 0.6) overlap_threshold
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix of SUM projection file", value = ".tif", persist=false) suffix


run("Fresh Start");
processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i])) //also process files in directories within the input directory
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix)){
			run("Fresh Start"); //clean up, closing open images, reset ROI manager, clears results window
			inputFile = input + File.separator + list[i]; 
			open(inputFile); 
			file_name_without_extension = File.nameWithoutExtension; 
			rename(file_name_without_extension); 
			//print(file_name_without_extension); 
			//directory = File.getParent(inputFile); 
			//print("Directory: " + directory); 
			
			convert_instance_mask_to_ROIs(file_name_without_extension, input, output);
			label_image = segment_bacteria_with_StarDist(file_name_without_extension, output, probability_threshold, overlap_threshold);
			number_of_bacteria = count_instance_mask_objects_in_ROI(output, file_name_without_extension, label_image); 
			assemble_and_save_results(number_of_bacteria, output); 
			 
	}
	run("Fresh Start"); //clean up, closing open images, reset ROI manager, clears results window
}

waitForUser("Done!"); 

function convert_instance_mask_to_ROIs(file_name_without_extension, input, output){
	cells_label_file_name = file_name_without_extension + "_cp_masks.png"; //open Cellpose cell instance mask
	open(input + File.separator + cells_label_file_name);
	run("glasbey_on_dark");
	run("LabelMap to ROI Manager (2D)"); // generate ROIs from instance mask
	roiManager("Save", output + File.separator + file_name_without_extension + "_cells.zip"); //save ROI set
}

function segment_bacteria_with_StarDist(file_name_without_extension, output, probability_threshold, overlap_threshold){
	selectWindow(file_name_without_extension); 
	run("Duplicate...", "duplicate channels=3"); // extract channel 3, assumed to be the bacteria channel
	bacteria_image = getTitle(); 
	
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'" + bacteria_image + "', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'100.0', 'probThresh':'" + probability_threshold + "', 'nmsThresh':'" + overlap_threshold + "', 'outputType':'Label Image', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'true', 'showCsbdeepProgress':'true', 'showProbAndDist':'false'], process=[false]"); //identify bacteria using StarDist plugin
	selectImage("Label Image");
	saveAs("Tiff", output + File.separator + file_name_without_extension + "_bact_labelmask.tif");
	label_image = getTitle();
	return label_image; 
}

function count_instance_mask_objects_in_ROI(output, file_name_without_extension, label_image){
	nBins = 256; 
	number_of_bacteria = newArray(roiManager("count")); 
	
	for (i = 0; i < roiManager("count") ; i++) {
		//print("i: " + i);
		values = 0; 
		counts = 0; 
		selectWindow(label_image); 
		roiManager("Select", i);
		roiManager("Rename", i);
		roiManager("Save", output + File.separator + file_name_without_extension + "_cells.zip");
		getHistogram(values, counts, nBins); // get histogram of "gray values" of instance mask. Each bin with >0 frequency is one cell plus one for background
		IDs = counts; 
		IDs_no_zero = Array.deleteValue(IDs, 0); // remove all bins that have 0 frequency
		//Array.show("title", IDs, IDs_no_zero);
		
		number_of_bacteria[i] = IDs_no_zero.length - 1; // measure length of array and subtract one (background) to get number of bacteria within the cell
		//Array.print(number_of_bacteria); 
	}
	//Array.print(number_of_bacteria); 
	return number_of_bacteria; 
}

function assemble_and_save_results(number_of_bacteria, output){
	run("Clear Results");
	//write results into a new results window 
	for (i = 0; i < roiManager("count"); i++) {
		setResult("Cell ID", i, i+1);
		setResult("Bacteria Count", i, number_of_bacteria[i]);
	}
	saveAs("Results", output + File.separator + file_name_without_extension + "_bacteria_count_results.csv");
	close("Results");
}