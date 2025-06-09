/*


												- Written by Marie Held [mheldb@liverpool.ac.uk] June 2024
												  Liverpool CCI (https://cci.liverpool.ac.uk/)
________________________________________________________________________________________________________________________

BSD 2-Clause License

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*
*/

run("Fresh Start");

//get input parameters
#@ String(value="Please specify the files and parameters to be applied for the processing. \n The input file is expected to have at least two channels: fluorophore 1, fluorophore 2. \n The file can contain other channels which will be ignored.", visibility="MESSAGE") message
#@ File (label = "File for processing:", style = "open") segmentation_inputFile
#@ File (label = "ROI set exported from TrackMate:", style = "open") TrackMate_inputFile
#@ Integer(label="Ch1 Prominence: " , value = 1000) Ch1_prominence
#@ Integer(label="Ch2_prominence: " , value = 1000) Ch2_prominence
#@ Double(label="Intensity histogram kurtosis cutoff value: ", value = 3.0) kurtosis_cutoff
#@ String(label = "Limit number of processed ROIs to 10% (of sum of ROIS across all time frames) for testing purposes?", choices={"No","Yes"}, style="radioButtonHorizontal") limit_ROIS_choice


file_path = File.getDirectory(segmentation_inputFile);
segmentation_title = File.getName(segmentation_inputFile);

run("Bio-Formats Windowless Importer", "open=[" + file_path + File.separator + segmentation_title + "]");
roiManager("Open", TrackMate_inputFile); 
original_ID = getImageID();

function file_name_remove_extension(segmentation_title){
	dotIndex = lastIndexOf(segmentation_title, "." ); 
	file_name_without_extension = substring(segmentation_title, 0, dotIndex );
	//print( "Name without extension: " + file_name_without_extension );
	return file_name_without_extension;
}

original_name_no_ext = file_name_remove_extension(segmentation_title);
directory = File.directory;

selectImage(original_ID);
rename(original_name_no_ext); 

run("Set Measurements...", "area mean standard modal min centroid center shape skewness kurtosis solidity stack display redirect=None decimal=3");

channel_index = newArray(); 
spot_count = newArray(); 
frame_index = newArray(); 
roi_name_array = newArray(); 
skewness = newArray(); 
kurtosis = newArray(); 
area = newArray; 
solidity = newArray();
feret = newArray();
max_grey = newArray();
median_grey = newArray();
mean_grey = newArray();
classification = newArray();

counter = 0; 
counter_another_one = 0; 
ROI_count = roiManager("count");

if (roiManager("count") > 0) {
	for (ch = 1; ch <=2; ch++){
	    //selectWindow("C" + ch + "-" + original_name_no_ext);
	    selectImage(original_ID);
	    if (limit_ROIS_choice == "No") { 
	    	loopStop = ROI_count; 
	    }
	    else if (limit_ROIS_choice == "Yes"){
	    	loopStop = Math.floor(ROI_count/10); 
	    }
	    //print("Loop Stop: " + loopStop); 
	    for (k = 0; k < loopStop; k++){
	    	//print("Frame " + i + ", Channel " + ch + ", ROI " + k);
	    	Stack.getPosition(channel, slice, frame);
	    	roiManager("select", k);
	    	Stack.setChannel(ch); 
	    	//run("Draw", "slice");
	    	List.setMeasurements; 
			skewness[counter] = getValue("Skew");
			kurtosis[counter] = getValue("Kurt"); 
			area[counter] = getValue("Area");
			solidity[counter] = getValue("Solidity");
			if (kurtosis[counter] > kurtosis_cutoff) {
				classification[counter] = "spots";
			} 
			else{
				classification[counter] = "homogeneous";
			}
			feret[counter] = getValue("Feret");
			max_grey[counter] = getValue("Maximum");
			median_grey[counter] = getValue("Median");
			mean_grey[counter] = getValue("Mean");
	    	run("Find Maxima...", "noise=Ch" + ch + "_prominence output=[Count]");
	    	channel_index[counter] = ch;
			spot_count[counter] = getResult("Count", counter);
			frame_index[counter] = frame;
			roi_name_array[counter]  = Roi.getName;
			run("Find Maxima...", "noise=Ch" + ch + "_prominence output=[Point Selection]");
			run("Add Selection...");
			counter++; 
			//print("Counter: " + counter); 
	    }	
	}
	run("Clear Results");
	roiManager("reset");
	run("Select None"); 	
}

Table.create(original_name_no_ext + "_Measurements");
Table.setColumn("Channel index", channel_index);
Table.setColumn("Frame index", frame_index);
Table.setColumn("ROI name", roi_name_array);
Table.setColumn("Spot count", spot_count); 
Table.setColumn("Classification", classification); 
Table.setColumn("Skewness", skewness); 
Table.setColumn("Kurtosis", kurtosis); 
Table.setColumn("Area", area); 
Table.setColumn("Solidity", solidity); 
Table.setColumn("Feret diameter", feret); 
Table.setColumn("Max Grey Value", max_grey); 
Table.setColumn("Median Grey Valye", median_grey); 
Table.setColumn("Mean Grey Value", mean_grey); 

Table.save(directory + File.separator + original_name_no_ext + "_Measurements.csv");
run("Clear Results");

waitForUser("Done!");