run("Select None");
roiManager("reset");
run("Hide Overlay");

kurtosis_cutoff = 3; 
prominence = newArray(200, 50); 

raw_ID = getImageID();
raw_title = getTitle();
file_path = getDirectory("image");
file_base_title = file_name_remove_extension(raw_title); 


run("Duplicate...", "duplicate");
duplicate_title = getTitle();
run("Split Channels");
selectImage("C1-" + duplicate_title);
run("Bleach Correction", "correction=[Exponential Fit]");
bleach_corrected_Ch1_ID = getImageID();
Stack.getDimensions(width, height, channels, slices, frames);
selectImage("y = a*exp(-bx) + c");
close();
selectImage("C2-" + duplicate_title);
run("Bleach Correction", "correction=[Exponential Fit]");
bleach_corrected_Ch2_ID = getImageID();
Stack.getDimensions(width, height, channels, slices, frames);
selectImage("y = a*exp(-bx) + c");
close();
run("Merge Channels...", "c2=[DUP_C1-" + duplicate_title + "] c6=[DUP_C2-" + duplicate_title + "] create");
bleach_corrected_ID = getImageID();
Stack.getDimensions(width, height, channels, slices, frames);

/*
run("Bleach Correction", "correction=[Exponential Fit]");
bleach_corrected_ID = getImageID();
Stack.getDimensions(width, height, channels, slices, frames);
selectImage("y = a*exp(-bx) + c");
close();
*/



selectImage(raw_ID);
run("Duplicate...", "title=First-frame duplicate frames=1");
first_frame_title = getTitle();
run("Split Channels");
Ch1_title = "C1-"+ first_frame_title;  
Ch2_title = "C2-"+ first_frame_title; 
imageCalculator("Add create 32-bit", Ch1_title , Ch2_title);
added_channels_ID = getImageID();
selectWindow(Ch1_title);
close();
selectWindow(Ch2_title);
close();

selectImage(added_channels_ID);
run("Enhance Contrast", "saturated=0.35");
run("Median...", "radius=3");

run("Duplicate...", " ");
//run("Maximum...", "radius=2");
added_channels_ID_duplicate = getImageID();

//run("Threshold...");
setAutoThreshold("Huang dark");
setOption("BlackBackground", true);
run("Convert to Mask");
mask_ID = getImageID();

selectImage(added_channels_ID);
run("Find Maxima...", "prominence=200 output=[Segmented Particles]");
segmented_particles_ID = getImageID();
setPasteMode("AND");
selectImage(segmented_particles_ID);
run("Copy");
selectImage(mask_ID);
run("Paste");


//selectImage(segmented_particles_ID);
//close(); 
//selectImage(added_channels_ID_duplicate);
//close();

selectImage(mask_ID);
run("Select None");
run("Convert to Mask");
run("Analyze Particles...", "size=0.3-20 exclude clear add");
roiManager("Save", file_path + File.separator + file_base_title + "_bacteria.zip");

run("Set Measurements...", "area mean shape kurtosis redirect=None decimal=3");

spot_counting(bleach_corrected_ID, channels, frames,kurtosis_cutoff, prominence); 

waitForUser("Done"); 

// ##################### SPOT COUNTING ##################### 

function spot_counting(bleach_corrected_ID, channels, frames, kurtosis_cutoff, prominence){
	i=0; 
	selectImage(bleach_corrected_ID);

	roiManager("Deselect");
	roiManager("multi-measure measure_all");
	
	number_of_bacteria = roiManager("count");
	//number_of_bacteria = 5; 	
	
	for (n = 0; n < number_of_bacteria; n++) {
	    roiManager("Select", n);	    
		//run("Enlarge...", "enlarge=3 pixel");
		roiManager("Update");
	    roiManager("rename", "Object_" + IJ.pad((n), 3));
	}
	roiManager("Save", file_path + File.separator + file_base_title + "_bacteria.zip");

    ROI_number = newArray(number_of_bacteria * frames * channels);
    Frame_number = newArray(number_of_bacteria * frames * channels);
    Channel_number = newArray(number_of_bacteria * frames * channels);
    grey_mean = newArray(number_of_bacteria * frames * channels);
    object_class = newArray(number_of_bacteria * frames * channels); 
    kurtosis = newArray(number_of_bacteria * frames * channels);
    spot_count = newArray(number_of_bacteria * frames * channels);

	selectImage(bleach_corrected_ID);
    run("Remove Overlay");
    
    for(channel_index = 1; channel_index <= channels; channel_index++){    
        for(frame_index = 1; frame_index <= frames; frame_index++){
            for (object_index = 0; object_index < number_of_bacteria;  object_index++) {
                roiManager("Select", object_index); 
                selectImage(bleach_corrected_ID);
                Stack.setChannel(channel_index);
                Stack.setFrame(frame_index);                 
                 
                //run("Find Maxima...", "noise="+prominence+" output=[Count]");
                run("Find Maxima...", "noise="+prominence[(channel_index-1)]+" output=[Count]");
                spot_count[i] = getResult("Count", ((number_of_bacteria * frames * channels) + i));
                //print("spot count: " + spot_count[i]); 
                
                kurtosis[i] = getResult("Kurt", i);
                if (kurtosis[i] > kurtosis_cutoff) {
                    object_class[i] = "spots";
                } 
                else{
                    object_class[i] = "homogeneous";
                }
                grey_mean[i] = getResult("Mean", i);
                ROI_number[i] = object_index; 
                Frame_number[i] = frame_index; 
                Channel_number[i] = channel_index;                 
                i++; 
            }  
        }    
    }
    //Array.show("Measurements", ROI_number, Frame_number, Channel_number, grey_mean, object_class, kurtosis, spot_count);
    assemble_and_save_results(number_of_bacteria, frames, grey_mean, kurtosis, object_class); 
    return number_of_bacteria; 
}

function file_name_remove_extension(raw_title ){
	dotIndex = lastIndexOf(raw_title , "." ); 
	file_name_without_extension = substring(raw_title , 0, dotIndex );
	//print( "Name without extension: " + file_name_without_extension );
	return file_name_without_extension;
}

function assemble_and_save_results(number_of_bacteria, frames, grey_mean, kurtosis, object_class){
	k = 0; 
	run("Clear Results");
	//write results into a new results window 
	for (object_index = 1; object_index <= number_of_bacteria; object_index++){   
    	for(channel_index = 1; channel_index <= channels; channel_index++){
        	for(frame_index = 1; frame_index <= frames; frame_index++){
                setResult("Object index", k, "Object_" + IJ.pad((object_index), 3));
                setResult("Frame", k, frame_index);         
                setResult("Channel", k, channel_index);
                setResult("Spot count", k, spot_count[k]); 
                setResult("Grey mean", k, grey_mean[k]); 
                setResult("Kurtosis", k, kurtosis[k]);
                setResult("Object Class", k , object_class[k]); 
                k++; 
            }
        }
    }
    
	saveAs("Results", file_path + File.separator + file_base_title + "_spot_count_results.csv");
	close("Results");
}