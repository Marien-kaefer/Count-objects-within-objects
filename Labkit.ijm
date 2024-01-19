folder = "Z:/private/Marie/Image Analysis/2023-06-14-LIU-Mengru-count-spots-in-bacteria/input/Labkit"
for (i = 0; i < 3; i++) {
   open(folder  + File.separator + i + "-2" + ".czi");
   file = getTitle();
   run("Segment Image With Labkit", "segmenter_file=" + file + " [Z:/private/Marie/Image Analysis/2023-06-14-LIU-Mengru-count-spots-in-bacteria/input/Labkit/bacteria_TL_segmentation.classifier] use_gpu=false");
   saveAs("Tiff", folder + "segmentation_" + i + ".tif");
}
close("*");


















run("Segment Images in Directory with Labkit", "logger=org.scijava.log.DefaultLogger@27b26b56 input_directory=[Z:\\private\\Marie\\Image Analysis\\2023-06-14-LIU-Mengru-count-spots-in-bacteria\\input\\Labkit] file_filter=*.czi output_directory=[Z:\\private\\Marie\\Image Analysis\\2023-06-14-LIU-Mengru-count-spots-in-bacteria\\input\\Labkit] output_file_suffix=_segmentation.tif use_gpu=false");

