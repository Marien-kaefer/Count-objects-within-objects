openImages = getList("image.titles");


selectWindow("0-1.czi");
run("Split Channels");
selectWindow("0-2.czi");
run("Size...", "width=2560 height=2560 depth=1 constrain average interpolation=None");
run("Merge Channels...", "c1=C1-0-1.czi c2=C2-0-1.czi c4=0-2.czi create");
setSlice(3);
//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");
