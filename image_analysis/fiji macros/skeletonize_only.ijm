// ==== Setup ====
inputDir = getDirectory("Choose input folder with binary masks");
outputDir = getDirectory("Choose output folder for results");
roiDir = outputDir + "ROIs/";
shapeDataDir = outputDir + "ShapeData/";
skeletonDir = outputDir + "Skeletons/";

File.makeDirectory(roiDir);
File.makeDirectory(shapeDataDir);
File.makeDirectory(skeletonDir);

// ==== Parameters ====
minArea = 200;
maxArea = 1.0E7;

// ==== Begin ====
setBatchMode(true);

fileList = getFileList(inputDir);
for (i = 0; i < fileList.length; i++) {
    name = fileList[i];
    if (!(endsWith(name, ".tif") || endsWith(name, ".tiff"))) continue;

    path = inputDir + name;
    open(path);
    if (!isOpen(name)) {
        print("Failed to open: " + name);
        continue;
    }

    origTitle = getTitle();
    print("Processing: " + origTitle);

    run("Set Measurements...", "area shape feret fit redirect=None decimal=3");
    run("Analyze Particles...", "size=" + minArea + "-" + maxArea + " show=Nothing clear add");

    baseName = replace(name, ".tif", "");
    baseName = replace(baseName, ".tiff", "");
    roiCount = roiManager("count");

    if (roiCount > 0) {
        roiManager("Measure");
        saveAs("Results", shapeDataDir + baseName + "_results.csv");
        roiManager("Save", roiDir + baseName + "_rois.zip");

        for (r = 0; r < roiCount; r++) {
            // Always start with a fresh canvas
            if (isOpen("ROI_temp")) {
                selectWindow("ROI_temp");
                close();
            }
            newImage("ROI_temp", "8-bit black", getWidth(), getHeight(), 1);

            roiManager("Select", r);
            selectWindow("ROI_temp");
            run("Fill", "slice");

            run("8-bit");
            setAutoThreshold("Default");
            run("Convert to Mask");
            run("Skeletonize");

            saveAs("Tiff", skeletonDir + baseName + "_roi" + r + "_skeleton.tif");

            close("ROI_temp");
        }
    }

    if (isOpen(origTitle)) {
        selectWindow(origTitle);
        close();
    }

    roiManager("Reset");
}

run("Close All");
setBatchMode(false);
