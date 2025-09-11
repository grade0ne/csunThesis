// ==== Setup ====
inputDir = getDirectory("Choose input folder with binary masks");
outputDir = getDirectory("Choose output folder for results");
roiDir = outputDir + "ROIs/";
skeletonDir = outputDir + "SkeletonData/";
shapeDataDir = outputDir + "ShapeData/";
File.makeDirectory(roiDir);
File.makeDirectory(skeletonDir);
File.makeDirectory(shapeDataDir);

// ==== Parameters ====
minArea = 200;
maxArea = 1.0E7;

// ==== Utility ====
function closeAllResultsTables() {
    titles = getList("window.titles");
    for (i = 0; i < titles.length; i++) {
        title = titles[i];
        if (indexOf(title, "Results") >= 0 || indexOf(title, "branchinfo") >= 0) {
            selectWindow(title);
            close();
        }
    }
}

// ==== Begin Batch ====
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

        newImage("ROI_temp", "8-bit black", getWidth(), getHeight(), 1);
        setForegroundColor(255, 255, 255);
        setBackgroundColor(0, 0, 0);

        for (r = 0; r < roiCount; r++) {
            if (isOpen("ROI_temp")) {
                selectWindow("ROI_temp");
                run("Select All");
                run("Clear");
            }

            roiManager("Select", r);
            if (isOpen("ROI_temp")) {
                selectWindow("ROI_temp");
                run("Fill", "slice");

                run("8-bit");
                setAutoThreshold("Default");
                run("Convert to Mask");
                run("Skeletonize");

                // ✅ Run Analyze Skeleton silently
                run("Analyze Skeleton (2D/3D)", "prune=none");

                // ✅ Save any table containing "branchinfo"
                titles = getList("window.titles");
                for (j = 0; j < titles.length; j++) {
                    t = titles[j];
                    if (indexOf(t, "branchinfo") >= 0) {
                        selectWindow(t);
                        saveAs("Results", skeletonDir + baseName + "_branchinfo_roi" + r + ".csv");
                        close();
                        break;
                    }
                }

                closeAllResultsTables();
            }
        }

        if (isOpen("ROI_temp")) {
            selectWindow("ROI_temp");
            close();
        }
    }

    if (isOpen(origTitle)) {
        selectWindow(origTitle);
        close();
    }

    roiManager("Reset");
}

// Final cleanup
closeAllResultsTables();
run("Close All");
setBatchMode(false);
