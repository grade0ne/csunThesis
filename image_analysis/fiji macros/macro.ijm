// ==== Setup ====
inputDir = getDirectory("Choose input folder with binary masks");
outputDir = getDirectory("Choose output folder for results");
roiDir = outputDir + "ROIs/";
File.makeDirectory(roiDir);

// ==== Analysis Parameters ====
minArea = 200;
maxArea = 1.0E7;
minAR = 2.5;
maxAR = 10.0;

// ==== Begin Batch ====
fileList = getFileList(inputDir);
for (i = 0; i < fileList.length; i++) {
    name = fileList[i];
    
    // Skip non-image files
    if (!(endsWith(name, ".tif") || endsWith(name, ".tiff"))) {
        print("Skipping: " + name);
        continue;
    }

    path = inputDir + name;
    open(path);
    
    // Ensure something opened
    if (!isOpen(name)) {
        print("Failed to open: " + name);
        continue;
    }

    origTitle = getTitle();
    print("Processing: " + origTitle);

    // --- Binary Cleanup ---
    run("8-bit");
    run("Make Binary");
    run("Fill Holes");
    run("Close-"); // morphological operation, NOT window close

    // --- Analyze Particles ---
    run("Set Measurements...", "area shape feret fit redirect=None decimal=3");
    run("Analyze Particles...", "size=" + minArea + "-" + maxArea +
        " show=Nothing display clear add");

    // --- File Name Base (strip extension) ---
    baseName = name;
    if (endsWith(name, ".tif")) {
        baseName = replace(name, ".tif", "");
    } else if (endsWith(name, ".tiff")) {
        baseName = replace(name, ".tiff", "");
    }

    // --- Save Results and ROIs ---
    roiCount = roiManager("count");
    if (roiCount > 0) {
        roiManager("Measure");
        saveAs("Results", outputDir + baseName + "_results.csv");
        roiManager("Save", roiDir + baseName + "_rois.zip");
    }

    // --- Cleanup ---
    if (isOpen(origTitle)) {
        selectWindow(origTitle);
        close();
    }
    roiManager("Reset");
}
