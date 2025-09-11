inputDir = "G:/Research/imageAnalysis/segmented images/";
fileName = "c7-17-1_0002_sm_SEG.tiff"; // Replace with one of your test files

open(inputDir + fileName);
origTitle = getTitle();
print("Opened: " + origTitle);

// Convert to 16-bit so label values are preserved
run("16-bit");

// Get dimensions
width = getWidth();
height = getHeight();

// Create binary image
newImage("rotifer_only", "8-bit black", width, height, 1);

// Copy only label-1 pixels into the binary image
for (y = 0; y < height; y++) {
    for (x = 0; x < width; x++) {
        val = getPixel(x, y);
        if (val == 1) {
            setPixel(x, y, 255);
        }
    }
}

// Select the new image and threshold it explicitly
selectWindow("rotifer_only");
run("Threshold...", "set min=255 max=255");
setOption("BlackBackground", true);
run("Convert to Mask");

// Binary cleanup
run("Fill Holes");
run("Close");

// Set measurements
run("Set Measurements...", "area shape feret fit redirect=None decimal=3");
run("Analyze Particles...", "size=200-Infinity show=Nothing display clear add");

// AR filtering
roiCount = roiManager(
