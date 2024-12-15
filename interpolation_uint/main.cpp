// Copy below to use interpolation functions. 
// Parameters can be changed in 
// img is interlaced, grayscale, uint8 image.  

bool lensCorr = true;

// Determine rotaion coordinates
parameters = coordinateLSQ(img);
if (parameters[0] == 0) {
    // Default paramters if LSQ doesnt work
    parameters = { 376, -240, -float(1.4968 * pi / 180.0) };
}
    
// Simple image rotation
Mat img_shift = imShiftSimple(img, parameters);

// Interpolation
Mat img_coord = imShiftCoord(img, parameters, lensCorr);
Mat img_intp = interpolateNN(img_coord, 4, lensCorr);
//Mat img_intp = imread("intp_tmp.png", IMREAD_GRAYSCALE);
imwrite("intp_tmp.png", img_intp);  // Save file for quick load when testing other functions 

// Print flux values of objects in image
blobFlux(img_intp, 40, 15);

// Intensity plot of given coordinate (x,y)
intensityPlot(img_intp, 498, 304);
