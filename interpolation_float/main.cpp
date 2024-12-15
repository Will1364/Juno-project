// Copy below to use interpolation functions. 
// Parameters can be changed in 
// img is interlaced, grayscale, CV_32FC1 image with values [0, 1]
// 


bool lensCorr = true;

// Determine rotaion coordinates
array<float, 3> parameters = { 0, 0, 0 };
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
imwrite("intp_tmp.png", img_intp * 255.0);  // Save file for quick load when testing other functions 

// Print flux values of objects in image
blobFlux(img_intp, 40/255.0, 15/255.0);

// Intensity plot of given coordinate (x,y)
intensityPlot(img_intp, 498, 304);

// Create upscaled version of star (~260 * 260 image)
Mat img_starScaled = starPlot(img_intp, 183, 401, 10);
