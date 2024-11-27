// Copy below to use interpolation functions. 
// Parameters can be changed in 
// img is interlaced, grayscale, uint8 image.  

// Determine rotaion coordinates
array<float, 3> parameters = coordinateLSQ(img);
//array<float, 3> parameters = { 373.468, -240.483, - float(1.50362 * pi / 180) };

// Simple image rotation
// Mat img_shift = imShiftSimple(img, parameters);

// Interpolation
Mat img_coord = imShiftCoord(img, parameters);
Mat img_intp = interpolateNN(img_coord, 8);
// Mat img_intp = imread("intp_tmp.png", IMREAD_GRAYSCALE);
imwrite("intp_tmp.png", img_intp);  // Save file for quick load when testing other functions 

// Flux
blobFlux(img_intp, 40, 15);

// Plot intensities
intensityPlot(img_intp, 182, 404);
intensityPlot(img_intp, 518, 373);

    
