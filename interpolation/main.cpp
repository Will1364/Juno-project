// Copy below to use interpolation functions. 
// Parameters can be changed in 
// img is interlaced, grayscale, uint8 image.  

// Determine rotaion coordinates
array<float, 3> parameters = coordinateLSQ(img);
   
// Simple image despin
Mat img_shift = imShiftSimple(img, parameters);

// Despin by interpolation
Mat img_coord = imShiftCoord(img, parameters);
Mat img_intp = interpolateNN(img_coord, 8);
