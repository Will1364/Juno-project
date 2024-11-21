// Bellow is the steps to use the interpolation functions:
// img is interlaced, grayscale, uint8 image.  

// Determine rotaion coordinates
array<float, 3> parameters = coordinateLSQ(img);
   
// Simple image despin
Mat img_shift = imShiftSimple(img, parameters);

// Despin by interpolation
Mat img_coord = imShiftCoord(img, parameters);
Mat img_intp = interpolateNN(img_coord, 8);
