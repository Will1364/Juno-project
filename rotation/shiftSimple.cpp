//
// Simple shift of image
// Input: Interlaced image, array with rotation point and angle eg. array<float, 3> parameters = {373.468, -240.482, -0.0262431};
// Line for header file: Mat imShiftSimple(Mat img);
// 

Mat imShiftSimple(Mat img, array<float, 3> rot_prm) {

	float angle = rot_prm[2];        // Estimated rotation in test image, approximately 1.4 degrees
	      
	float x0 = rot_prm[0];							// Estimated origin of test image
	float y0 = rot_prm[1];

	int width = img.cols;					// Parameter setup
	int height = img.rows;

	int x_shift = 0;
	int y_shift = 0;

	Mat img_shift(height, width, 0);

	// Loop through every pixel in the image
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			// 
			if (i % 2 == 1) {
				img_shift.at<uchar>(i, j) = img.at<uchar>(i, j);
			}
			else if (i % 2 == 0) {

				x_shift = ((j - x0) * cos(angle) - (i - y0) * sin(angle)) + x0;
				y_shift = ((j - x0) * sin(angle) + (i - y0) * cos(angle)) + y0;
				y_shift = y_shift + y_shift % 2;

				if (x_shift >= 0 and x_shift < width and y_shift >= 0 and y_shift < height) {
					// Checks if shifted coordinates is within image, else sets val == 0
					img_shift.at<uchar>(i, j) = img.at<uchar>(y_shift, x_shift);
					// img_used.at<uchar>(y_shift, x_shift) = 255;
				}
				else {
					img_shift.at<uchar>(i, j) = 0;
				}
			}
		}
	}

	return img_shift;

}
