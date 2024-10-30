//
// Simple shift of image
// Input: Interlaced image
// Line for header file: Mat imShiftSimple(Mat img);
// 

Mat imShiftSimple(Mat img) {

	float angle = -1.4 * 3.1415/180;    // Estimated rotation in test image, approximately 1.4 degrees
	      
	int x0 = 352;							          // Estimated origin of test image
	int y0 = -297;

	int width = img.cols;					      // Parameter setup
	int height = img.rows;

	int x_shift = 0;    
	int y_shift = 0;

	Mat img_shift(height, width, 0);    // New image 

	// Loop through every pixel in the image
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {

			// Odd lines are preserved
			if (i % 2 == 1) {

				img_shift.at<uchar>(i, j) = img.at<uchar>(i, j);

			}

      // Even lines are shifted (using basic rotaion matrix)
			else if (i % 2 == 0) {

				x_shift = ((j - x0) * cos(angle) - (i - y0) * sin(angle)) + x0;
				y_shift = ((j - x0) * sin(angle) + (i - y0) * cos(angle)) + y0;
				y_shift = y_shift + y_shift % 2;      // Forces even line in y_shift -> doesn't duplicate data

				if (x_shift >= 0 and x_shift < width and y_shift >= 0 and y_shift < height) {
					// Checks if shifted coordinates is within image, else sets val == 0
					img_shift.at<uchar>(i, j) = img.at<uchar>(y_shift, x_shift);

				}
				else {
					img_shift.at<uchar>(i, j) = 0;
				}

			}

		}



	}

	return img_shift;

}
