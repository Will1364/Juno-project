//
// Interpolation finding N nearest neighbours set by arr_size 
// Input: Ouput of imShiftCoord function ( image with [x,y,val] channels, see function bellow )
// Output: Grayscale image
// 

Mat interpolateNearestN(Mat img) {
	// Image input in format of imShiftCoord output [x,y,val] channels
	cout << "interpolating...";

	int width = img.cols;					// Parameter setup
	int height = img.rows;

	int padding = 0;

	int width_i = width - padding;
	int height_i = height - padding;
	int search = 50;

	int s_x = 0;
	int s_y = 0;

	double tempx = 0;
	double tempy = 0;
	double dist2 = 0;

	const int arr_size = 3;
	float min_d[arr_size];
	int min_i[arr_size];
	int min_j[arr_size];

	double val = 0;
	double tmp = 0;
	int used = 0;

	Mat interpolated(height_i, width_i, CV_8UC1);

	Mat vals_used(height, width, CV_8UC1);
	vals_used = 255;

	// Loop through entire image
	for (int i = 0; i < height_i; i++) {
		for (int j = 0; j < width_i; j++) {

			for (int k = 0; k < arr_size; k++) {
				min_d[k] = 9990;
				min_i[k] = 0;
				min_j[k] = 0;
			}

			s_y = i + padding - (search / 2);
			s_x = j + padding - (search / 2);

			// Fills search window with distance to point in a search x search grid around 
			for (int k = 0; k < search; k++) {
				for (int l = 0; l < search; l++) {

					if (s_y + k > 0 and s_y + k < height and s_x + l > 0 and s_x + l < width) {
						// Test if search coordinate is within window and != 0
						tempx = img.at<Vec3d>(s_y + k, s_x + l)[0] - j;
						tempy = img.at<Vec3d>(s_y + k, s_x + l)[1] - i;

						dist2 = tempx * tempx + tempy * tempy;

						for (int m = 0; m < arr_size; m++) {
							if (min_d[m] > dist2) {
								for (int n = arr_size - 1; n > m; n--) {
									min_d[n] = min_d[n - 1];
									min_i[n] = min_i[n - 1];
									min_j[n] = min_j[n - 1];
								}

								min_d[m] = dist2;
								min_i[m] = s_y + k;
								min_j[m] = s_x + l;

								break;
							}
						}
					}
				}
			}

			if (min_i[0] % 2 == 0) {
				used += 1;
			}

			if (min_d[0] < 1) {
				 
				val = img.at<Vec3d>(min_i[0], min_j[0])[2];

				interpolated.at<uchar>(i, j) = val;
				vals_used.at<uchar>(min_i[0], min_j[0]) = 0;

			}
			else {
				interpolated.at<uchar>(i, j) = 0;
			}

			val = 0;
			tmp = 0;

		}
	}

	// Show image with the values used in interpolation
	//namedWindow("Image intp shift vals used", WINDOW_AUTOSIZE);
	//imshow("Image intp shift vals used", vals_used);
	//imwrite("img_intp_used.png", vals_used);

	cout << "\rinterpolated    " << endl << "fraction used from each segment: " << float(used) / (height_i * width_i) << endl << endl;
	return interpolated;
}


//
// imShiftCoord function
// Rotates image based on rot_prm (x, y, theta), and saves the rotated coordinates
// Input: Interlaced image and rotation parameters (as outputted from LSQ function) 
// Output: Matrix with [x, y, val] channels 
// For use with interpolation function
Mat imShiftCoord(Mat img, array<float, 3> rot_prm) {
	// Creates a Mat object with [x,y,val] (insted of RGB) 
	// Same math as is imShiftSimple, but rotates each part of the image at half the angle 

	double angle = rot_prm[2];   
	angle = angle / 2;

	float x0 = rot_prm[0];							
	float y0 = rot_prm[1];

	int width = img.cols;					// Parameter setup
	int height = img.rows;

	double x_coord = 0;
	double y_coord = 0;

	int x_shift = 0;
	int y_shift = 0;

	Mat img_shift(height, width, CV_64FC3);

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {

			
			if (i % 2 == 1) {

				x_coord = ((j - x0) * cos(angle) - (i - y0) * sin(angle)) + x0;
				x_shift = x_coord;

				y_coord = ((j - x0) * sin(angle) + (i - y0) * cos(angle)) + y0;
				y_shift = y_coord + (1 - int(y_coord) % 2);

				img_shift.at<Vec3d>(i, j)[0] = x_coord;
				img_shift.at<Vec3d>(i, j)[1] = y_coord;
				
				img_shift.at<Vec3d>(i, j)[2] = img.at<uchar>(i, j);
				
			}

			else if (i % 2 == 0) {

				x_coord = ((j - x0) * cos(angle) + (i - y0) * sin(angle)) + x0;
				x_shift = x_coord;

				y_coord = (-(j - x0) * sin(angle) + (i - y0) * cos(angle)) + y0;
				y_shift = y_coord + int(y_coord) % 2;

				img_shift.at<Vec3d>(i, j)[0] = x_coord;
				img_shift.at<Vec3d>(i, j)[1] = y_coord;

				img_shift.at<Vec3d>(i, j)[2] = img.at<uchar>(i, j);
			}
		}
	}

	return img_shift;
}
