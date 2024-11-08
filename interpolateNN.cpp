//
// Interpolation finding N nearest neighbours set by arr_size 
// Input: Ouput of imShiftCoord function ( image with [x,y,val] channels )
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

	//namedWindow("Image intp shift vals used", WINDOW_AUTOSIZE);
	//imshow("Image intp shift vals used", vals_used);
	//imwrite("img_intp_used.png", vals_used);

	cout << "\rinterpolated    " << endl << "fraction used from each segment: " << float(used) / (height_i * width_i) << endl << endl;
	return interpolated;
}
