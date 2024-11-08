//
// Calculates rotation point and angle using least squares
// input: interlaced image
// output: array with 3 floats: {x, y, theta} 
// In addition to the two functions beneath (splitbinary and blobCoordinates), dillate/errode functions are needed (see seperate document).  
// 

array<float, 3> coordinateLSQ(Mat img) {

	// Input: image from JUNO starcam, 
	// Return: array with (x,y,theta) for rotation

	pair<Mat, Mat> img_binary = splitBinary(img);

	// cout << "Stars in first image: " << endl;
	pair<Mat, Mat> coord_two = blobCoordinates(img_binary.first);

	// cout << endl << "Stars in second image: " << endl;
	pair<Mat, Mat> coord_one = blobCoordinates(img_binary.second);

	int N = 0;
	int N2 = 0;

	array<float, 3> parameters = { 0, 0, 0 };
	

	// Check that the same amount of stars is found: 
	for (int i = 0; i < coord_one.first.rows; i++) {
		if (coord_one.first.at<float>(i, 0) < 0) {
			break;
		}
		N++;
	}

	for (int i = 0; i < coord_two.first.rows; i++) {
		if (coord_two.first.at<float>(i, 0) < 0) {
			break;
		}
		N2++;
	}

	if (N != N2) {
		cout << "Error LSQ: different amount of stars found... " << endl; 
		return parameters; 
	}

	cout << "# stars found:  " << N << endl;

	// Initialize variable: 
	Mat x1(N, 1, CV_32FC1);
	Mat y1(N, 1, CV_32FC1);
	Mat x2(N, 1, CV_32FC1);
	Mat y2(N, 1, CV_32FC1);

	Mat d(2 * N, 1, CV_32FC1);

	Mat g0(2 * N, 1, CV_32FC1);
	
	Mat dg_dt(2 * N, 1, CV_32FC1);
	Mat dg_dx(2 * N, 1, CV_32FC1);
	Mat dg_dy(2 * N, 1, CV_32FC1);

	Mat G(2 * N, 3, CV_32FC1);

	double m_theta = 1;
	double m_x = 1;
	double m_y = 1;

	Mat m(3, 1, CV_32FC1);
	m.at<float>(0, 0) = m_theta;
	m.at<float>(1, 0) = m_x;
	m.at<float>(2, 0) = m_y;


	// Move data into variables 
	for (int i = 0; i < N; i++) {
		x1.at<float>(i, 0) = (coord_one.first.at<float>(i, 0));
		y1.at<float>(i, 0) = (coord_one.second.at<float>(i, 0));

		x2.at<float>(i, 0) = (coord_two.first.at<float>(i, 0));
		y2.at<float>(i, 0) = (coord_two.second.at<float>(i, 0));
	}

	// Put values into d and ones
	for (int i = 0; i < 2 * N; i++) {
		if (i < N) {
			d.at<float>(i, 0) = x2.at<float>(i, 0);
		}
		else {
			d.at<float>(i, 0) = y2.at<float>(i - N, 0);
		}

	}

	for (int j = 0; j < 20; j++) {

		// Calculate g0
		for (int i = 0; i < 2 * N; i++) {
			if (i < N) {
				g0.at<float>(i, 0) = (x1.at<float>(i, 0) - m_x) * cos(m_theta) - (y1.at<float>(i, 0) - m_y) * sin(m_theta) + m_x;
			}
			else {
				g0.at<float>(i, 0) = (x1.at<float>(i - N, 0) - m_x) * sin(m_theta) + (y1.at<float>(i - N, 0) - m_y) * cos(m_theta) + m_y;
			}
		}

		// Calculate dg_dtheta, dg_dx0 and dg_dy0
		for (int i = 0; i < 2 * N; i++) {
			if (i < N) {
				dg_dt.at<float>(i, 0) = -(x1.at<float>(i, 0) - m_x) * sin(m_theta) - (y1.at<float>(i, 0) - m_y) * cos(m_theta);
				dg_dx.at<float>(i, 0) = 1 - cos(m_theta);
				dg_dy.at<float>(i, 0) = sin(m_theta);
			}
			else {
				dg_dt.at<float>(i, 0) = (x1.at<float>(i - N, 0) - m_x) * cos(m_theta) - (y1.at<float>(i - N, 0) - m_y) * sin(m_theta);
				dg_dx.at<float>(i, 0) = -sin(m_theta);
				dg_dy.at<float>(i, 0) = 1 - cos(m_theta);
			}
		}

		// Values into G
		G.col(0) = (dg_dt + 0);
		G.col(1) = (dg_dx + 0);
		G.col(2) = (dg_dy + 0);

		m = m + ((G.t() * G).inv()) * G.t() * (d - g0);

		// cout << "G.Gt() inv: " << endl << ((G * G.t()).inv()).size() << endl;


		m_theta = m.at<float>(0, 0);
		m_x = m.at<float>(1, 0);
		m_y = m.at<float>(2, 0);

	}

	parameters[0] = m.at<float>(1, 0);
	parameters[1] = m.at<float>(2, 0);
	parameters[2] = m.at<float>(0, 0);


	cout << "Rotation parameters: (" << parameters[0] << ", " << parameters[1] << "), theta = " << parameters[2] << endl; 
	cout << endl;

	return parameters;
}

pair<Mat, Mat> blobCoordinates(Mat img) {
	// Input: Image with every other line removed
	// output: Two matrices with x and y coordinates of found blobs
	// Conditions for area can be changed near the bottom of the function
	
	Mat img_tmp;
	img.copyTo(img_tmp);		// Temporary image for analysis
	

	int width = img.cols;		// Parameter setup
	int height = img.rows;

	// Remove speckles
	img_tmp = erodeLR(img_tmp);
	
	// Remove left side (due to noise
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (j < 20) {

				img_tmp.at<uchar>(i, j) = 0;
			}

		}
	}

	// dilate back to make stars uniform 
	img_tmp = dilateCross(img_tmp);
	img_tmp = dilateCross(img_tmp);

	int border = 0;
	int area = 0;

	int x = 0;
	int y = 0;

	int x_sum = 0;
	int y_sum = 0;

	int x_min = 0;
	int x_max = 0;

	int y_min = 0;
	int y_max = 0;

	int x_0 = 0;		// Blob start coordinate for connectedness
	int y_0 = 0;

	int direction = 0;
	bool skip_left = false;

	int val = 0;

	Mat x_pos(20, 1, CV_32FC1);
	Mat y_pos(20, 1, CV_32FC1); 

	int count = 0;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {

			if (img_tmp.at<uchar>(i, j) == 255) {
				
				x_min = j;
				x_max = j;

				y_min = i;
				y_max = i;

				while (!(x == j and y == i)) {
					if (border == 0) {
						x = j;
						y = i;
					}

					border += 1;

					// Check for min max values
					if (x < x_min) {
						x_min = x;
					}
					else if (x > x_max) {
						x_max = x;
					}
					if (y < y_min) {
						y_min = y;
					}
					else if (y > y_max) {
						y_max = y;
					}

					// Find next pixel, goes through 
					// Directions: 0 -> up, 1 -> right, 2 -> down, 3 -> left
					//cout << "step ";

					for (int k = 0; k < 4; k++) {

						skip_left = false;

						// Up
						val = img_tmp.at<uchar>(y - 1, x);
						if (direction == 0 and val == 255) {

							// Update coordinate and direciton, then break
							y = y - 1;
							direction = 3;
							skip_left = true;  // added so that it doesnt trigger left right away
							//cout << "up" << endl;
							break;
						}

						// Right
						val = img_tmp.at<uchar>(y, x + 1);
						if (direction == 1 and val == 255) {

							x = x + 1;
							direction = 0;
							//cout << "right" << endl;
							break;
						}

						// Down
						val = img_tmp.at<uchar>(y + 1, x);
						if (direction == 2 and val == 255) {

							y = y + 1;
							direction = 1;
							//cout << "down" << endl;
							break;
						}

						// Left
						val = img_tmp.at<uchar>(y, x - 1);
						if (direction == 3 and val == 255 and !skip_left) {

							x = x - 1;
							direction = 2;
							//cout << "left" << endl;
							break;
						}

						direction = (direction + 1) % 4;

					}


				}

				// Remove blob from image
				for (int k = y_min; k < y_max + 1; k++) {
					for (int l = x_min; l < x_max + 1; l++) {

						val = img_tmp.at<uchar>(k, l);
						if (val == 255) {
							// Add area
							area += 1;
							x_sum += l;
							y_sum += k;

							img_tmp.at<uchar>(k, l) = 0;
						}




					}

				}

				if (area > 30 and area < 150 and count < 20) {

					x_pos.at<float>(count, 0) = float(x_sum) / area;
					y_pos.at<float>(count, 0) = float(y_sum) / area;

					// cout << "" << int(x_pos.at<float>(count, 0)) << ", " << int(y_pos.at<float>(count, 0)) << ", area: " << area << endl;

					count++;
					if (count == 20) {
						cout << "Too many stars, using first 20... " << endl;
					}
				}
				



				x_sum = 0;
				y_sum = 0;
				area = 0;
				border = 0;

			}

			

		}
	}


	return make_pair(x_pos, y_pos);
}

pair<Mat, Mat> splitBinary(Mat img) {
	// Splits interlaced image into two binaries, keeping the same size. 
	int width = img.cols;		// Parameter setup
	int height = img.rows;

	Mat img_bn1, img_bn2;
	img_bn1 = img > 230;
	img_bn2 = img > 230;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (i % 2 == 0) {
				img_bn2.at<uchar>(i, j) = 0;
			}

			else {
				img_bn1.at<uchar>(i, j) = 0;
			}

		}
	}

	return make_pair(img_bn1, img_bn2);

}
