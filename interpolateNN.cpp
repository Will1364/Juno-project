//
// Interpolation finding N nearest neighbours set by arr_size 
// Input: Ouput of imShiftCoord function ( image with [x,y,val] channels, see function bellow )
// Output: Grayscale image
// 

Mat interpolateNearestN(Mat img, int arr_size) {
	// Image input in format of imShiftCoord output [x,y,val] channels
	cout << "interpolating";


	int width = img.cols;					// Parameter setup
	int height = img.rows;

	int padding = 0;

	int width_i = width - padding;
	int height_i = height - padding;

	int search = 100;		// NxN search area around each pixel. 50 seems to be a good value

	int s_x = 0;
	int s_y = 0;

	//double tempx = 0;
	//double tempy = 0;
	double dist = 0;
	
	vector<float> min_d(arr_size, 0);
	vector<float> min_v(arr_size, 0);
	
	vector<int> min_i(arr_size, 0);
	vector<int> min_j(arr_size, 0);

	double val = 0;
	double tmp = 0;   
	

	Mat interpolated(height_i, width_i, CV_8UC1);

	Mat vals_used(height, width, CV_8UC1);
	vals_used = 250;	// To display number of pixels closer than 1 

	int used = 0;		// Counts used pixels from even rows
	int unused = 0;		// Counts number of unused pixels
	int progress = 1;	// Used to display progress

	// Loop through entire image
	for (int i = 0; i < height_i; i++) {
		for (int j = 0; j < width_i; j++) {

			for (int k = 0; k < arr_size; k++) {
				min_d[k] = 9990;
				min_i[k] = 0;
				min_j[k] = 0;
				min_v[k] = 0;
			}

			s_y = i + padding - (search / 2);
			s_x = j + padding - (search / 2);

			// Fills search window with distance to point in a search x search grid around 
			for (int k = 0; k < search; k++) {
				for (int l = 0; l < search; l++) {

					if (s_y + k > 0 and s_y + k < height and s_x + l > 0 and s_x + l < width) {
						// Test if search coordinate is within window
						
						dist = sqrt(pow(img.at<Vec3d>(s_y + k, s_x + l)[0] - j, 2) + pow(img.at<Vec3d>(s_y + k, s_x + l)[1] - i, 2));

						for (int m = 0; m < arr_size; m++) {
							if (min_d[m] > dist) {
								for (int n = arr_size - 1; n > m; n--) {
									min_d[n] = min_d[n - 1];
									min_i[n] = min_i[n - 1];
									min_j[n] = min_j[n - 1];
								}

								min_d[m] = dist;
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

				for (int k = 0; k < arr_size; k++) {
					min_v[k] = img.at<Vec3d>(min_i[k], min_j[k])[2];
				}
				
				// val = img.at<Vec3d>(min_i[0], min_j[0])[2];			// Nearest neighbour
				val = intpCoreDist(min_d, min_v);					// Inverse distance weighting 
				// val = intpCoreRadialBasis(j, i, min_d, min_i, min_j, min_v, img);		// Radial basis interpolation
				
				
				

				interpolated.at<uchar>(i, j) = int(val);
				vals_used.at<uchar>(min_i[0], min_j[0]) -= 50;

			}
			else {
				interpolated.at<uchar>(i, j) = 0;
				unused++;
			}

			if (i > progress * height / 20) {
				progress++;
				cout << ".";
			}

			val = 0;
			tmp = 0;

		}
	}

	//namedWindow("Image intp shift vals used", WINDOW_AUTOSIZE);
	//imshow("Image intp shift vals used", vals_used);
	//imwrite("img_intp_used.png", vals_used);

	cout << "\rinterpolated          " << endl;
	cout << "fraction used from each segment: " << float(used) / (height_i * width_i) << endl;
	cout << "fraction unused: " << float(unused) / (height_i * width_i) << endl << endl;
	
	return interpolated;
}

// 
// Interpolation cores -> for use with interpolation funtion
// 

double intpCoreDist(vector<float> min_d, vector<float> min_v) {
	// Inverse distance weight for interpolation

	double sumW = 0;
	double sumWu = 0;
	int p = 2;


	for (int k = 0; k < min_d.size(); k++) {
		if (min_d[k] < 1) {
			sumWu += min_v[k] / pow(min_d[k], p);
			sumW += 1 / pow(min_d[k], p);
		}
	}
	if (sumW != 0) {
		return sumWu/sumW;
	}
	else {
		return 0.0;
	}
	
}

double intpCoreRadialBasis(int x_, int y_, vector<float> min_d, vector<int> min_i, vector<int> min_j, vector<float> min_v, Mat img) {
	// Local Radial basis interpolation based on -----
	// Uses n nearest points 
	// 

	float eps = 0.5;		// Shape parameter -> used to scale the thingy  
	int N = min_d.size();

	double r;


	Mat f(N, 1, CV_64FC1);
	Mat phi(N, N, CV_64FC1);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			
			r = sqrt(pow(img.at<Vec3d>(min_i[i], min_j[i])[0] - img.at<Vec3d>(min_i[j], min_j[j])[0], 2) + pow(img.at<Vec3d>(min_i[i], min_j[i])[1] - img.at<Vec3d>(min_i[j], min_j[j])[1], 2));
			
			phi.at<double>(i, j) = exp(-pow(eps * abs(r), 2));		// Gaussian 
			// phi.at<double>(i, j) = pow(r, 2) * log(r * eps);				// Thin plate spline
			// phi.at<double>(i, j) = 1 / (1 + pow(eps * r, 2));			// Inverse quadratic
		}

		f.at<double>(i, 0) = min_v[i];

	}

	Mat w = phi.inv() * f;
	double s = 0;				// Intp value
	//double sum_phi = 0;			// Normalization factor

	for (int i = 0; i < N; i++) {
		r = min_d[i];
		s += w.at<double>(i, 0) * exp(-pow(eps * r, 2));				// Gaussian
		// s += w.at<double>(i, 0) * pow(r, 2) * log(r * eps);				// Thin plate spline
		// s += w.at<double>(i, 0) * 1 / (1 + pow(eps * r, 2));			// Inverse quadratic

		//sum_phi += exp(-pow(eps * r, 2));			
	}

	/*
	if (x_ == 259 and y_ == 87) {
		cout << "s = " << s << endl;
		cout << "phi: " << endl << phi.inv() << endl;
		cout << "min_d: " << endl << min_d[0] << ", " << min_d[1] << ", " << ", " << endl;
	}
	*/

	

	if (s > 0 and s < 256) {
		return s;
	}
	else if (s < 0) {
		return 0;
	}
	else if (s > 255) {
		return 255;
	}
	
	return 0;
}


//
// imShiftCoord function
// Rotates image based on rot_prm (x, y, theta), and saves the rotated coordinates
// Input: Interlaced image and rotation parameters (as outputted from LSQ function) 
// Output: Matrix with [x, y, val] channels 
// For use with interpolation function
// 

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
