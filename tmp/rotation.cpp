#include "header.h"

// Constants 
const double pi = 3.14159265358979;

Mat imReadGS(string img_name) {
	Mat img;

	string path = "G:/My Drive/UNI/30330 Image Analysis with Microcomputer/30330 Project/pngs/";

	string path_full = path + img_name;

	img = imread(path_full, IMREAD_GRAYSCALE);

	cout << "Image loaded. \n\n";
	cout << "filename:    " << img_name << "\n";
	cout << "full path:   " << path_full << "\n";
	cout << "Image size:  " << img.size() << "\n";
	cout << endl;
	return img;

}

array<int, 256> histArrayGS(Mat img) {

	unsigned char value = 0; 			// index value for the histogram
	array<int, 256> histogram; 			// histogram array - remember to set to zero initially
	int width = img.cols; 				// INIT USING IMAGE VALUES
	int height = img.rows; 				// INIT USING IMAGE VALUES

	int k = 256;						// Loop variable -> goes over every image value

	// Subtracts one from k, then checks if it is greater than zero
	while (k-- > 0) {
		histogram[k] = 0; 				// reset histogram entry
	}

	// Loops over image and adds one to histogram bin for given pixel value
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {

			value = img.at<uchar>(i, j);
			histogram[value] += 1;

		}
	}



	return histogram;
}

void histDisplay(array<int, 256> histogram, string windowName) {

	Mat histImg(260, 260, CV_8U, 1);

	int k = 256;
	int n = *max_element(histogram.begin(), histogram.end());
	int y0, y1;

	Point pt1, pt2;

	int offset = 5;

	// Draw vertical lines
	int l = 10;
	int dcol;
	int nrows;
	while (l-- > 0) {

		nrows = histImg.rows;
		dcol = int(float(histImg.cols) / 10);

		pt1 = Point2d(dcol * l + offset, 0);
		pt2 = Point2d(dcol * l + offset, nrows);
		line(histImg, pt1, pt2, 100);


	}


	// Draw histogram lines
	while (k-- > 1) {

		y0 = 255 - int(float(histogram[k - 1]) / n * 255.0);
		y1 = 255 - int(float(histogram[k]) / n * 255.0);

		pt1 = Point2d(offset + k, y0);
		pt2 = Point2d(offset + k + 1, y1);

		line(histImg, pt1, pt2, 255);

	}

	namedWindow(windowName, WINDOW_AUTOSIZE);
	imshow(windowName, histImg);

	return;
}

pair<Mat, Mat> splitImage(const Mat& OGimage) {
	int width = OGimage.cols; //width of image
	int height = OGimage.rows; // hight of image

	int n = height / 2;
	int maxn;
	Mat splitImg1(n, width, OGimage.type());  // creates empty image with half height and full width of original image.
	Mat splitImg2(n, width, OGimage.type());  // The new images are the same type as the original
	// creating empty Mat structures  with specified sizes are more efficient than ones with unspecified dimensions

	for (int i = 0; i < n; i++) {

		OGimage.row(i * 2).copyTo(splitImg1.row(i));      // Copy even rows to splitImg1 in a safe way
		OGimage.row(i * 2 + 1).copyTo(splitImg2.row(i));  // Copy odd rows to splitImg2


		/*
		splitImg1.row(i) = OGimage.row(i * 2);			// simpler way to copy pixels, though dangerous as
		splitImg2.row(i) = OGimage.row(i * 2 + 1);		// the Mat structures will share memory with the original this way (a seperate copy is not actualy created/shallow copy)
		*/

		maxn = i * 2 + 1;
	}
	cout << "maxn = " << maxn << "\n";

	return make_pair(splitImg1, splitImg2);

}

pair<Mat, Mat> splitBinary(Mat img) {

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

Mat erodeLR(Mat img) {

	int val;
	Mat img_tmp;
	img.copyTo(img_tmp);

	int width = img.cols;		// Parameter setup
	int height = img.rows;

	for (int i = 0; i < height; i++) {
		for (int j = 1; j < width - 1; j++) {

			val = img.at<uchar>(i, j);

			if (val == 0) {

				img_tmp.at<uchar>(i, j - 1) = 0;
				img_tmp.at<uchar>(i, j) = 0;
				img_tmp.at<uchar>(i, j + 1) = 0;

			}

		}

	}

	return img_tmp;

}

Mat dilateCross(Mat img) {

	int val;
	Mat img_tmp;
	img.copyTo(img_tmp);

	int width = img.cols;		// Parameter setup
	int height = img.rows;

	for (int i = 1; i < height - 1; i++) {
		for (int j = 1; j < width - 1; j++) {

			val = img.at<uchar>(i, j);

			//cout << val << endl;

			if (val == 255) {

				//cout << "E" << endl;

				img_tmp.at<uchar>(i, j - 1) = 255;
				img_tmp.at<uchar>(i - 1, j) = 255;
				img_tmp.at<uchar>(i, j + 1) = 255;
				img_tmp.at<uchar>(i + 1, j) = 255;

			}

		}

	}

	return img_tmp;

}

Mat imShiftSimple(Mat img, array<float, 3> rot_prm) {

	float angle = rot_prm[2];        // Estimated rotation in test image, approximately 1.4 degrees
	      
	float x0 = rot_prm[0];							// Estimated origin of test image
	float y0 = rot_prm[1];

	int width = img.cols;					// Parameter setup
	int height = img.rows;

	int x_shift = 0;
	int y_shift = 0;

	Mat img_shift(height, width, 0);
	Mat vals_used(height, width, CV_8UC1);
	vals_used = 255;

	// Loop through every pixel in the image
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {

			// 
			if (i % 2 == 1) {

				img_shift.at<uchar>(i, j) = img.at<uchar>(i, j);
				vals_used.at<uchar>(i, j) = 0;
			}

			else if (i % 2 == 0) {

				x_shift = ((j - x0) * cos(angle) - (i - y0) * sin(angle)) + x0;
				y_shift = ((j - x0) * sin(angle) + (i - y0) * cos(angle)) + y0;
				y_shift = y_shift + y_shift % 2;

				if (x_shift >= 0 and x_shift < width and y_shift >= 0 and y_shift < height) {
					// Checks if shifted coordinates is within image, else sets val == 0

					img_shift.at<uchar>(i, j) = img.at<uchar>(y_shift, x_shift);
					vals_used.at<uchar>(y_shift, x_shift) = 0;

				}
				else {
					img_shift.at<uchar>(i, j) = 0;
				}

			}

		}

	}

	
	//namedWindow("Vals used in rotation", WINDOW_AUTOSIZE);
	//imshow("Vals used in rotation", vals_used);
	

	return img_shift;

}

Mat imShiftCoord(Mat img, array<float, 3> rot_prm) {
	// Creates a Mat object with [x,y,val] (insted of RGB) 
	// Same math as is imShiftSimple

	double angle = rot_prm[2];   // Estimated rotation in test image, approximately 1.4 degrees
	angle = angle / 2;

	float x0 = rot_prm[0];							// Estimated origin of test image
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

pair<Mat, Mat> blobCoordinates(Mat img) {
	// Input: Image with every other line removed

	Mat img_tmp;
	img.copyTo(img_tmp);		// Temporary image for analysis

	int width = img.cols;		// Parameter setup
	int height = img.rows;

	// Remove speckles
	img_tmp = erodeLR(img_tmp);
	
	// Remove left side (due to noise
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (j < 30) {

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

	for (int i = 1; i < height-1; i++) {
		for (int j = 1; j < width-1; j++) {

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

					//cout << "" << int(x_pos.at<float>(count, 0)) << ", " << int(y_pos.at<float>(count, 0)) << ", area: " << area << endl;

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

array<float, 3> coordinateLSQ(Mat img) {

	// Input: image from JUNO starcam, 
	// Return: array with (x,y,theta) for rotation

	pair<Mat, Mat> img_binary = splitBinary(img);

	//cout << "Stars in first image: " << endl;
	pair<Mat, Mat> coord_two = blobCoordinates(img_binary.first);

	//cout << endl << "Stars in second image: " << endl;
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
		cout << "Error LSQ: different amount of objects found, returning zeroes... " << endl; 
		return parameters; 
	}

	cout << "objects found:  " << N << endl;

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


	cout << "rotation parameters: (" << parameters[0] << ", " << parameters[1] << ", " << parameters[2] * 180/pi << ") " << endl;
	cout << endl;

	return parameters;
}

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

	
	if (x_ == 259 and y_ == 87) {
		cout << "s = " << s << endl;
		cout << "phi: " << endl << phi.inv() << endl;
		cout << "min_d: " << endl << min_d[0] << ", " << min_d[1] << ", " << ", " << endl;
	}
	

	

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
