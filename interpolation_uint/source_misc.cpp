#include "header.h"

// 
// Various functions used for handling rotaion of image and lens correction.  
// Text for header file: 
/*

// Constants
const double pi = 3.14159265358979;

// Split image
pair<Mat, Mat> splitBinary(Mat img, int threshold);

// Dilation and errosion
Mat dilateCross(Mat img);
Mat erodeLR(Mat img);

// Image shifting
Mat imShiftSimple(Mat img, array<float, 3> rot_prm);
Mat imShiftCoord(Mat img, array<float, 3> rot_prm, bool lensCorrection);     // Creates a 3d array with [x,y,val] along the 3rd axis

// Figures
void intensityPlot(Mat img, int x_peak, int y_peak);

// Other functions
pair<Mat, Mat> blobCoordinates(Mat img);
void blobFlux(Mat img, int threshold, int noise_floor);
array<float, 3> coordinateLSQ(Mat img, bool lensCorrection = false);

Mat eulerRot(int ax1, int ax2, int ax3, float ang1, float ang2, float ang3, bool useDeg = false);
pair<float, float> angleCoord(float x, float y);

*/
// 


pair<Mat, Mat> splitBinary(Mat img, int threshold) {
	// Split image into two binaries
	// For use with coordinateLSQ function
	// 
	// Input: Intervowen image
	// Output: Pair with two binary images 
	//


	// Parameters
	int width = img.cols;		
	int height = img.rows;

	Mat img_bn1, img_bn2;
	img_bn1 = img > threshold;
	img_bn2 = img > threshold;

	// Go through and split image
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
		
			// Even rows
			if (i % 2 + border_pad % 2 == 0) {
				img_bn2.at<uchar>(i, j) = 0;
			}

			// Odd rows
			else {
				img_bn1.at<uchar>(i, j) = 0;
			}


		}
	}

	return make_pair(img_bn1, img_bn2);

}

pair<Mat, Mat> splitUint(Mat img) {
	// Split image into two binaries
	// For use with coordinateLSQ function
	// 
	// Input: Intervowen image
	// Output: Pair with two binary images 
	//


	// Parameters
	int width = img.cols;
	int height = img.rows;

	Mat img_bn1, img_bn2;
	img.copyTo(img_bn1);
	img.copyTo(img_bn2);

	// Go through and split image
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {

			// Even rows
			if (i % 2 + border_pad % 2 == 0) {
				img_bn2.at<uchar>(i, j) = 0;
			}

			// Odd rows
			else {
				img_bn1.at<uchar>(i, j) = 0;
			}


		}
	}

	return make_pair(img_bn1, img_bn2);

}


Mat erodeLR(Mat img) {
	// Errodes image in left and right direction
	// For use with coordinateLSQ/blobCoordinets as a simple noise filter for binary images

	// Parameters 
	int val;
	Mat img_tmp;
	img.copyTo(img_tmp);

	int width = img.cols;		// Parameter setup
	int height = img.rows;

	// Go through image and errode
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
	// Dilation in a small + (up, right, down, left)
	// For use with coordinateLSQ/blobCoordinets to combine stars (stripes in binary images)

	// Parameter setup
	int val;
	Mat img_tmp;
	img.copyTo(img_tmp);

	int width = img.cols;		
	int height = img.rows;

	// Dilation
	for (int i = 1; i < height - 1; i++) {
		for (int j = 1; j < width - 1; j++) {

			val = img.at<uchar>(i, j);

			if (val == 255) {

				img_tmp.at<uchar>(i, j - 1) = 255;
				img_tmp.at<uchar>(i - 1, j) = 255;
				img_tmp.at<uchar>(i, j + 1) = 255;
				img_tmp.at<uchar>(i + 1, j) = 255;

			}
		}
	}

	return img_tmp;
}

pair<Mat, Mat> blobCoordinates(Mat img, int threshold) {
	// Finds blobs (stars) and determines their center coordinate and area. Returns list of coordinates (x, y). 
	// Does simple errosion/dilation to remove noise/combine stripes to full stars, 
	// and 4-connectedness to determine position/area of stars. 
	// Used in coordinateLSQ.
	// 
	// Can optionally print coordinates and area of each found blob. 
	// 
	// Change area conditions bellow if errors occour (or add more/less errosion and dilation) 
	// 
	// Input
	// img:		Binary image with every other line removed
	// 
	// Output
	// Pair with two arrays of coordinates, x and y, for each star 
	//

	// Variables 
	int area_min = 25;		// Minimum area for object to count
	int area_max = 450;		// Maximum area
	int count_max = 20;		// Maximum number of object to use 
	bool use_max_val = false;  // Tells it to use the max val in a star instead of COM
	bool use_COM_weight = true;

	// Paramter setup
	Mat img_tmp = img > threshold;


	// img.copyTo(img_tmp);		// Temporary image for analysis

	int width = img.cols;		// Parameter setup
	int height = img.rows;

	// Remove speckles by erosion
	img_tmp = erodeLR(img_tmp);
	
	// Remove left side (due to noise)
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
	img_tmp = dilateCross(img_tmp);
	img_tmp = dilateCross(img_tmp);
	
	// More parameters
	int border = 0;
	double area = 0;

	int x = 0;
	int y = 0;

	double x_sum = 0;
	double y_sum = 0;

	int x_min = 0;
	int x_max = 0;

	int y_min = 0;
	int y_max = 0;

	int x_0 = 0;		// Blob start coordinate for connectedness
	int y_0 = 0;

	int direction = 0;
	bool skip_left = false;

	int val = 0;		// Binary 

	int star_val = 0;		// For finding maxval
	int star_peak_x = 0;
	int star_peak_y = 0;
	int star_peak = 0;

	double flux = 0;

	// Matrices to hold coordinates. Max set at 20 since it's probably way more than needed. 
	int count = 0;		// Used to count number of stars found
	
	Mat x_pos(count_max + 1, 1, CV_32FC1);
	Mat y_pos(count_max + 1, 1, CV_32FC1);

	

	// Scan through binary image 
	for (int i = 1; i < height-1; i++) {
		for (int j = 1; j < width-1; j++) {

			// Four connectedness when object is found:
			if (img_tmp.at<uchar>(i, j) == 255) {

				// Min/max x of the four connectedness, to cut out object afterwards
				x_min = j;
				x_max = j;

				y_min = i;
				y_max = i;

				// Loop until back at start
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

					// Find next pixel
					// Directions: 0 -> up, 1 -> right, 2 -> down, 3 -> left
					for (int k = 0; k < 4; k++) {

						skip_left = false;

						// Up
						if (y > 0) {
							val = img_tmp.at<uchar>(y - 1, x);
							if (direction == 0 and val == 255) {

								// Update coordinate and direciton, then break
								y = y - 1;
								direction = 3;
								skip_left = true;  // added so that it doesnt trigger left right away
								//cout << "up" << endl;
								break;
							}
						}


						// Right
						if (x < width - 1) {
							val = img_tmp.at<uchar>(y, x + 1);
							if (direction == 1 and val == 255) {

								x = x + 1;
								direction = 0;
								//cout << "right" << endl;
								break;
							}
						}


						// Down
						if (y < height - 1) {
							val = img_tmp.at<uchar>(y + 1, x);
							if (direction == 2 and val == 255 and y < height) {

								y = y + 1;
								direction = 1;
								//cout << "down" << endl;
								break;
							}
						}


						// Left
						if (x > 0) {
							val = img_tmp.at<uchar>(y, x - 1);
							if (direction == 3 and val == 255 and !skip_left and x > 0) {

								x = x - 1;
								direction = 2;
								//cout << "left" << endl;
								break;
							}
						}


						// Rotate start direction 90 degrees
						direction = (direction + 1) % 4;

					}

				}

				// Remove blob from image
				for (int k = y_min; k < y_max + 1; k++) {
					for (int l = x_min; l < x_max + 1; l++) {

						val = img_tmp.at<uchar>(k, l);
						star_val = img.at<uchar>(k, l);


						

						// Compute area and x_sum/y_sum for center value
						if (val == 255) {
							// Add area

							if (use_COM_weight) {
								area += 1;

								x_sum += l * star_val;
								y_sum += k * star_val;

								flux += star_val;
							}
							else {
								area += 1;

								x_sum += l;
								y_sum += k;
							}

							

							img_tmp.at<uchar>(k, l) = 0;
						}

						if (star_val > star_peak) {
							star_peak = star_val;
							star_peak_x = l;
							star_peak_y = k;
						}

					}

				}

				// Save coordinate if area of object is within conditions
				if (area > area_min and area < area_max and count < count_max) {

					if (use_max_val) {
						x_pos.at<float>(count, 0) = star_peak_x;
						y_pos.at<float>(count, 0) = star_peak_y;
					}
					else if (use_COM_weight) {
						x_pos.at<float>(count, 0) = x_sum / flux;
						y_pos.at<float>(count, 0) = y_sum / flux;
					}
					else {
						x_pos.at<float>(count, 0) = float(x_sum) / area;
						y_pos.at<float>(count, 0) = float(y_sum) / area;
					}
					

					//cout << "[x,y]: " << (x_pos.at<float>(count, 0)) << ", " << (y_pos.at<float>(count, 0)) << ", area: " << area << endl;

					count++;
					if (count == count_max) {
						cout << "blobCoordinates: Too many stars, using first " << count_max << "..." << endl;
						return make_pair(x_pos, y_pos);
					}
				}
				
				// Reset parameters
				x_sum = 0;
				y_sum = 0;
				area = 0;
				border = 0;
				star_peak = 0;
				flux = 0;

			}

		}
	}

	
	// Return set of coordinates
	return make_pair(x_pos, y_pos);
}

array<float, 3> coordinateLSQ(Mat img, bool lensCorrection) {
	// Least squares solution to find rotation parameters for image. 
	// Based on non-linear least squares.
	// 
	// Depend on the following functions: 
	// splitBinary
	// blobCoordinates
	//		erodeLR
	//		dilateCross
	// 
	// Inputs 
	// img:				Image, interlaced
	// lensCorrection:  true/false if lens correction should be enabled
	//					currently doesn't work properly
	// 
	// Ouput
	// Roation parameters as array [x, y, theta]
	// 

	// Variables 
	int threshold = 40;
	int dist_max = 30;		// Used to remove objects in only one image if different amount of objects are found

	cout << "LSQ: determining rotation paramters..." << endl;

	// Parameter setup
	pair<Mat, Mat> img_binary = splitUint(img); 




	//cout << "Stars in first image: " << endl;
	pair<Mat, Mat> coord_two = blobCoordinates(img_binary.first, threshold);

	//cout << endl << "Stars in second image: " << endl;
	pair<Mat, Mat> coord_one = blobCoordinates(img_binary.second, threshold);

	
	int N = 0;
	int N2 = 0;

	int tmp = 0;
	float dist = 0;
	float dist_min = 9999;
	

	array<float, 3> parameters = { 0, 0, 0 };
	
	// Lens correction parameters
	float dx = 8.6;				// Physical pixel size [µm]
	float dy = 8.3;
	float pix_ratio = dy / dx;	// Pixel ratio
	int x_l0 = 383 - border_pad;				// Lens center (from Matlab script)
	int y_l0 = 257 - border_pad;
	int efl = 20006;			// Effective focal length 
	double kappa = 3.3e-08;		// distortion correction 
	double f = efl / dx;		// focal length
	double r = 0;


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

	
	if (N == 0 or N2 == 0) {
		cout << "Error LSQ: No objects found in image, returning zeroes..." << endl << endl;
		return parameters; 
	}

	if (N != N2) {
		

		cout << "Error LSQ: different amount of objects found in binary images, " << N << "/" << N2 << ", trying to remove invalid values" << endl;

		// Loop through distance test twice 
		// (In case there are more invalid coords in set two but, but still invalid coords in set 1)
		for (int l = 0; l < 2; l++) {

			// If more objects are found in first image: 
			if (N > N2) {

				// Remove corrdinates for objects with no nearby objects in other image
				for (int i = 0; i < N; i++) {

					dist_min = 9999;

					for (int j = 0; j < N2; j++) {

						// Determine closest coordinate set: 
						dist = pow((coord_one.first.at<float>(i, 0) - coord_two.first.at<float>(j, 0)), 2) + pow((coord_one.second.at<float>(i, 0) - coord_two.second.at<float>(j, 0)), 2);
						dist = sqrt(dist);

						if (dist < dist_min) {
							dist_min = dist;
						}

					}

					if (dist_min > dist_max) {

						for (int k = i; k < N; k++) {
							coord_one.first.at<float>(k, 0) = coord_one.first.at<float>(k + 1, 0);
							coord_one.second.at<float>(k, 0) = coord_one.second.at<float>(k + 1, 0);
						}

					}

				}

				// Calculate number of valid coordinates again
				N = 0;
				for (int i = 0; i < coord_one.first.rows; i++) {
					if (coord_one.first.at<float>(i, 0) < 0) {
						break;
					}
					N++;
				}


			}


			if (N2 > N) {
				// Remove objects with no nearby object in other image
				tmp = 0;
				for (int i = 0; i < N; i++) {

					dist_min = 9999;

					for (int j = 0; j < N2; j++) {

						// Determine closest coordinate set: 
						dist = pow((coord_one.first.at<float>(j, 0) - coord_two.first.at<float>(i, 0)), 2) + pow((coord_one.second.at<float>(j, 0) - coord_two.second.at<float>(i, 0)), 2);
						dist = sqrt(dist);

						if (dist < dist_min) {
							dist_min = dist;
						}

					}

					if (dist_min > dist_max) {

						for (int k = i; k < N + 1; k++) {
							coord_two.first.at<float>(k, 0) = coord_two.first.at<float>(k + 1, 0);
							coord_two.second.at<float>(k, 0) = coord_two.second.at<float>(k + 1, 0);
						}

					}

				}

				// Calculate number of valid coordinates again
				N2 = 0;
				for (int i = 0; i < coord_two.first.rows; i++) {
					if (coord_two.first.at<float>(i, 0) < 0) {
						break;
					}
					N2++;
				}


			}
		}

		if (N != N2) {
			cout << "Error LSQ: Couldn't match objects, returning emtpy array..." << endl << endl;
			return parameters;
		}
		
	}

	cout << "valid objects found: " << N << endl;

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

	double m_theta = -1.5 * pi/180;
	double m_x = 380;
	double m_y = -250;

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

	// Apply lens correction to coordinates 
	if (lensCorrection) {

		float xt = 0;
		float yt = 0;

		for (int i = 0; i < N; i++) {
			// Correction for set 1
			xt = x1.at<float>(i, 0) - x_l0;
			yt = y1.at<float>(i, 0) - y_l0;
			yt = yt * dy / dx;

			r = sqrt((xt - x_l0) * (xt - x_l0) + (yt - y_l0) * (yt - y_l0));
			r = 1 + r * r * kappa;

			x1.at<float>(i, 0) = xt * r + x_l0;
			y1.at<float>(i, 0) = yt * r + y_l0 * dy/dx;


			// Correction for set 2
			xt = x2.at<float>(i, 0);
			yt = y2.at<float>(i, 0);

			r = sqrt((xt - x_l0) * (xt - x_l0) + (yt - y_l0) * (yt - y_l0));
			r = 1 + r * r * kappa;

			x2.at<float>(i, 0) = xt * r;
			y2.at<float>(i, 0) = yt * r;
		}


	}

	// Setup data vecter
	for (int i = 0; i < 2 * N; i++) {
		if (i < N) {
			d.at<float>(i, 0) = x2.at<float>(i, 0);
		}
		else {
			d.at<float>(i, 0) = y2.at<float>(i - N, 0);
		}

	}

	// Non linear least squares computation
	for (int j = 0; j < 50; j++) {

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

		/*
		if (abs(m_x - m.at<float>(1, 0)) / (m.at<float>(1, 0)) < 1e-6) {
			
			m_theta = m.at<float>(0, 0);
			m_x = m.at<float>(1, 0);
			m_y = m.at<float>(2, 0);

			cout << "broke at itteration: " << j << endl;

			break;
		}
		*/

		m_theta = m.at<float>(0, 0);
		m_x = m.at<float>(1, 0);
		m_y = m.at<float>(2, 0);

	}

	parameters[0] = m.at<float>(1, 0);
	parameters[1] = m.at<float>(2, 0);
	parameters[2] = m.at<float>(0, 0);


	// Print parameters to console
	cout << "rotation parameters: (" << parameters[0] << ", " << parameters[1] << ", " << parameters[2] * 180/pi << ") " << endl;
	cout << endl;

	return parameters;
}

Mat imShiftSimple(Mat img, array<float, 3> rot_prm) {
	// Simple rotation of image. Forces use of even/odd lines
	// 
	// Inputs:
	// img:		Interlaced image
	// rot_prm: Rotation parameters [x0, y0, theta] (rotaion point and angle) 
	//			corresponds to output of coordinateLSQ
	// 
	// Output:
	// Despinned image
	// 

	// Parameter setup
	float angle = -rot_prm[2];

	float x0 = rot_prm[0];							// Estimated origin of test image
	float y0 = rot_prm[1];

	int width = img.cols;					// Parameter setup
	int height = img.rows;

	int x_shift = 0;
	int y_shift = 0;

	int count = 0;

	Mat img_shift(height, width, 0);

	// Debugging paramters 
	Mat vals_used(height, width, CV_8UC1);		// Mat used to mark used pixels for debugging
	vals_used = 255;

	cout << "Performing 'simple shift'..." << endl;

	// Loop through every pixel in the image and calculate shift for every other row: 
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {

			// Keeps odd rows
			if (i % 2 + border_pad % 2 == 0) {

				img_shift.at<uchar>(i, j) = img.at<uchar>(i, j);
				vals_used.at<uchar>(i, j) = 0;
				
			}

			else if (i % 2 + border_pad % 2 == 1) {

				x_shift = ((j - x0) * cos(angle) - (i - y0) * sin(angle)) + x0;
				y_shift = ((j - x0) * sin(angle) + (i - y0) * cos(angle)) + y0;
				y_shift = y_shift + (1 - y_shift % 2);

				// Checks if shifted coordinates is within input image, else sets val to 0
				if (x_shift >= 0 and x_shift < width and y_shift >= 0 and y_shift < height) {

					img_shift.at<uchar>(i, j) = img.at<uchar>(y_shift, x_shift);
					vals_used.at<uchar>(y_shift, x_shift) = 0;
					
				}
				else {
					img_shift.at<uchar>(i, j) = 0;
					count++;
				}

			}

		}

	}

	cout << "done. Fraction of unused values: " << float(count)/(width * height) << endl << endl;

	namedWindow("Vals used in rotation", WINDOW_AUTOSIZE);
	imshow("Vals used in rotation", vals_used);

	// Return despinned image
	return img_shift;
}

Mat imShiftCoord(Mat img, array<float, 3> rot_prm, bool lensCorrection) {
	// Used to determine despinned/rotated image coordinates
	// Creates a Mat object with [x,y,val]. For use with interpolation function. 
	// Rotates each image by half the rotation angle to minimize distortion at edges 
	// 
	// Inputs:
	// img:				Interlaced image
	// rot_prm:			Rotation parameters, [x0, y0, theta], same output format as coordinateLSQ
	// lensCorrection:	Enable/disable lens correction (true/false)
	// 

	// Lens correction parameters
	float dx = 8.6;			// Physical pixel size [µm]
	float dy = 8.3;

	float pix_ratio = dy / dx; // Pixel ratio

	int x_l0 = 383 - border_pad;			// Radial distortion center (from Matlab script)
	int y_l0 = 257 - border_pad;

	int efl = 20006;		// Effective focal length
	double kappa = 3.3e-08; // Radial distortion correction 

	double f = efl / dx;	// focal length

	double r = 0;


	// Parameters
	double angle = rot_prm[2];				// Estimated rotation in test image, approximately 1.4 degrees
	angle = angle / 2;

	float x0 = rot_prm[0] - x_l0 * lensCorrection;					// Estimated origin of test image
	float y0 = rot_prm[1] - y_l0 * lensCorrection;

	float xt = 0;							// Temporary coord
	float yt = 0;

	int width = img.cols;					// Parameter setup
	int height = img.rows;

	double x_coord = 0;
	double y_coord = 0;

	Mat img_shift(height, width, CV_64FC3);

	

	// Rotation of the coordinates
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {

			xt = j;
			yt = i;

			// Radial lens correction
			if (lensCorrection) {

				xt -= x_l0;
				yt -= y_l0;

				yt = yt * dy / dx;

				r = sqrt(xt * xt + yt * yt);
				r = 1 + r * r * kappa;

				xt = xt * r;
				yt = yt * r;

				//xt += x_l0;
				//yt += y_l0 * dy / dx;

			}

			if (i % 2 == 1) {

				x_coord = ((xt - x0) * cos(angle) - (yt - y0) * sin(angle)) + x0;
				y_coord = ((xt - x0) * sin(angle) + (yt - y0) * cos(angle)) + y0;

				img_shift.at<Vec3d>(i, j)[0] = x_coord;
				img_shift.at<Vec3d>(i, j)[1] = y_coord;

				img_shift.at<Vec3d>(i, j)[2] = img.at<uchar>(i, j);

			}

			else if (i % 2 == 0) {

				x_coord = ((xt - x0) * cos(angle) + (yt - y0) * sin(angle)) + x0;
				y_coord = (-(xt - x0) * sin(angle) + (yt - y0) * cos(angle)) + y0;

				img_shift.at<Vec3d>(i, j)[0] = x_coord;
				img_shift.at<Vec3d>(i, j)[1] = y_coord;

				img_shift.at<Vec3d>(i, j)[2] = img.at<uchar>(i, j);

			}

		}

	}

	// Returns image with channels [x, y, z] 
	return img_shift;
}

void blobFlux(Mat img, int threshold, int noise_floor) {
	// Finds blobs (stars) and determines their center coordinate and area. Returns list of coordinates (x, y). 
	// Does simple errosion/dilation to remove noise/combine stripes to full stars, 
	// and 4-connectedness to determine position/area of stars. 
	// Used in coordinateLSQ.
	// 
	// Can optionally print coordinates and area of each found blob. 
	// 
	// Change area conditions bellow if errors occour (or add more/less errosion and dilation) 
	// 
	// Input
	// img:		Interpolated image
	// 
	// Output
	// Prints coordinate, flux and area to console
	//

	// Variables
	int min_flux = 300;
	int min_area = 100;

	// Paramter setup
	Mat img_tmp = img > threshold;	// Binary image for blob analysis
	
	int width = img.cols;		// Parameter setup
	int height = img.rows;

	// Remove speckles by erosion
	img_tmp = erodeLR(img_tmp);
	img_tmp = erodeLR(img_tmp);

	
	// Remove left side (since it is mostly noise)
	/*
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < 30; j++) {

			img_tmp.at<uchar>(i, j) = 0;

		}
	}
	*/

	// dilate back to make stars uniform, add extra padding
	img_tmp = dilateCross(img_tmp);
	img_tmp = dilateCross(img_tmp); 
	img_tmp = dilateCross(img_tmp); 
	img_tmp = dilateCross(img_tmp); 


	Mat img_used; 
	img_tmp.copyTo(img_used);

	// More parameters
	int border = 0;
	int area = 0;

	int x = 0;
	int y = 0;

	float x_sum = 0;
	float y_sum = 0;

	int x_min = 0;
	int x_max = 0;

	int y_min = 0;
	int y_max = 0;

	int x_0 = 0;		// Blob start coordinate for connectedness
	int y_0 = 0;

	int x_peak = 0;		// Find maxpoint in blob 
	int y_peak = 0;
	int val_peak = 0;

	int direction = 0;
	bool skip_left = false;

	double flux = 0;

	int val = 0;

	cout << "Object fluxes: " << endl;

	// Scan through binary image 
	for (int i = 1; i < height - 1; i++) {
		for (int j = 1; j < width - 1; j++) {

			// Four connectedness when object is found:
			if (img_tmp.at<uchar>(i, j) == 255) {

				// Min/max x of the four connectedness, to cut out object afterwards
				x_min = j;
				x_max = j;

				y_min = i;
				y_max = i;

				// Loop until back at start
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

					// Find next pixel
					// Directions: 0 -> up, 1 -> right, 2 -> down, 3 -> left
					for (int k = 0; k < 4; k++) {

						skip_left = false;

						// Up
						if (y > 0) {
							val = img_tmp.at<uchar>(y - 1, x);
							if (direction == 0 and val == 255) {

								// Update coordinate and direciton, then break
								y = y - 1;
								direction = 3;
								skip_left = true;  // added so that it doesnt trigger left right away
								//cout << "up" << endl;
								break;
							}
						}
						

						// Right
						if (x < width - 1) {
							val = img_tmp.at<uchar>(y, x + 1);
							if (direction == 1 and val == 255) {

								x = x + 1;
								direction = 0;
								//cout << "right" << endl;
								break;
							}
						}
						

						// Down
						if (y < height - 1) {
							val = img_tmp.at<uchar>(y + 1, x);
							if (direction == 2 and val == 255 and y < height) {

								y = y + 1;
								direction = 1;
								//cout << "down" << endl;
								break;
							}
						}
						

						// Left
						if (x > 0) {
							val = img_tmp.at<uchar>(y, x - 1);
							if (direction == 3 and val == 255 and !skip_left and x > 0) {

								x = x - 1;
								direction = 2;
								//cout << "left" << endl;
								break;
							}
						}
						

						// Rotate start direction 90 degrees
						direction = (direction + 1) % 4;

					}

				}

				
				// Calculate flux and area 
				for (int k = y_min; k < y_max + 1; k++) {
					for (int l = x_min; l < x_max + 1; l++) {

						if (img_tmp.at<uchar>(k, l) == 255) {
							
							val = img.at<uchar>(k, l);
							

							// Compute area and x_sum/y_sum for center value
							if (val > 1 * noise_floor) {
								
								flux += img.at<uchar>(k, l) - noise_floor;
								img_used.at<uchar>(k, l) = 120;

								// Add area
								area += 1;
								x_sum += l;
								y_sum += k;

								// img_used.at<uchar>(k, l) = 120;


								if (val > val_peak) {
									val_peak = val;
									x_peak = l;
									y_peak = k;
								}
							}

							// img_tmp.at<uchar>(k, l) = 0;		// Remove blob from binary image
						}

					}

				}


				if (flux > min_flux and area > min_area) {
					
					//intensityPlot(img, x_peak, y_peak);  // Plot intensity distribution

					cout << "flux: " << flux;
					// cout << ", [x,y]: " << x_sum / area << ", " << y_sum / area;
					cout << ", [x,y]: " << x_peak << ", " << y_peak;
					cout << ", area: " << area;
					cout << endl;

					
				}
				
				// Remove blob from image
				for (int k = y_min; k < y_max + 1; k++) {
					for (int l = x_min; l < x_max + 1; l++) {

						img_tmp.at<uchar>(k, l) = 0;		// Remove blob from binary image

						/*
						if (flux > min_flux and area > min_area) {
							val = img.at<uchar>(k, l);
							if (val > noise_floor) {
								img_used.at<uchar>(k, l) = 120;  // 
							}
						}
						*/

					}

				}

				img_used.at<uchar>(y_peak, x_peak) = 170;

				// Reset parameters
				x_sum = 0;
				y_sum = 0;

				area = 0;
				flux = 0;

				border = 0;

				val_peak = 0;

			}

		}
	}
	
	namedWindow("Image flux vals used", WINDOW_AUTOSIZE);
	imshow("Image flux vals used", img_used);

	imwrite("img_flux_used.png", img_used);

	cout << endl << endl;

	return void();

}

Mat eulerRot(int ax1, int ax2, int ax3, float ang1, float ang2, float ang3, bool useDeg) {
	// Creates rotation matrix based on euler angles
	// Inputs
	// ax: Axis to rotate around 
	// ang: Angle of rotation
	// 
	// Output
	// Rotation matrix as Mat object (float) 
	//
	
	// Parameter setup
	array<float, 3> ax = { ax1, ax2, ax3 };
	array<float, 3> ang = { ang1, ang2, ang3 };

	float c = 0;
	float s = 0;

	// Convert from degress
	if (useDeg) {
		for (int i = 0; i < 3; i++) {
			ang[i] = ang[i] * pi / 180;
		}
	}
	
	// Create empty matrix
	Mat Rot(3, 3, CV_32FC1);
	Mat tmp(3, 3, CV_32FC1);

	Rot = 0;
	tmp = 0;

	// add ones to the diagonal: 
	for (int i = 0; i < 3; i++) {
		Rot.at<float>(i, i) = 1;
	}
	
	// Calculate rotation matrix
	for (int i = 0; i < 3; i++) {
		c = cos(ang[i]);
		s = sin(ang[i]);

		// Calculate individual roation
		if (ax[i] == 1) {
			tmp.at<float>(0, 0) = 1;
			tmp.at<float>(1, 1) = c;
			tmp.at<float>(1, 2) = s;
			tmp.at<float>(2, 1) = -s;
			tmp.at<float>(2, 2) = c;
		}
		else if (ax[i] == 2) {
			tmp.at<float>(0, 0) = c;
			tmp.at<float>(0, 2) = -s;
			tmp.at<float>(1, 1) = 1;
			tmp.at<float>(2, 0) = s;
			tmp.at<float>(2, 2) = c;
		}
		else if (ax[i] == 3) {
			tmp.at<float>(0, 0) = c;
			tmp.at<float>(0, 1) = s;
			tmp.at<float>(1, 0) = -s;
			tmp.at<float>(1, 1) = c;
			tmp.at<float>(2, 2) = 1;
		}

		// Add to Rot matrix and reset tmp
		Rot = tmp * Rot;
		tmp = 0;

	}

	// Check and remove really small values (usually rounding errors
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (abs(Rot.at<float>(i, j)) < 1e-6) {
				Rot.at<float>(i, j) = 0;
			}
		}
	}

	return Rot;
}


pair<float, float> angleCoord(float x, float y) {
	// Calculate 'undistorted theta and rho mapping'
	// Input: Image coordinate 
	// Based on matlab code
	// Not optimal for use with large number of points (since rotation matrix is calculated each time) 


	// parameters
	float dx = 8.6;					// Physical pixel size [µm]
	float dy = 8.3;
	float pix_ratio = dy / dx;		// Pixel ratio
	int x_l0 = 383 - border_pad;	// Lens center (from Matlab script)
	int y_l0 = 257 - border_pad;
	int efl = 20006;				// Effective focal length 
	double kappa = 3.3e-08;			// distortion correction 
	double f = efl / dx;			// focal length
	
	// Convert to virtual CCD coordinates
	float xp = (x - x_l0) / f;
	float yp = (y - y_l0) * pix_ratio / f;

	double r = sqrt(xp * xp + yp * yp + 1);
	double theta = 0;
	double rho = 0;
	
	Mat vec(1, 3, CV_32FC1);
	Mat Rot = eulerRot(1, 2, 3, 167.387, 0.4255, -0.3141, true);		// Spacecraft rotation matrix, values from given script

	// Vector coordinates
	vec.at<float>(0, 0) = xp / r;
	vec.at<float>(0, 1) = yp / r;
	vec.at<float>(0, 2) = 1 / r;

	vec = vec * Rot;

	// Determine theta,rho angles
	theta = atan2(vec.at<float>(0, 1), vec.at<float>(0, 0));
	rho = acos(-vec.at<float>(0, 2));
	
	cout << "[theta, rho]: " << theta * float(180) / pi << ", " << rho * float(180) / pi << endl;

	return make_pair(theta, rho);

}


array<float, 3> coordinateLSQMulti(string im_names, bool lensCorrection) {
	// An expansion of coordinateLSQ to work with multiple images. 
	// Reads specified images from given folder. 
	// Names should be given as: "img01.png,img02.png,img03.png," 
	// The commma at the end is necessary, otherwise it skips the last image. 
	// 
	// Least squares solution to find rotation parameters for image. 
	// Based on non-linear least squares.
	// 
	// Depend on the following functions: 
	// splitBinary
	// blobCoordinates
	//		erodeLR
	//		dilateCross
	// 
	// Inputs 
	// img:				Image, interlaced
	// lensCorrection:  true/false if lens correction should be enabled
	//					currently doesn't work properly
	// 
	// Ouput
	// Roation parameters as array [x, y, theta]
	// 

	// Variables 
	int threshold = 40;
	int dist_max = 30;		// Used to remove objects in only one image if different amount of objects are found

	cout << "MultiLSQ: determining rotation paramters..." << endl;

	Mat img; 

	// Parameter setup
	pair<Mat, Mat> img_binary;// = splitBinary(img, threshold);

	//cout << "Stars in first image: " << endl;
	pair<Mat, Mat> coord_two;// = blobCoordinates(img_binary.first);

	//cout << endl << "Stars in second image: " << endl;
	pair<Mat, Mat> coord_one;// = blobCoordinates(img_binary.second);

	int N = 0;
	int N2 = 0;
	int N_points = 0;

	int tmp = 0;
	float dist = 0;
	float dist_min = 9999;


	array<float, 3> parameters = { 0, 0, 0 };

	Mat coord_one_first_old; 
	Mat coord_one_second_old;
	Mat coord_one_first_new;
	Mat coord_one_second_new;

	Mat coord_two_first_old;
	Mat coord_two_second_old;
	Mat coord_two_first_new;
	Mat coord_two_second_new;

	// Lens correction parameters
	float dx = 8.6;				// Physical pixel size [µm]
	float dy = 8.3;
	float pix_ratio = dy / dx;	// Pixel ratio
	int x_l0 = 383 - border_pad;				// Lens center (from Matlab script)
	int y_l0 = 257 - border_pad;
	int efl = 20006;			// Effective focal length 
	double kappa = 3.3e-08;		// distortion correction 
	double f = efl / dx;		// focal length
	double r = 0;

	string delimiter = ",";  
	string token = ""; 
	string path = "G:/My Drive/UNI/30330 Image Analysis with Microcomputer/30330 Project/pngs_lsq/"; 

	while (im_names.find(delimiter) < 99) {
		token = im_names.substr(0, im_names.find(delimiter));
		im_names.erase(0, im_names.find(delimiter) + delimiter.length());
		cout << "Searcing image: " << token << endl;

		N = 0;
		N2 = 0;

		coord_one_first_new.copyTo(coord_one_first_old);
		coord_one_second_new.copyTo(coord_one_second_old);
		coord_two_first_new.copyTo(coord_two_first_old);
		coord_two_second_new.copyTo(coord_two_second_old);

		// Load image
		img = imread(path + token, IMREAD_GRAYSCALE);
		//img = RemoveSmallParticleNoise(img);

		// Split image 
		img_binary = splitBinary(img, threshold);

		// Find coordinates of objects 
		coord_two = blobCoordinates(img_binary.first, threshold);
		coord_one = blobCoordinates(img_binary.second, threshold);

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

		if (N == 0 or N2 == 0) {
			cout << "Error LSQ: No objects found in image, returning zeroes..." << endl << endl;
			return parameters;
		}

		if (N != N2) {

			cout << "   Warning: different amount of objects found in binary images, " << N << "/" << N2 << ", trying to remove invalid values" << endl;

			// Loop through distance test twice 
			// (In case there are more invalid coords in set two but, but still invalid coords in set 1)
			for (int l = 0; l < 2; l++) {

				// If more objects are found in first image: 
				if (N > N2) {

					// Remove corrdinates for objects with no nearby objects in other image
					for (int i = 0; i < N; i++) {

						dist_min = 9999;

						for (int j = 0; j < N2; j++) {

							// Determine closest coordinate set: 
							dist = pow((coord_one.first.at<float>(i, 0) - coord_two.first.at<float>(j, 0)), 2) + pow((coord_one.second.at<float>(i, 0) - coord_two.second.at<float>(j, 0)), 2);
							dist = sqrt(dist);

							if (dist < dist_min) {
								dist_min = dist;
							}

						}

						if (dist_min > dist_max) {

							for (int k = i; k < N; k++) {
								coord_one.first.at<float>(k, 0) = coord_one.first.at<float>(k + 1, 0);
								coord_one.second.at<float>(k, 0) = coord_one.second.at<float>(k + 1, 0);
							}

						}

					}

					// Calculate number of valid coordinates again
					N = 0;
					for (int i = 0; i < coord_one.first.rows; i++) {
						if (coord_one.first.at<float>(i, 0) < 0) {
							break;
						}
						N++;
					}


				}


				if (N2 > N) {
					// Remove objects with no nearby object in other image
					tmp = 0;
					for (int i = 0; i < N; i++) {

						dist_min = 9999;

						for (int j = 0; j < N2; j++) {

							// Determine closest coordinate set: 
							dist = pow((coord_one.first.at<float>(j, 0) - coord_two.first.at<float>(i, 0)), 2) + pow((coord_one.second.at<float>(j, 0) - coord_two.second.at<float>(i, 0)), 2);
							dist = sqrt(dist);

							if (dist < dist_min) {
								dist_min = dist;
							}

						}

						if (dist_min > dist_max) {

							for (int k = i; k < N + 1; k++) {
								coord_two.first.at<float>(k, 0) = coord_two.first.at<float>(k + 1, 0);
								coord_two.second.at<float>(k, 0) = coord_two.second.at<float>(k + 1, 0);
							}

						}

					}

					// Calculate number of valid coordinates again
					N2 = 0;
					for (int i = 0; i < coord_two.first.rows; i++) {
						if (coord_two.first.at<float>(i, 0) < 0) {
							break;
						}
						N2++;
					}


				}
			}

			if (N != N2) {
				cout << "Error LSQ: Couldn't match objects, returning emtpy array..." << endl << endl;
				return parameters;
			}

		}

		
		// Add new values to mat objects; 
		coord_one_first_new = Mat::zeros(N_points + N, 1, CV_32FC1);
		coord_one_second_new = Mat::zeros(N_points + N, 1, CV_32FC1);
		coord_two_first_new = Mat::zeros(N_points + N, 1, CV_32FC1);
		coord_two_second_new = Mat::zeros(N_points + N, 1, CV_32FC1);

		for (int i = 0; i < N_points; i++) {
			// Move old coords into new set
			coord_one_first_new.at<float>(i, 0) = coord_one_first_old.at<float>(i, 0);
			coord_one_second_new.at<float>(i, 0) = coord_one_second_old.at<float>(i, 0);
			coord_two_first_new.at<float>(i, 0) = coord_two_first_old.at<float>(i, 0);
			coord_two_second_new.at<float>(i, 0) = coord_two_second_old.at<float>(i, 0);
		}
		for (int i = N_points; i < N_points + N; i++) {
			// Move coord from new image into new set
			coord_one_first_new.at<float>(i, 0) = coord_one.first.at<float>(i - N_points, 0);
			coord_one_second_new.at<float>(i, 0) = coord_one.second.at<float>(i - N_points, 0);
			coord_two_first_new.at<float>(i, 0) = coord_two.first.at<float>(i - N_points, 0);
			coord_two_second_new.at<float>(i, 0) = coord_two.second.at<float>(i - N_points, 0);
		}

		N_points += N;

	}


	coord_one.first = coord_one_first_new;
	coord_one.second = coord_one_second_new;
	coord_two.first = coord_two_first_new;
	coord_two.second = coord_two_second_new;
	N = N_points;

	cout << "valid objects found: " << N_points << endl;

	// Initialize least squares variables: 
	Mat x1(N, 1, CV_64FC1);
	Mat y1(N, 1, CV_64FC1);
	Mat x2(N, 1, CV_64FC1);
	Mat y2(N, 1, CV_64FC1);

	Mat d(2 * N, 1, CV_64FC1);

	Mat g0(2 * N, 1, CV_64FC1);

	Mat dg_dt(2 * N, 1, CV_64FC1);
	Mat dg_dx(2 * N, 1, CV_64FC1);
	Mat dg_dy(2 * N, 1, CV_64FC1);

	Mat G(2 * N, 3, CV_64FC1);

	double m_theta = -1.4968 * pi / 180;
	double m_x = 0;// 370;
	double m_y = 0;// -240;

	Mat m(3, 1, CV_64FC1);
	m.at<double>(0, 0) = m_theta;
	m.at<double>(1, 0) = m_x;
	m.at<double>(2, 0) = m_y;


	// Move data into variables 
	for (int i = 0; i < N; i++) {
		x1.at<double>(i, 0) = (coord_one.first.at<float>(i, 0));
		y1.at<double>(i, 0) = (coord_one.second.at<float>(i, 0));

		x2.at<double>(i, 0) = (coord_two.first.at<float>(i, 0));
		y2.at<double>(i, 0) = (coord_two.second.at<float>(i, 0));
	}

	// Apply lens correction to coordinates 
	if (lensCorrection) {

		double xt = 0;
		double yt = 0;

		for (int i = 0; i < N; i++) {
			// Correction for set 1
			xt = x1.at<double>(i, 0) - x_l0;
			yt = y1.at<double>(i, 0) - y_l0;
			yt = yt * dy / dx;

			r = sqrt((xt - x_l0) * (xt - x_l0) + (yt - y_l0) * (yt - y_l0));
			r = 1 + r * r * kappa;

			x1.at<double>(i, 0) = xt * r + x_l0;
			y1.at<double>(i, 0) = yt * r + y_l0 * dy / dx;


			// Correction for set 2
			xt = x2.at<double>(i, 0);
			yt = y2.at<double>(i, 0);

			r = sqrt((xt - x_l0) * (xt - x_l0) + (yt - y_l0) * (yt - y_l0));
			r = 1.0 + r * r * kappa;

			x2.at<double>(i, 0) = xt * r;
			y2.at<double>(i, 0) = yt * r;
		}


	}

	// Setup data vecter
	for (int i = 0; i < 2 * N; i++) {
		if (i < N) {
			d.at<double>(i, 0) = x2.at<double>(i, 0);
		}
		else {
			d.at<double>(i, 0) = y2.at<double>(i - N, 0);
		}

	}

	float step = 0.1;
	// Non linear least squares computation
	for (int j = 0; j < 100; j++) {

		// Calculate g0
		for (int i = 0; i < 2 * N; i++) {
			if (i < N) {
				g0.at<double>(i, 0) = (x1.at<double>(i, 0) - m_x) * cos(m_theta) - (y1.at<double>(i, 0) - m_y) * sin(m_theta) + m_x;
			}
			else {
				g0.at<double>(i, 0) = (x1.at<double>(i - N, 0) - m_x) * sin(m_theta) + (y1.at<double>(i - N, 0) - m_y) * cos(m_theta) + m_y;
			}
		}

		// Calculate dg_dtheta, dg_dx0 and dg_dy0
		for (int i = 0; i < 2 * N; i++) {
			if (i < N) {
				dg_dt.at<double>(i, 0) = -(x1.at<double>(i, 0) - m_x) * sin(m_theta) - (y1.at<double>(i, 0) - m_y) * cos(m_theta);
				dg_dx.at<double>(i, 0) = 1 - cos(m_theta);
				dg_dy.at<double>(i, 0) = sin(m_theta);
			}
			else {
				dg_dt.at<double>(i, 0) = (x1.at<double>(i - N, 0) - m_x) * cos(m_theta) - (y1.at<double>(i - N, 0) - m_y) * sin(m_theta);
				dg_dx.at<double>(i, 0) = -sin(m_theta);
				dg_dy.at<double>(i, 0) = 1 - cos(m_theta);
			}
		}

		// Values into G
		G.col(0) = (dg_dt + 0);
		G.col(1) = (dg_dx + 0);
		G.col(2) = (dg_dy + 0);

		m = m + step * ((G.t() * G).inv()) * G.t() * (d - g0);

		// cout << "G.Gt() inv: " << endl << ((G * G.t()).inv()).size() << endl;

		
		if (abs(m_x - m.at<double>(1, 0)) / (m.at<double>(1, 0)) < 1e-9) {

			m_theta = m.at<double>(0, 0);
			m_x = m.at<double>(1, 0);
			m_y = m.at<double>(2, 0);

			cout << "broke at itteration: " << j << endl;

			break;
		}
		

		m_theta = m.at<double>(0, 0);
		m_x = m.at<double>(1, 0);
		m_y = m.at<double>(2, 0);

	}

	parameters[0] = m.at<double>(1, 0);
	parameters[1] = m.at<double>(2, 0);
	parameters[2] = m.at<double>(0, 0);


	// Print parameters to console
	cout << "rotation parameters: (" << parameters[0] << ", " << parameters[1] << ", " << parameters[2] * 180 / pi << ") " << endl;
	cout << endl;

	return parameters;
	 
}
