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
Mat imShiftCoord(Mat img, array<float, 3> rot_prm, bool lensCorrection = false);     // Creates a 3d array with [x,y,val] along the 3rd axis

// Other functions
pair<Mat, Mat> blobCoordinates(Mat img);
array<float, 3> coordinateLSQ(Mat img, bool lensCorrection = false);

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
			if (i % 2 == 0) {
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

pair<Mat, Mat> blobCoordinates(Mat img) {
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

	// Paramter setup
	Mat img_tmp;
	img.copyTo(img_tmp);		// Temporary image for analysis

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

	// More parameters
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

	// Matrices to hold coordinates. Max set at 20 since it's probably way more than needed. 
	Mat x_pos(20, 1, CV_32FC1);
	Mat y_pos(20, 1, CV_32FC1); 

	int count = 0;		// Used to count number of stars found

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

						// Rotate start direction 90 degrees
						direction = (direction + 1) % 4;

					}

				}

				// Remove blob from image
				for (int k = y_min; k < y_max + 1; k++) {
					for (int l = x_min; l < x_max + 1; l++) {

						val = img_tmp.at<uchar>(k, l);

						// Compute area and x_sum/y_sum for center value
						if (val == 255) {
							// Add area
							area += 1;
							x_sum += l;
							y_sum += k;

							img_tmp.at<uchar>(k, l) = 0;
						}

					}

				}

				// Save coordinate if area of object is within conditions
				if (area > 30 and area < 150 and count < 20) {

					x_pos.at<float>(count, 0) = float(x_sum) / area;
					y_pos.at<float>(count, 0) = float(y_sum) / area;

					//cout << "" << int(x_pos.at<float>(count, 0)) << ", " << int(y_pos.at<float>(count, 0)) << ", area: " << area << endl;

					count++;
					if (count == 20) {
						cout << "Too many stars, using first 20... " << endl;
					}
				}
				
				// Reset parameters
				x_sum = 0;
				y_sum = 0;
				area = 0;
				border = 0;

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
	int threshold = 230;

	// Parameter setup
	pair<Mat, Mat> img_binary = splitBinary(img, threshold); 

	//cout << "Stars in first image: " << endl;
	pair<Mat, Mat> coord_two = blobCoordinates(img_binary.first);

	//cout << endl << "Stars in second image: " << endl;
	pair<Mat, Mat> coord_one = blobCoordinates(img_binary.second);

	int N = 0;
	int N2 = 0;

	array<float, 3> parameters = { 0, 0, 0 };
	
	// Lens correction parameters
	float dx = 8.6;				// Physical pixel size [µm]
	float dy = 8.6;
	float pix_ratio = dy / dx;	// Pixel ratio
	int x_l0 = 383;				// Lens center (from Matlab script)
	int y_l0 = 257;
	int efl = 20006;			// Focal length also?  
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

	if (N != N2) {
		cout << "Error LSQ: different amount of objects found in binary images, returning empty array... " << endl; 
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
	float angle = rot_prm[2];

	float x0 = rot_prm[0];							// Estimated origin of test image
	float y0 = rot_prm[1];

	int width = img.cols;					// Parameter setup
	int height = img.rows;

	int x_shift = 0;
	int y_shift = 0;

	Mat img_shift(height, width, 0);

	// Debugging paramters 
	Mat vals_used(height, width, CV_8UC1);		// Mat used to mark used pixels for debugging
	vals_used = 255;

	// Loop through every pixel in the image and calculate shift for every other row: 
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {

			// Keeps odd rows
			if (i % 2 == 1) {

				img_shift.at<uchar>(i, j) = img.at<uchar>(i, j);
				vals_used.at<uchar>(i, j) = 0;
			}

			else if (i % 2 == 0) {

				x_shift = ((j - x0) * cos(angle) - (i - y0) * sin(angle)) + x0;
				y_shift = ((j - x0) * sin(angle) + (i - y0) * cos(angle)) + y0;
				y_shift = y_shift + y_shift % 2;

				// Checks if shifted coordinates is within input image, else sets val to 0
				if (x_shift >= 0 and x_shift < width and y_shift >= 0 and y_shift < height) {

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

	// Parameters
	double angle = rot_prm[2];   // Estimated rotation in test image, approximately 1.4 degrees
	angle = angle / 2;

	float x0 = rot_prm[0];							// Estimated origin of test image
	float y0 = rot_prm[1];

	float xt = 0;			// Temporary coord
	float yt = 0;

	int width = img.cols;					// Parameter setup
	int height = img.rows;

	double x_coord = 0;
	double y_coord = 0;

	Mat img_shift(height, width, CV_64FC3);

	// Lens correction parameters
	float dx = 8.6;			// Physical pixel size [µm]
	float dy = 8.6;

	float pix_ratio = dy / dx; // Pixel ratio

	int x_l0 = 383;			// Radial distortion center (from Matlab script)
	int y_l0 = 257;

	int efl = 20006;		// Effective focal length
	double kappa = 3.3e-08; // Radial distortion correction 

	double f = efl / dx;	// focal length

	double r = 0;


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

				xt += x_l0;
				yt += y_l0 * dy / dx;

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
