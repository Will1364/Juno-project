#include "header.h"

// 
// Functions for interpolating images 
// InterpolateNN does the interpolation, based on chosen interpolation 'core' 
// Each core is a different function for determining the interpolated value. 
// 
// Text for header file: 
/*

// Interpolation
Mat interpolateNN(Mat img, int arr_size);

// Interpolation cores
double intpCoreNNPlus(vector<float> min_d, vector<float> min_v, int max_dist);
double intpCoreDist(vector<float> min_d, vector<float> min_v);
double intpCoreRadialBasis(int x_, int y_, vector<float> min_d, vector<int> min_i, vector<int> min_j, vector<float> min_v, Mat img);

*/
//


Mat interpolateNN(Mat img, int arr_size) {
	// Finds the nearest n (arr_size) neighbours for each point, and runs through an interpolation core, 
	// a seperate function for finding an 'optimal' value for the pixel. 
	// Interpolation 'core' (type of interpolation), can be changed further down. 
	// 
	// Inputs
	// img: Image with 3 channels, [x, y, val], as doubles.  
	// corresponds to output of 'imShiftCoord' function
	// arr_size: Number of neighbours to interpolate over
	// 
	// Output:
	// Interpolated image, of original size minus padding in each direction
	// has 8 bit unsigned int as data type
	// 

	// Variables : 
	int search = 50;		// NxN search area around each pixel. 50 is a good balance between compute time and extra buffer 
	int padding = 0;		// Number of pixels to remove from the border. 
	float max_dist = 2;		// The furthest away 'closest' pixel can be 
	float min_dist = 0.05;  // If min distance is closer than this, uses uninterpolated value. Set to zero or negative to disable

	// Parameter setup
	int width = img.cols;
	int height = img.rows;

	int width_i = width - padding;
	int height_i = height - padding;

	int s_x = 0;
	int s_y = 0;

	double dist = 0;

	vector<float> min_d(arr_size, 0);
	vector<float> min_v(arr_size, 0);

	vector<int> min_i(arr_size, 0);
	vector<int> min_j(arr_size, 0);

	double val = 0;
	double tmp = 0;

	Mat interpolated(height_i, width_i, CV_8UC1);

	// Bellow are 'niceness' parameters, used for debugging / to track progress
	Mat vals_used(height, width, CV_8UC1);		// Mat used to display pixels used. Can be commented out without affecting the interpolation
	vals_used = 250;

	int used = 0;		// Counts used pixels from even rows
	int unused = 0;		// Counts number of unused pixels
	int progress = 1;	// Used to display progress
	int count_mindist = 0;


	// Interpolation function
	cout << "interpolating";
	// Loops through all image values
	for (int i = 0; i < height_i; i++) {
		for (int j = 0; j < width_i; j++) {

			// Reset paramters 
			for (int k = 0; k < arr_size; k++) {
				min_d[k] = 9990;
				min_i[k] = 0;
				min_j[k] = 0;
				min_v[k] = 0;
			}

			// Search center coordinate (in input image)  
			s_y = i + padding - (search / 2);
			s_x = j + padding - (search / 2);

			// Fills search window with distance to point in a search x search grid around 
			// Finds nearest arr_size values within search window 
			for (int k = 0; k < search; k++) {
				for (int l = 0; l < search; l++) {

					// Test if search coordinate is within input image 
					if (s_y + k > 0 and s_y + k < height and s_x + l > 0 and s_x + l < width) {

						// Calculate distance
						dist = sqrt(pow(img.at<Vec3d>(s_y + k, s_x + l)[0] - j, 2) + pow(img.at<Vec3d>(s_y + k, s_x + l)[1] - i, 2));


						// Update arrays if with closest values 
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

			
			// Count to check that even and uneven rows from intput image are used roughly equally
			if (min_i[0] % 2 == 0) {
				used += 1;
			}
			
			// Uses uninterpolated value if pixel is closer than min_dist
			if (min_d[0] < min_dist) {
				interpolated.at<uchar>(i, j) = img.at<Vec3d>(min_i[0], min_j[0])[2];
				count_mindist++;
			}

			else if (min_d[0] < max_dist) {

				// Save values for closest pixels
				for (int k = 0; k < arr_size; k++) {
					min_v[k] = img.at<Vec3d>(min_i[k], min_j[k])[2];
				}

				// 
				// Interpolation cores: Seperate functions where interpolation calculations are done
				// Currently made: 
				// Nearest neighbour: Just closest value, val = min_v[0]
				// Inverse distance weighting: intpCoreDist(min_d, min_v), weight by 1/(d**p), with p being weighing parameter
				// Radial basis interpolation: intpCoreRadialBasis(j, i, min_d, min_i, min_j, min_v, img), see function description
				// 

				// val = min_v[0];															// Nearest neighbour
				// val = intpCoreNNPlus(min_d, min_v, max_dist);							// NN+; Takes average of all value within radius of max_dist (1, set in core).  					
				val = intpCoreDist(min_d, min_v);										// Inverse distance weighting 
				// val = intpCoreRadialBasis(j, i, min_d, min_i, min_j, min_v, img);		// Radial basis interpolation


				// Put value into interpolated image 
				interpolated.at<uchar>(i, j) = int(val);

				// Count pixels used (optional, can be commented out)
				vals_used.at<uchar>(min_i[0], min_j[0]) -= 50;

			}
			else {
				// If closest pixel is further away 
				interpolated.at<uchar>(i, j) = 0;
				unused++;
			}

			// Progress bar, can be commented out
			if (i > progress * height / 10) {
				progress++;
				cout << ".";
			}

			// Reset parameters
			val = 0;
			tmp = 0;

		}
	}

	// display val_used image
	//namedWindow("Image intp shift vals used", WINDOW_AUTOSIZE);
	//imshow("Image intp shift vals used", vals_used);
	//imwrite("img_intp_used.png", vals_used);

	// Output text at function end 
	cout << "\rinterpolated                       " << endl;
	cout << "fraction used from each segment: " << float(used) / (height_i * width_i) << endl;
	cout << "fraction unused: " << float(unused) / (height_i * width_i) << endl << endl;
	cout << "Number of pixels closer than min_dist (" << min_dist << "): " << count_mindist << endl;
	cout << endl;

	return interpolated;
}

double intpCoreDist(vector<float> min_d, vector<float> min_v) {
	// Inverse distance weight for interpolation
	// Weights each value with 1/(d**p), d being distance to pixel, p being power 
	// Larger p means pixel further away count less. 

	// Variables
	float max_dist = 10;	// Maximum distance for pixel to be allowed
	float r0 = 0.2;			// Distance scaling, to prevet really small values
	int p = 1.2;			// Power for distance weight. Closer to one -> distance matters less 

	double W = 0;
	double sumW = 0;		// Sum of weights, 1/(d**p)
	double sumWu = 0;		// Sum of val * weights
	

	for (int k = 0; k < min_d.size(); k++) {
		if (min_d[k] < max_dist) {
			W = 1 / pow(min_d[k]/r0, p);
			sumW += W;
			sumWu += min_v[k] * W;
		}
	}
	return sumWu / sumW;
}

double intpCoreRadialBasis(int x_, int y_, vector<float> min_d, vector<int> min_i, vector<int> min_j, vector<float> min_v, Mat img) {
	// Local Radial basis interpolation based on distance function. 
	// Based on chpt. 3.7.1 in https://numerical.recipes/book.html
	// Variant that uses nearest N values for each image. 
	// Based on matrix equation: phi @ w = f
	// 
	// Basis functions that are implemented bellow (comment out the ones not used): 
	// Gaussian: exp(-(r/r0)**2)
	// Inverse quadratic: 1/(1+(r/r0)**2)
	// Thin plate spline: r**2 log(r / r0)     !Doesn't work
	// 


	// Variables: 
	float r0 = 2.0;		// Shape parameter -> used to scale the basis function  

	// Parameter setup
	int N = min_d.size();
	double r;

	Mat f(N, 1, CV_64FC1);		// Function values
	Mat phi(N, N, CV_64FC1);	// Radial function, with phi_ij = b(p_i,p_j), b being some function and p being points

	// 
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {

			r = sqrt(pow(img.at<Vec3d>(min_i[i], min_j[i])[0] - img.at<Vec3d>(min_i[j], min_j[j])[0], 2) + pow(img.at<Vec3d>(min_i[i], min_j[i])[1] - img.at<Vec3d>(min_i[j], min_j[j])[1], 2));

			phi.at<double>(i, j) = exp(-pow(r / r0, 2));		// Gaussian 
			// phi.at<double>(i, j) = pow(r, 2) * log(r / r0);				// Thin plate spline
			// phi.at<double>(i, j) = 1 / (1 + pow( r / r0, 2));			// Inverse quadratic

		}

		f.at<double>(i, 0) = min_v[i];

	}

	Mat w = phi.inv() * f;		// Determine weights
	double s = 0;				// Intp value

	// Calculate interpolation value based on weights
	for (int i = 0; i < N; i++) {
		r = min_d[i];

		s += w.at<double>(i, 0) * exp(-pow(r / r0, 2));					// Gaussian
		// s += w.at<double>(i, 0) * pow(r, 2) * log(r / r0);			// Thin plate spline
		// s += w.at<double>(i, 0) * 1 / (1 + pow(r / r0, 2));			// Inverse quadratic		
	}

	// If statements used to construct values since the function is a bit unstable. 
	// This is probably due to the points being way too close when the lines 'cross' 
	if (s > 0 and s < 256) {
		return s;
	}
	else if (s > 255) {
		return 255;
	}
	return 0;
}

double intpCoreNNPlus(vector<float> min_d, vector<float> min_v, int max_dist) {
	// Extended Nearest neightbours: Takes average of all pixels within max_dist 
	// For max_dist > 1, it acts as an averaging filter, consider using 
	// inverse distance weighting (intpCoreDist) instead. 
	// 

	// Parameters
	double val = 0;
	int count = 0;
	int size = min_d.size();
	// max_dist = 1;	// Overwrite max_dist

	// Loop through array, count and add up values for d < max_dist
	for (int i = 0; i < size; i++) {
		if (min_d[i] < max_dist) {
			val += min_v[i];
			count++;
		}
	}


	return val/count;
}
