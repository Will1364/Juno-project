#include "header.h"

// 
// Functions for figures.  
// Text for header file: 
/*

// Figures
Mat imBrightened(Mat img, float scale);
Mat starPlot(Mat img, int x, int y, int scale);
void intensityPlot(Mat img, int x_peak, int y_peak);

*/

Mat imBrightened(Mat img, float scale) {
	// Multiplies whole image with scale
	// Parameters
	Mat img_scaled;
	img.copyTo(img_scaled);

	int width = img.cols;		// Parameter setup
	int height = img.rows;
	float val = 0;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {

			val = float(img.at<float>(i, j)) * scale;

			if (val < 1.0) {
				img_scaled.at<float>(i, j) = val;
			}
			else {
				img_scaled.at<float>(i, j) = 1.0;
			}

		}
	}

	return img_scaled;
}

void intensityPlot(Mat img, int x_peak, int y_peak) {
	// Plots the intensity in x and y direction for a specified point (x_peak, y_peak). 
	// Restricts size to 
	// 

	// Variables 
	int height = 260;		// Size of plot
	int width = 260 * 2;

	int noiseFloor = 20;	// Noise floor -> used to determine 'end' of object 
	int padding = 5;				// Number of points beneath noise floor before ending 
	int max_dist = 20;				// Maximum  distance from peak (-> maximum 40 point wide figure)

	// Parameters
	Mat imgPlot(height, width, CV_8UC1);
	imgPlot = 0;

	int x_min = 0;
	int y_min = 0;
	int x_max = 0;
	int y_max = 0;

	float val = 0;
	int tmp = 0;

	// Drawing parameters
	int pointWidth = 0;
	int y0 = 0;
	int y1 = 0;
	Point pt1, pt2;

	int arr_size_x = 0;
	int arr_size_y = 0;
	int arr_width = 0;

	// Find x_min
	for (int j = x_peak; tmp < padding + 1; j--) {
		val = img.at<float>(y_peak, j) * 255;

		if (val > noiseFloor and tmp > 0) {
			tmp--;
		}
		else if (val < noiseFloor) {
			tmp++;
			if (tmp == padding) {
				x_min = j;
			}
		}
		if (x_peak - j > max_dist or j == 0) {
			x_min = j;
			break;
		}
	}
	tmp = 0;

	// Find x_max 
	for (int j = x_peak; tmp < padding + 1; j++) {

		val = img.at<float>(y_peak, j) * 255;

		if (val > noiseFloor and tmp > 0) {
			tmp--;
		}
		else if (val < noiseFloor) {
			tmp++;
			if (tmp == padding) {
				x_max = j;
			}
		}
		if (j - x_peak > max_dist or j == img.cols) {
			x_max = j;
			break;
		}
	}
	tmp = 0;


	// Find y_min
	for (int i = y_peak; tmp < padding + 1; i--) {
		val = img.at<float>(i, x_peak) * 255;

		if (val > noiseFloor and tmp > 0) {
			tmp--;
		}
		else if (val < noiseFloor) {
			tmp++;
			if (tmp == padding) {
				y_min = i;
			}
		}
		if (y_peak - i > max_dist or i == 0) {
			y_min = i;
			break;
		}
	}
	tmp = 0;

	// Find y_max
	for (int i = y_peak; tmp < padding + 1; i++) {
		val = img.at<float>(i, x_peak) * 255;

		if (val > noiseFloor and tmp > 0) {
			tmp--;
		}
		else if (val < noiseFloor) {
			tmp++;
			if (tmp == padding) {
				y_max = i;
			}
		}
		if (i - y_peak > max_dist or i == img.rows) {
			y_max = i;
			break;
		}
	}
	tmp = 0;


	arr_size_x = x_max - x_min + 1;
	arr_size_y = y_max - y_min + 1;
	arr_width = max(x_peak - x_min, y_peak - y_min) + max(x_max - x_peak, y_max - y_peak);	// Combined width 

	vector<int> x_arr(arr_size_x, 0);
	vector<int> y_arr(arr_size_y, 0);

	// Fill x-array
	for (int j = x_min; j <= x_max; j++) {
		x_arr[j - x_min] = img.at<float>(y_peak, j) * 255;
	}

	// Fill y-array
	for (int i = y_min; i <= y_max; i++) {
		y_arr[i - y_min] = img.at<float>(i, x_peak) * 255;
	}

	// Width of points
	pointWidth = (width) / arr_width;


	// Draw x-grid 
	for (int i = 0; i <= arr_width; i++) {
		pt1 = Point2d((i)*pointWidth, 0);
		pt2 = Point2d((i)*pointWidth, height);

		line(imgPlot, pt1, pt2, 40);
	}
	//brighter line at peak position
	tmp = max(x_peak - x_min, y_peak - y_min);
	pt1 = Point2d((tmp)*pointWidth, 0);
	pt2 = Point2d((tmp)*pointWidth, height);
	line(imgPlot, pt1, pt2, 60);

	// Draw noise floor
	//pt1 = Point2d(0, 260 - noiseFloor);
	//pt2 = Point2d(width, 260 - noiseFloor);
	//line(imgPlot, pt1, pt2, 40);

	// Draw max line
	pt1 = Point2d(0, height - 255);
	pt2 = Point2d(width, height - 255);
	line(imgPlot, pt1, pt2, 60);



	// Plot y-direction
	tmp = max(x_peak - x_min, y_peak - y_min) - (y_peak - y_min);
	for (int i = 0; i < arr_size_y; i++) {
		if (i == 0) {
			y0 = height;
			y1 = height - y_arr[i];

			pt1 = Point2d((i + tmp + 1) * pointWidth, y0);
			pt2 = Point2d((i + tmp + 1) * pointWidth, y1);
		}
		else if (i == arr_size_y - 1) {
			y0 = height;
			y1 = height - y_arr[i];

			pt1 = Point2d((i + tmp) * pointWidth, y0);
			pt2 = Point2d((i + tmp) * pointWidth, y1);
		}
		else {
			y0 = height - y_arr[i];
			y1 = height - y_arr[i + 1];

			pt1 = Point2d((i + tmp) * pointWidth, y0);
			pt2 = Point2d((i + tmp + 1) * pointWidth, y1);
		}


		line(imgPlot, pt1, pt2, 200);

	}

	// Plot x-direction
	tmp = max(x_peak - x_min, y_peak - y_min) - (x_peak - x_min);
	for (int j = 0; j < arr_size_x; j++) {
		if (j == 0) {
			y0 = height;
			y1 = height - x_arr[j];

			pt1 = Point2d((j + tmp + 1) * pointWidth, y0);
			pt2 = Point2d((j + tmp + 1) * pointWidth, y1);
		}
		else if (j == arr_size_x - 1) {
			y0 = height;
			y1 = height - x_arr[j];

			pt1 = Point2d((j + tmp) * pointWidth, y0);
			pt2 = Point2d((j + tmp) * pointWidth, y1);
		}
		else {
			y0 = height - x_arr[j];
			y1 = height - x_arr[j + 1];

			pt1 = Point2d((j + tmp) * pointWidth, y0);
			pt2 = Point2d((j + tmp + 1) * pointWidth, y1);
		}



		line(imgPlot, pt1, pt2, 100);

	}


	string windowName = "Intensity plot for object at [" + to_string(x_peak) + ", " + to_string(y_peak) + "]";
	namedWindow(windowName, WINDOW_AUTOSIZE);
	imshow(windowName, imgPlot);
	imwrite("intensityPlot.png", imgPlot);

	return void();

}

Mat starPlot(Mat img, int x, int y, int scale) {
	// Variables 
	int height = 260 + (scale - 260 % scale);		// Size of plot, adjusted to be divisble by scale
	int width = 260 + (scale - 260 % scale);

	int pix_width = width / scale;

	int x_c = x - pix_width / 2;
	int y_c = y - pix_width / 2;

	Mat img_scaled(height, width, CV_32FC1);
	img_scaled = 0;

	for (int i = 0; i < pix_width; i++) {
		for (int j = 0; j < pix_width; j++) {

			for (int m = 0; m < scale; m++) {
				for (int n = 0; n < scale; n++) {
					img_scaled.at<float>(i * scale + m, j * scale + n) = img.at<float>(y_c + i, x_c + j);
				}
			}


		}
	}

	string windowName = "Star at [" + to_string(x) + ", " + to_string(y) + "]";
	namedWindow(windowName, WINDOW_AUTOSIZE);
	imshow(windowName, img_scaled);

	return img_scaled;
}
