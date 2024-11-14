#pragma once
#include <iostream>
#include <opencv2/opencv.hpp>
//#include <utility>  // For std::pair

using namespace std;		// Tells the program that the std library doesn't need an identifier (Skip using 'std::xxx')
using namespace cv;


// Image reading 
Mat imReadGS(string img_name);

// Histogram functions
array<int, 256> histArrayGS(Mat img);
void histDisplay(array<int, 256> histogram, string windowName = "Histogram");

// Split image
pair<Mat, Mat> splitImage(const Mat& OGimage);
pair<Mat, Mat> splitBinary(Mat img);

// Dilation and errosion
Mat dilateCross(Mat img);
Mat erodeLR(Mat img);

// Image shifting
Mat imShiftSimple(Mat img, array<float, 3> rot_prm);	
Mat imShiftCoord(Mat img, array<float, 3> rot_prm);     // Creates a 3d array with [x,y,val] along the 3rd axis


// Interpolations 
Mat interpolateNearestN(Mat img, int arr_size);
// Interpolation core
double intpCoreDist(vector<float> min_d, vector<float> min_v);
// double intpCoreRadialBasis(vector<float> min_d, vector<float> min_v); 
double intpCoreRadialBasis(int x_, int y_, vector<float> min_d, vector<int> min_i, vector<int> min_j, vector<float> min_v, Mat img);

// Other functions
pair<Mat, Mat> blobCoordinates(Mat img);
array<float, 3> coordinateLSQ(Mat img);