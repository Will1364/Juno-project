#include "HeaderFile.h"

Mat RemoveSmallParticleNoise(Mat image) {
	
	//Make a copy of the existing image to compare to the original
	Mat image2 = image.clone();
	copyMakeBorder(image2, image2, 2, 2, 1, 1, BORDER_CONSTANT, Scalar(0, 0, 0));

	//Set some threshold values
	int t1 = 30; //Sets the difference between a pixel and its imidiate neighbour - 35 works best for the unscaled images
	int t2 = 2; //Sets the difference between the 2 neighbours of a pixel

	int width = image2.cols;
	int height = image2.rows;

	//The loop runs through every pixel in the image and compares different scenarios for noise
	//For each scenario there is a sub-scenario for being at the edge of the image
	//If the conditions are met, the particle is assigned the value of the average of the surrounding pixels
	//Since a lot of the noise is made of streched lines rather than single pixels it is prefered to compare the pixels values to only some of the surrounding 8 pixels rather than all of them
	//It is still compared to enough neighbours to make sure the stars are not removed
	for (int i = 2; i < height - 2; i++) {
		for (int j = 1; j < width - 1; j++) {
			//A separate scenario is made for noise followed by a dark echo - here both pixels should be considered simultaneously and both pixels are assigned the value of the average of the surrounding pixels
			if ((image2.at<uchar>(i, j) - image2.at<uchar>(i, j - 1)) >= t1 && (image2.at<uchar>(i, j - 1) - image2.at<uchar>(i, j + 1) >= t2)) {
				if (i == 2) {
					if (j == 1) { image2.at<uchar>(i, j) = (image2.at<uchar>(i + 2, j) + image2.at<uchar>(i + 2, j + 1)) / 2; }
					else image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i + 2, j - 1) + image2.at<uchar>(i + 2, j) + image2.at<uchar>(i + 2, j + 1)) / 4;

					image2.at<uchar>(i, j + 1) = image2.at<uchar>(i, j);
				}
				if (i == height - 2) {
					if (j == 1) { image2.at<uchar>(i, j) = (image2.at<uchar>(i - 2, j) + image2.at<uchar>(i - 2, j + 1)) / 2; }
					else image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i - 2, j - 1) + image2.at<uchar>(i - 2, j) + image2.at<uchar>(i - 2, j + 1)) / 4;

					image2.at<uchar>(i, j + 1) = image2.at<uchar>(i, j);
				}
				else {
					if (j == 1) { image2.at<uchar>(i, j) = (image2.at<uchar>(i - 2, j) + image2.at<uchar>(i + 2, j) + image2.at<uchar>(i + 2, j + 1) + image2.at<uchar>(i - 2, j + 1)) / 4; }
					else image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i-2, j + 1) + image2.at<uchar>(i - 2,j))/ 3;// + image2.at<uchar>(i - 2, j - 1) + image2.at<uchar>(i + 2, j - 1) + image2.at<uchar>(i - 2, j) + image2.at<uchar>(i + 2, j) + image2.at<uchar>(i + 2, j + 1) + image2.at<uchar>(i + 2, j - 1)) / 7;

					image2.at<uchar>(i, j + 1) = image2.at<uchar>(i, j);
				}
			}

			if (i == 2 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j - 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i + 2, j)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i + 2, j - 1)) >= t1 && (image2.at<uchar>(i, j - 1) - image2.at<uchar>(i, j + 1) <= t2)) {
				if (j == 1) {
					image2.at<uchar>(i, j) = (image2.at<uchar>(i, j + 1) + image2.at<uchar>(i + 2, j) + image2.at<uchar>(i + 2, j + 1)) / 3;
				}
				else image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i + 2, j - 1) + image2.at<uchar>(i + 2, j)) / 3;
			}
			else if (i == height - 2 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j - 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i - 2, j)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i - 2, j - 1)) >= t1 && (image2.at<uchar>(i, j - 1) - image2.at<uchar>(i, j + 1) <= t2)) {
				if (j == 1) {
					image2.at<uchar>(i, j) = (image2.at<uchar>(i, j + 1) + image2.at<uchar>(i - 2, j) + image2.at<uchar>(i - 2, j + 1)) / 3;
				}
				else image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i - 2, j - 1) + image2.at<uchar>(i - 2, j)) / 3;
			}
			//else if ((image2.at<uchar>(i, j) - image2.at<uchar>(i, j - 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i + 2, j)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i + 2, j - 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i - 2, j)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i - 2, j - 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i - 2, j + 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i + 2, j + 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i - 2, j + 2)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i + 2, j +2)) >= t1 && (image2.at<uchar>(i, j - 1) - image2.at<uchar>(i, j + 1) <= t2)) { //her
				//if (j == 1) {
					//image2.at<uchar>(i, j) = (image2.at<uchar>(i - 2, j) + image2.at<uchar>(i + 2, j) + image2.at<uchar>(i + 2, j + 1) + image2.at<uchar>(i - 2, j + 1) + image2.at<uchar>(i, j + 1)) / 5;
				//}
				//else image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i - 2, j - 1) + image2.at<uchar>(i + 2, j - 1) + image2.at<uchar>(i - 2, j) + image2.at<uchar>(i + 2, j)) / 5;
			//}

			else if (j <= width - 2 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j - 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i + 2, j + 2)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j + 2)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i - 2, j)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i + 2, j - 1)) >= t1 && (image2.at<uchar>(i, j - 1) - image2.at<uchar>(i, j + 1) <= t2)) {
				if (j == 1) {
					image2.at<uchar>(i, j) = (image2.at<uchar>(i + 2, j + 2) + image2.at<uchar>(i, j + 2) + image2.at<uchar>(i - 2, j)) / 3;
				}
				image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i + 2, j + 2) + image2.at<uchar>(i, j + 2) + image2.at<uchar>(i - 2, j)) / 4;
			}
			else if (j<=width-2 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j + 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i + 2, j)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i + 2, j + 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i - 2, j)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i - 2, j + 1)) >= t1 && (image2.at<uchar>(i, j - 1) - image2.at<uchar>(i, j + 1) <= t2)) {
				image2.at<uchar>(i, j) = (image2.at<uchar>(i, j + 1) + image2.at<uchar>(i + 2, j) + image2.at<uchar>(i + 2, j + 1) + image2.at<uchar>(i - 2, j) + image2.at<uchar>(i - 2, j + 1)) / 5;
			}

			else if ((image2.at<uchar>(i, j) - image2.at<uchar>(i, j - 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j + 1) >= t1 && (image2.at<uchar>(i, j - 1) - image2.at<uchar>(i, j + 1) <= t2))) {
				if (j == 1) {
					image2.at<uchar>(i, j) = (image2.at<uchar>(i, j + 1));
				}
				image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i, j + 1) / 2);
			}
			else if (j <= width - 2 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j - 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j + 2) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i-2,j-1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i - 2, j + 1)) >= t1 && (image2.at<uchar>(i, j - 1) - image2.at<uchar>(i, j + 1) <= t2))) {
				if (j == 1) {
					image2.at<uchar>(i, j) = (image2.at<uchar>(i, j + 2));
				}
				image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i, j + 2) / 2);
			}
			
			//else if (j <= width - 2 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j - 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i + 2, j + 2)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i - 2, j)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i - 2, j + 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i-2, j - 1)) >= t1 && (image2.at<uchar>(i, j - 1) - image2.at<uchar>(i, j + 1) <= t2)) { //her
				//if (j == 1) {
					//image2.at<uchar>(i, j) = (image2.at<uchar>(i - 2, j) + image2.at<uchar>(i + 2, j) + image2.at<uchar>(i + 2, j + 1) + image2.at<uchar>(i - 2, j + 1) + image2.at<uchar>(i, j + 1)) / 5;
				//}
				//else image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i - 2, j) + image2.at<uchar>(i + 2, j + 1) + image2.at<uchar>(i - 2, j + 1)) / 4;
			//}
			else if (j <= width - 2 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j - 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j +2)) >= t1 && (image2.at<uchar>(i, j - 1) - image2.at<uchar>(i, j + 1) <= t2)) { //her
				if (j == 1) {
					image2.at<uchar>(i, j) = (image2.at<uchar>(i - 2, j) + image2.at<uchar>(i + 2, j) + image2.at<uchar>(i + 2, j + 1) + image2.at<uchar>(i - 2, j + 1) + image2.at<uchar>(i, j + 1)) / 5;
				}
				else image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i - 2, j) + image2.at<uchar>(i + 2, j + 1) + image2.at<uchar>(i - 2, j + 1)) / 4;
			}

		}
	}
	if ((image2.at<uchar>(355, 54) - image2.at<uchar>(355, 56)) > 2) {
		cout << "hiiii" << endl;
	}
	cout << (image2.at<uchar>(355, 54) - image2.at<uchar>(355, 56)) << endl;
	cout << "hi" << endl; 
	image2 = image2.rowRange(2, image2.rows - 2).colRange(1, image2.cols - 1);
	
	Mat ImDif;
	absdiff(image, image2, ImDif);

	threshold(ImDif, ImDif, 2, 255, THRESH_BINARY);

	/*namedWindow("Dif", WINDOW_NORMAL);
	resizeWindow("Dif", image2.rows * 2, image2.cols * 2);
	imshow("Dif", ImDif);
	waitKey(0);*/

	namedWindow("Original image", WINDOW_NORMAL);
	resizeWindow("Original image", image2.rows * 2, image2.cols * 2);
	imshow("Original image", image);
	waitKey(0);

	namedWindow("New image", WINDOW_NORMAL);
	resizeWindow("New image", image2.rows * 2, image2.cols * 2);
	imshow("New image", image2);
	waitKey(0);


	return image2;
}
//Failed script
/*Mat RemoveLargeParticleNoise(Mat image) {
	//image.convertTo(image, CV_32FC1);
	Mat imagex=image.clone();
	Mat imagex2 = imagex.clone();
	cout << "here" << endl;

	int width = imagex.cols;   // Image width
	int height = imagex.rows;  // Image height

	for (int i = 2; i < height - 2; i++) {
		for (int j = width - 2; j >1; j--) {
			//cout << "here" << endl;
			imagex.at<uchar>(i, j) = (image.at<uchar>(i, j-1) - image.at<uchar>(i, j))/60;
		}
	}

	//threshold(imagex, imagex, 2, 255, THRESH_BINARY);

	imshow("1derv", imagex);
	waitKey(0);

	/*for (int i = 2; i < height - 2; i++) {
		for (int j = 2; j < width - 1; j++) {
			imagex2.at<float>(i, j) = (imagex.at<float>(i, j) - imagex.at<float>(i, j - 1))/100;
		}
	}

	imshow("2derv", imagex2);
	waitKey(0);


	/*int Hm = 400;
	unsigned char value = 0; // index value for the histogram
	int histogram[400]; // histogram array -remember to set to zero initially
	
	int k = Hm;
	while (k-- > 0) {
		histogram[k] = 0; // reset histogram entry
	}
	for (int i = 2; i < height-2; i++) {
		for (int j = 1; j < width-1; j++) {
			value = imagex.at<float>(i, j); //uchar is used for single-channel pixels
			//value = image.at<Vec3b>(i, j) [channel] ;
				//cout << value << endl;
			histogram[value] += 1;

		}
	}


	int bin_w = cvRound((double)width / Hm);
	Mat histImage(height, width, CV_8UC1, Scalar(255, 255, 255));
	//Mat histogram;

	//normalize(image, , 0, histImage.rows, NORM_MINMAX, -1, Mat());


	float maxVal = *max_element(histogram, histogram + Hm); // Find max value in array
	for (int i = 0; i < Hm; i++) {
		histogram[i] = (histogram[i] / maxVal) * height;  // Scale values to fit the image height
	}


	for (int i = 1; i < Hm; i++) {
		line(histImage, Point(bin_w * (i - 1), height - cvRound(histogram[i - 1])),
			Point(bin_w * i, height - cvRound(histogram[i])),
			Scalar(0, 0, 0), 2, 8, 0);  // Draw lines in black
	}
	namedWindow("Histogram", WINDOW_NORMAL);
	resizeWindow("Histogram", 600, 600);
	imshow("Histogram", histImage);

	waitKey(0);

	threshold(imagex, imagex, 100, 400, THRESH_BINARY);
	imshow("imagex2Bin", imagex);
	waitKey(0);

	return imagex;
}*/ 


//Remove the light column by comparing median of columns and local medians in 1 single image
//produced too much noise
Mat RemoveLightEdge(Mat image) { 
	Mat image2 = image.clone();

	int width = image2.cols;
	int height = image2.rows;
	int n = 0;

	vector<float> col;
	for (int i = 0; i < height; i++) {
		col.push_back(image2.at<uchar>(i, width-1));
		n++;
		//cout << i << "   " << col[i] << endl;
	}

	sort(col.begin(), col.end());
	float med=(col[n / 2]);

	vector<float> colerr;
	n = 0;
	int step = 2;
	vector<float> corr;
	int m = 0;

	for (int j =step;j<width-1; j = j + step * 2) {
		for (int i = step; i < height-1; i=i+step*2) {
			for (int k = i - step; k < i + step; k++) {
				for (int l = j - step; l < j + step; l++) {
					colerr.push_back(image2.at<uchar>(k, l) - med);
					corr.push_back(image2.at<uchar>(k, l) - med);
					n++;
					m++;
				}
			}
			sort(colerr.begin(), colerr.end());
			sort(corr.begin(), corr.end());
			float medErr = colerr[n / 2];
			float corErr = corr[m / 2];

			for (int k = i - step; k < i + step; k++) {
				for (int l = j - step; l < j + step; l++) {
					//image2.at<uchar>(k, l) = saturate_cast<uchar>(image2.at<uchar>(k, l) - medErr);
					if (image2.at<uchar>(k, l) - medErr >= med - 15) {
						image2.at<uchar>(k, l) = saturate_cast<uchar>(image2.at<uchar>(k, l) - medErr);
					}
					else {
						image2.at<uchar>(k, l) = saturate_cast<uchar>(image2.at<uchar>(k, l) - corErr);
					}
					m = 0;
					corr.clear();

				}
			}
		}
		colerr.clear();
		n = 0;
	}
	namedWindow("New image2", WINDOW_NORMAL);
	resizeWindow("New image2", image2.rows * 2, image2.cols * 2);
	imshow("New image2", image2);
	waitKey(0);

	return image2;
}

/*Mat MakeHistogram(Mat image) {
	int width = image.cols;   // Image width
	int height = image.rows;  // Image height
	int Hm = 256;
	unsigned char value = 0; // index value for the histogram
	int histogram[256]; // histogram array -remember to set to zero initially

	int k = Hm;
	while (k-- > 0) {
		histogram[k] = 0; // reset histogram entry
	}
	for (int i = 2; i < height - 2; i++) {
		for (int j = 1; j < width - 1; j++) {
			value = image.at<uchar>(i, j); //uchar is used for single-channel pixels
			//value = image.at<Vec3b>(i, j) [channel] ;
				//cout << value << endl;
			histogram[value] += 1;

		}
	}


	int bin_w = cvRound((double)width / Hm);
	Mat histImage(height, width, CV_8UC1, Scalar(255, 255, 255));
	//Mat histogram;

	//normalize(image, , 0, histImage.rows, NORM_MINMAX, -1, Mat());


	float maxVal = *max_element(histogram, histogram + Hm); // Find max value in array
	for (int i = 0; i < Hm; i++) {
		histogram[i] = (histogram[i] / maxVal) * height;  // Scale values to fit the image height
	}


	for (int i = 1; i < Hm; i++) {
		line(histImage, Point(bin_w * (i - 1), height - cvRound(histogram[i - 1])),
			Point(bin_w * i, height - cvRound(histogram[i])),
			Scalar(0, 0, 0), 2, 8, 0);  // Draw lines in black
	}
	namedWindow("Histogram", WINDOW_NORMAL);
	resizeWindow("Histogram", 600, 600);
	imshow("Histogram", histImage);

	waitKey(0);
	return histImage;
}*/

//comparison of profiles in graphs (not important)
Mat AnalyzeLight(Mat image) {
	Mat image2 = image.clone();
	
	int width = image2.cols;
	int height = image2.rows;
	int n = 0;

	vector<float> rowV;
	vector<float> colV;


	for (int i = 0; i < height; i++) {
		for (int j = 0; j < 20; j++) {
			rowV.push_back(image2.at<uchar>(i, j));
			n++;
		}
		sort(rowV.begin(), rowV.end());

		cout <<i<<"   "<< rowV[n / 2] << endl;
		colV.push_back(rowV[n / 2]);
		n = 0;
		rowV.clear();
	}

	Mat plotV = Mat::zeros(height, width, CV_8UC3);

	// Scaling factors
	float x_scale = width / static_cast<float>(colV.size() - 1);
	float y_max = *max_element(colV.begin(), colV.end());
	float y_scale = (height - 10) / y_max;

	// Draw the curve
	for (size_t i = 1; i < colV.size(); ++i) {
		Point pt1((i - 1) * x_scale, height - colV[i - 1] * y_scale);
		Point pt2(i * x_scale, height - colV[i] * y_scale);
		line(plotV, pt1, pt2, cv::Scalar(255, 0, 0), 2);
	}
	
	namedWindow("Vretical median profile", WINDOW_NORMAL);
	resizeWindow("Vretical median profile", image2.rows * 2, image2.cols * 2);
	imshow("Vretical median profile", plotV);
	waitKey(0);

	vector<float> rowH;
	vector<float> colH;

	for (int j = 0; j < width; j++) {
		for (int i = 0; i < height; i++) {
			colH.push_back(image2.at<uchar>(i, j));
			n++;
		}
		sort(colH.begin(), colH.end());

		cout << j << "   " << colH[n / 2] << endl;
		rowH.push_back(colH[n / 2]);
		n = 0;
		colH.clear();
	}

	Mat plotH = Mat::zeros(height, width, CV_8UC3);
	

	// Scaling factors
	float x_scaleH = width / static_cast<float>(rowH.size() - 1);
	float y_maxH = *max_element(rowH.begin(), rowH.end());
	float y_scaleH = (height - 10) / y_maxH;

	cout << "here" << endl;
	// Draw the curve
	for (size_t i = 1; i < rowH.size(); ++i) {
		Point pt1((i - 1) * x_scaleH, height - rowH[i - 1] * y_scaleH);
		Point pt2(i * x_scaleH, height - rowH[i] * y_scaleH);
		line(plotH, pt1, pt2, cv::Scalar(255, 0, 0), 2);
	}

	namedWindow("Horizontal median profile", WINDOW_NORMAL);
	resizeWindow("Horizontal median profile", image2.rows * 2, image2.cols * 2);
	imshow("Horizontal median profile", plotH);
	waitKey(0);

	return image;
}


//Noise removal with bluring (just for fun - not important)
Mat MedianFilter(Mat image) {
	//Attempt to remove noise with a median filter
	//It also works, but make the stars "fuzzier"
	Mat image2 = image.clone();
	//medianBlur(image, image2, 3);
	GaussianBlur(image, image2, Size(9, 9), 0);

	namedWindow("New image2", WINDOW_NORMAL);
	resizeWindow("New image2", image2.rows * 2, image2.cols * 2);
	imshow("New image2", image2);

	waitKey(0);
	return image2;
}


//Subtraction of column profile from scaled vs. unscaled (the darker) images
/*Mat RemoveLightColumnStar(Mat image, int type) {
	Mat MedSub;
	if (type == 6) {
		MedSub= imread("C:/Users/Bruger/Desktop/30330 project/MedianToSubtract1.bmp", 0);

		//MedSub = imread("C:/Users/Bruger/Desktop/30330 project/MedianToSubtractDark.bmp", 0);

	}
	else if (type == 1) {
		MedSub = imread("C:/Users/Bruger/Desktop/30330 project/MedianToSubtract1.bmp", 0);
		MedSub = MedSub / 3.2;
	}
	else {
		cout << "Wirte type either 1 or 6" << endl;
	}
	Mat image2 = image.clone();
	//GaussianBlur(MedSub, MedSub, Size(21, 21), 0);
	namedWindow("image with light", WINDOW_NORMAL);
	resizeWindow("image with light", image2.rows * 2, image2.cols * 2);
	imshow("image with light", MedSub);

	waitKey(0);
	image2 = image - MedSub;
	
	namedWindow("image without light", WINDOW_NORMAL);
	resizeWindow("image without light", image2.rows * 2, image2.cols * 2);
	imshow("image without light", image2);

	waitKey(0);
	
	return image2;

}*/


//Subtraction of column profile
Mat RemoveLightColumnProfile(Mat image) {
	Mat MedSub = imread("C:/Users/Bruger/Desktop/30330 project/MedianToSubtract1.bmp", 0); //profile made by comparing profiles of other images

	//Mat MedSub = imread("C:/Users/Bruger/Desktop/30330 project/MedianToSubtractDark.bmp", 0); //profile made by comparing dark images

	//Make a copy of the image to work with
	Mat image2 = image.clone();
	namedWindow("image with light", WINDOW_NORMAL);
	resizeWindow("image with light", image2.rows * 2, image2.cols * 2);
	imshow("image with light", MedSub);

	waitKey(0);

	//Simply subtract the median profile from the science image
	image2 = image - MedSub;

	namedWindow("image without light", WINDOW_NORMAL);
	resizeWindow("image without light", image2.rows * 2, image2.cols * 2);
	imshow("image without light", image2);

	waitKey(0);

	return image2;

}

//Subtraction of the lense reflection 
Mat RemoveLens(Mat image) {
	Mat lens= imread("C:/Users/Bruger/Desktop/30330 project/Lense2.bmp", 0);
	Mat image2;
	int width = image.cols;
	int height = image.rows;
	int n = 0;
	//The produced lense profile is much brighter than what is seen on the science immages. Therefore it should be scaled according to the image

	vector<float> IMG;
	vector<float> LENS;



	for (int i = 200; i < 300 ; i++) {
		for (int j = 300; j < 400 ; j++) {
			IMG.push_back(image.at<uchar>(i, j));
			LENS.push_back(lens.at<uchar>(i, j));
			n++;
		}
	}
	
	sort(IMG.begin(), IMG.end());
	sort(LENS.begin(), LENS.end());

	//cout << i << "   " << rowV[n / 2] << endl;
	float MedIm = IMG[n / 2];
	float MedL = LENS[n / 2];

	image2 = image - lens * (MedIm / MedL); //ratio is approx. 0.299
	cout << (MedIm / MedL) << endl;

	namedWindow("image with light", WINDOW_NORMAL);
	resizeWindow("image with light", image2.rows * 2, image2.cols * 2);
	imshow("image with light", lens*5);

	waitKey(0);

	namedWindow("image without light", WINDOW_NORMAL);
	resizeWindow("image without light", image2.rows * 2, image2.cols * 2);
	imshow("image without light", image2);

	waitKey(0);

	return image2;
}

