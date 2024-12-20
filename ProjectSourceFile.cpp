#include "HeaderFile.h"

Mat RemoveSmallParticleNoise(Mat image) {
	
	//Make a copy of the existing image to compare to the original
	Mat image2 = image.clone();
	copyMakeBorder(image2, image2, 2, 2, 1, 1, BORDER_CONSTANT, Scalar(0, 0, 0)); 

	//Set some threshold values
	int t1 = 55; //Sets the difference between a pixel and its imidiate neighbour 
	int t2 = 8; //Sets the difference between the 2 neighbours of a pixel

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
					else image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i-2, j + 1) + image2.at<uchar>(i - 2,j))/ 3;

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
			
			else if (j <= width - 2 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j - 1)) >= t1 && (image2.at<uchar>(i, j) - image2.at<uchar>(i, j +2)) >= t1 && (image2.at<uchar>(i, j - 1) - image2.at<uchar>(i, j + 1) <= t2)) { //her
				if (j == 1) {
					image2.at<uchar>(i, j) = (image2.at<uchar>(i - 2, j) + image2.at<uchar>(i + 2, j) + image2.at<uchar>(i + 2, j + 1) + image2.at<uchar>(i - 2, j + 1) + image2.at<uchar>(i, j + 1)) / 5;
				}
				else image2.at<uchar>(i, j) = (image2.at<uchar>(i, j - 1) + image2.at<uchar>(i - 2, j) + image2.at<uchar>(i + 2, j + 1) + image2.at<uchar>(i - 2, j + 1)) / 4;
			}

		}
	} 
	image2 = image2.rowRange(2, image2.rows - 2).colRange(1, image2.cols - 1);
	
	//If desired show image depicting residuals
	Mat ImDif;
	absdiff(image, image2, ImDif);

	threshold(ImDif, ImDif, 2, 255, THRESH_BINARY);

	/*namedWindow("Dif", WINDOW_NORMAL);
	resizeWindow("Dif", image2.rows * 2, image2.cols * 2);
	imshow("Dif", ImDif);
	waitKey(0);

	namedWindow("Original image", WINDOW_NORMAL);
	resizeWindow("Original image", image2.rows * 2, image2.cols * 2);
	imshow("Original image", image);
	waitKey(0);

	namedWindow("New image", WINDOW_NORMAL);
	resizeWindow("New image", image2.rows * 2, image2.cols * 2);
	imshow("New image", image2);
	waitKey(0);*/


	return image2;
}


//Remove the light column by comparing median of columns and local medians in 1 single image
Mat RemoveLightEdge(Mat image) { 
	//Make a copy of the image to work with
	Mat image2 = image.clone();

	int width = image2.cols;
	int height = image2.rows;
	int n = 0;

	//Create array of pixel values in the last column
	vector<float> col;
	for (int i = 0; i < height; i++) {
		col.push_back(image2.at<uchar>(i, width-1));
		n++;
	}
	//Find the median of the last column
	//Sort array values
	sort(col.begin(), col.end()); 
	//Find middle value
	float med=(col[n / 2]); 

	//create empty array for median errors for column groups
	vector<float> colerr;
	n = 0;
	//The half-size of the compare pixelcluster
	int step =3;
	//create empty array for median errors for local groups
	vector<float> locerr;
	int m = 0;

	//Go through every pixel in the image coumnwise
	for (int j =step;j<width-1; j = j + step * 2) {
		for (int i = step; i < height-1; i=i+step*2) {
			//For every pixel loop through the neighbour cluster
			for (int k = i - step; k < i + step; k++) {
				for (int l = j - step; l < j + step; l++) {
					colerr.push_back(image2.at<uchar>(k, l) - med);
					locerr.push_back(image2.at<uchar>(k, l) - med);
					n++;
					m++;
				}
			}
			//Find median for the column groups and the local groups
			sort(colerr.begin(), colerr.end());
			sort(locerr.begin(), locerr.end());
			float medErr = colerr[n / 2];
			float locErr = locerr[m / 2];

			//For every pixel loop through the neighbour cluster again
			for (int k = i - step; k < i + step; k++) {
				for (int l = j - step; l < j + step; l++) {
					//To avoid overcorrection, the limit for using column median is set to 15 below the median of the noise-free column
					if (image2.at<uchar>(k, l) - medErr >= med - 15) {
						image2.at<uchar>(k, l) = saturate_cast<uchar>(image2.at<uchar>(k, l) - medErr);
					}
					//If the image is overcorrected, the local median is used instead
					else {
						image2.at<uchar>(k, l) = saturate_cast<uchar>(image2.at<uchar>(k, l) - locErr);
					}
					m = 0;
					locerr.clear();

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
	//GaussianBlur(image, image2, Size(9, 9), 0);

	//Create emtpy array
	vector<float> medKern;
	int width = image2.cols;
	int height = image2.rows;
	int n = 0;
	int halfKernSize = 1;
	//Loop through the image
	for (int i = halfKernSize; i < height - halfKernSize; i++) {
		cout << i << endl;
		for (int j = halfKernSize; j < width - halfKernSize; j++) {
			//For each image find intensities of the neighbours
			for (int k = -halfKernSize; k <= halfKernSize; k++) {
				for (int l = -halfKernSize; l <= halfKernSize; l++) {
					medKern.push_back(image.at<uchar>(i + k, j + l));
					n++;
				}
			}
			//find the median
			sort(medKern.begin(), medKern.end());
			float median = medKern[n / 2];
			image2.at<uchar>(i, j) = saturate_cast<uchar>(median);
			n = 0;
			medKern.clear();
		}
	}

	/*namedWindow("New image2", WINDOW_NORMAL);
	resizeWindow("New image2", image2.rows * 2, image2.cols * 2);
	imshow("New image2", image2);
	*/
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

	//The images will be scaled according to the medians of the images

	for (int i = 200; i < 300 ; i++) {
		for (int j = 300; j < 400 ; j++) {
			IMG.push_back(image.at<uchar>(i, j));
			LENS.push_back(lens.at<uchar>(i, j));
			n++;
		}
	}
	
	sort(IMG.begin(), IMG.end());
	sort(LENS.begin(), LENS.end());

	
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


//Removes Jupiter for reflection modelling
Mat removeJup(Mat image) {
	//blur the image from the beginning to remove unwanted noise and artifacts
	GaussianBlur(image, image, Size(5, 5), 0);
	int width = image.cols;
	int height = image.rows;
	Mat image2 = image.clone();
	Mat image3 = image.clone();
	imshow("im", image);
	
	int sum = 0;
	//morphological watershed: set all values below a certain point to 0
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (image3.at<uchar>(i, j) < 43) { //27 fits for RightH, 32 fits for LeftH, 62 might fit LeftL
				image3.at<uchar>(i, j) = 0;
			}
		}
	}
	
	//Creation of different Gaussian filters
	float sig3 = 0.5;
	Mat GF3p5 = Mat::zeros(3, 3, CV_32F);
	for (int x = -1; x <= 1;x++) {
		for (int y = -1; y <= 1;y++) {
			GF3p5.at<float>(x + 1, y + 1) = 1 / (2 * pi * sig3) * exp(-(pow(x,2) + pow(y,2)) / (2 * sig3));

			cout << GF3p5.at<float>(x + 1, y + 1) << endl;
		}
	}
	float sig5 = 1;
	Mat GF5 = Mat::zeros(5, 5, CV_32F);
	for (int x = -2; x <= 2;x++) {
		for (int y = -2; y <= 2;y++) {
			GF5.at<float>(x + 2, y + 2) = 1 / (2 * pi * sig5) * exp(-(pow(x, 2) + pow(y, 2)) / (2 * sig5));
		}
	}

	float sig11 = 2.5;
	Mat GF11 = Mat::zeros(11, 11, CV_32F);
	for (int x = -5; x <= 5;x++) {
		for (int y = -5; y <= 5;y++) {
			GF11.at<float>(x + 5, y + 5) = 1 / (2 * pi * sig11) * exp(-(pow(x, 2) + pow(y, 2)) / (2 * sig11));
		}
	}
	
	//apply Gaussian filter to watershedded image
	for (int i = 1; i < height - 1; i++) {
		for (int j = 1; j < width - 1; j++) {
			float sum = 0.0; 
			for (int k = -1; k <= 1; k++) {
				for (int l = -1; l <= 1; l++) {
					int filterRow = k + 1; 
					int filterCol = l + 1; 
					sum += image3.at<uchar>(i + k, j + l) * GF3p5.at<float>(filterRow, filterCol);
				}
			}
			// Clamp the result to [0, 255] and store in the output image
			image3.at<uchar>(i, j) = cv::saturate_cast<uchar>(sum);
		}
	}

	namedWindow("LP2", WINDOW_NORMAL);
	resizeWindow("LP2", image2.rows * 2, image2.cols * 2);
	imshow("LP2", image3);
	waitKey();
	
	//Set upper limit for when to stop removing Jupiter and atmosphere
	int lim = 141;
	float blend_factor[141]; //41 for LeftL, 121 for RightH and mid, 101 for LeftL, 71 for mid
	for (int i = 0; i < lim; i++) {
		blend_factor[i] = static_cast<double>(i) / (lim - 1);
	}
	//For each pixel where Jupiter is detected replace pixel value with background estimate
	//For the extended columns a blend factor is introduces to blend into the original image
	for (int i = width-1; i > 0; i--){
		//cout << i << endl;
		for (int j = height-3; j > 0; j--) {
			if (i < 752) {
				//cout <<"J"<< j << endl;
				image.at<uchar>(j, i) = rand() % 10 + 23;
				if (image3.at<uchar>(j, i) == 0) {
					if (j > lim+1)
						for (int k = j - 1; k > j - lim; k--) {
							//cout << k << endl;
							image.at<uchar>(k, i) = (rand() % 10 + 23) * (1 - blend_factor[j - k - 1]) + image.at<uchar>(k, i) * (blend_factor[j - k - 1]); // 47 fits for LeftL, 15 fits for H, 25 for mid
						}
					break;
				}
			}
		}
	}

	//Apply another filter to the finished image
	for (int i = 5; i < height - 5; i++) {
		for (int j = 5; j < width - 5; j++) {
			float sum = 0.0; // Ensure sum is a float to handle the convolution properly
			for (int k = -5; k <= 5; k++) {
				for (int l = -5; l <= 5; l++) {
					// Calculate filter indices correctly
					int filterRow = k + 5; // Offset by 1 to map [-1, 0, 1] to [0, 1, 2]
					int filterCol = l + 5; // Same offset logic
					sum += image.at<uchar>(i + k, j + l) * GF11.at<float>(filterRow, filterCol);
				}
			}
			// Clamp the result to [0, 255] and store in the output image
			image.at<uchar>(i, j) = cv::saturate_cast<uchar>(sum);
		}
	}
	


	namedWindow("LP3", WINDOW_NORMAL);
	resizeWindow("LP3", image2.rows * 2, image2.cols * 2);
	imshow("LP3", image);
	waitKey();

	return image2;
}
