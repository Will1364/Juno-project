// this code is made with inspiration from https://insightfultscript.com/collections/programming/cpp/opencv-fourier-transform-cpp/


// special dependencies
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>




// To be inserted into the source file:

Mat fourier(const Mat& OGImage) {
	// -------settings-------------------------------------------------
	
	int mode = 1; 
	// 0 for just setting oscillation frequency to 0 
	// 1 for gaussian 

	string padding = "mirror";
	// options for padding are: "zeros", "clamp", "mirror", "wrap", "old version"
	
	int analysis = 0; 
	// if you want to see the magnitude profile and images set to 1
	
	int manual = true; 
	// set to false if you want the program to automatically detect the highest non-dc oscillation and remove it
	
	vector<int> fList = {307, 230, 153}; 
	// list of known frequencies to be removed in the manual setting

	int vvInterval = 0; 
	// half the vInterval width to be removed in y-direction is defined (Set to 0 for just changing values at the x-axis)

	
	const int halfFilterWidth = 30;  // controls width of gaussian filter
	const int cutoff =  10;		// controls agressiveness of gaussian filter
	
	// ------- program -----------------------------------------------------
	
	int rows = OGImage.rows;
	int cols = OGImage.cols;
	
	//cout << "Rows = " << rows << "  cols = " << cols << "\n";

	// First we need to resize the image so dimensions are optimal for DFT transform.
	// This can be done by padding all the edges around the image. The padding type is chosen in the settings
	

	int desiredRows = getOptimalDFTSize(OGImage.rows);
	int desiredCols = getOptimalDFTSize(OGImage.cols);
	Mat paddedImage(desiredRows, desiredCols, OGImage.type());
	
	//cout << "Desired rows = " << desiredRows << "  cols = " << desiredCols << "\n";

	int vBorder = (desiredRows - rows) / 2;
	int hBorder = (desiredCols - cols) / 2;
	

	// Here the padding is performed
	if (padding == "clamp") {
		copyMakeBorder(OGImage, paddedImage, vBorder, vBorder, hBorder, hBorder, BORDER_REPLICATE); 
	}
	else if (padding == "zeros") {
		Mat zeroRow = Mat::zeros(vBorder, cols, OGImage.type());
		vconcat(OGImage, zeroRow, paddedImage);
		vconcat(zeroRow, paddedImage, paddedImage);
		Mat zeroCol = Mat::zeros(paddedImage.rows, hBorder, OGImage.type()); 
		hconcat(paddedImage, zeroCol, paddedImage); 
		hconcat(zeroCol, paddedImage, paddedImage); 
	
	}
	else if (padding == "mirror") {
		copyMakeBorder(OGImage, paddedImage, vBorder, vBorder, hBorder, hBorder, BORDER_REFLECT101);
	}
	else if (padding == "wrap") {
		copyMakeBorder(OGImage, paddedImage, vBorder, vBorder, hBorder, hBorder, BORDER_WRAP);
	}
	else if (padding == "old version") {
		Mat zeroRow = Mat::zeros(desiredRows - rows, cols, OGImage.type()); 
		vconcat(OGImage, zeroRow, paddedImage); 
		Mat zeroCol = Mat::zeros(paddedImage.rows, desiredCols - cols, OGImage.type());   
		hconcat(paddedImage, zeroCol, paddedImage); 
	}
	else {
		cout << "No valid padding is chosen" << "\n";
		waitKey(0);
	}

	

	//imshow("padded image", paddedImage);
	//waitKey();

	//cout << "New rows = " << paddedImage.rows << "  new Cols = " << paddedImage.cols << "\n";
	
	
	// Now we will create the complex matrix which will contain the DFT output
	
	Mat real = Mat_<float>(paddedImage);
	Mat imaginary = Mat::zeros(paddedImage.size(),CV_32F); 
	// the output of the dft is a floating point value, so the two matrices are defined as float

	Mat channels[] = { real, imaginary };
	Mat complexImage; 
	merge(channels, 2, complexImage);

	// Doing the DFT !!! :-)
	dft(complexImage, complexImage, DFT_COMPLEX_OUTPUT);

	// finding the magnitude as sqrt(R^2 + I^2) for each element in the complex output
	split(complexImage, channels);
	Mat mag;
	magnitude(channels[0], channels[1], mag);


	double minVal, maxVal;
	minMaxLoc(mag, &minVal, &maxVal);
	//cout << "min mag is: " << minVal << "  max mag is: " << maxVal << "\n";

	if (minVal == 0) {
		mag += 1; // avoid accidentally deviding by zero when taking the natural logarithm
	}

	log(mag, mag);
	minMaxLoc(mag, &minVal, &maxVal);
	//cout << "min mag is: " << minVal << "  max mag is: " << maxVal << "\n";

	// values in floating images needs to be between 0-1 so we need to normalize
	normalize(mag, mag, 0, 1, cv::NORM_MINMAX); 
	minMaxLoc(mag, &minVal, &maxVal);
	//cout << "min mag is: " << minVal << "  max mag is: " << maxVal << "\n";

	// to make sure the coordinates (0,0) is at the centre of the image, the first quadrant is swapped with the third, an´d the second with the fourth
	mag = mag(Rect(0, 0, mag.cols & -2, mag.rows & -2)); // reduces the magnitude image to an even number of columns and rows

	int xCenter = mag.cols / 2;
	int yCenter = mag.rows / 2;
	mag = swapQuadrants(mag); 

	if (analysis == 1){
		imshow("fourier image", mag); 
		waitKey();
	}

	// now a graph of the x-axis is created so the oscilation frequency is more easily recognized
	Mat xAxis = mag.row(yCenter).clone();  //deep copy
	minMaxLoc(xAxis, &minVal, &maxVal);
	xAxis = (xAxis / maxVal)*500; // Normalize and rescale to fit window

	if (analysis == 1){
		int length = xAxis.cols;
		Mat graph(500, length, CV_8UC3, Scalar(255, 255, 255));
		
		float value;
		int roundedValue;
		int thickness = 1;     //line thickness is defined
		for (int i = 0; i < length-1; i++) {               // a line is drawn for each value/element in the graph
		
			value = xAxis.at<float>(0, i);
			roundedValue = static_cast<int>(round(value));
			Point p1(i , 500 - roundedValue);             // the two points are defined
		
			value = xAxis.at<float>(0, i + 1);
			roundedValue = static_cast<int>(round(value));
			Point p2((i + 1), 500 - roundedValue);
		
			line(graph, p1, p2, Scalar(0, 0, 0),               //here the actual line is drawn
				thickness, LINE_8);
		}
	
		namedWindow("graph", WINDOW_AUTOSIZE);
		imshow("graph",graph); 
		waitKey(); 
	}
	
	
	// the central peak is taken out of the profile to remove the DC contribution
	for (int i = 0; i < 40; i++) {				
		xAxis.at<float>(0, xCenter + i) = 0;
		xAxis.at<float>(0, xCenter - i) = 0;
	}
	Point minLoc;
	Point maxLoc;
	minMaxLoc(xAxis, &minVal, &maxVal, &minLoc, &maxLoc); // the largest oscillation is located
	
	//cout << maxLoc.x << "\n";
	int diff = xCenter - maxLoc.x;
	mag(Range(yCenter - vInterval, yCenter + vInterval), Range(xCenter + diff - 2, xCenter + diff + 2)) = 0; // target oscilation is set to 0 
	mag(Range(yCenter - vInterval, yCenter + vInterval), Range(xCenter - diff - 2, xCenter - diff + 2)) = 0; 

	if (analysis == 1){
		imshow("filterd fourier image", mag);
		waitKey();
	}
	
	// now we can filter out the target oscillation according to the mode selected
	complexImage = swapQuadrants(complexImage); 
	if (mode == 0) {
		if (manual == false) { // Applying the filter to frequencies detected automatically
			int diff = xCenter - maxLoc.x;
			//complexImage(Range(yCenter, yCenter + yCenter), Range(xCenter + ndiff - 5, xCenter + ndiff + 5)) = 0; // target oscilation is set to 0 
			//complexImage(Range(yCenter, yCenter + yCenter), Range(xCenter - ndiff - 5, xCenter - ndiff + 5)) = 0;  
			complexImage(Range(yCenter - vInterval, yCenter + vInterval), Range(xCenter + diff - 5, xCenter + diff + 5)) = 0; // target oscilation is set to 0  
			complexImage(Range(yCenter - vInterval, yCenter + vInterval), Range(xCenter - diff - 5, xCenter - diff + 5)) = 0; 
			mag(Range(yCenter - vInterval, yCenter + vInterval), Range(xCenter + diff - 5, xCenter + diff + 5)) = 0; // target oscilation is set to 0   
			mag(Range(yCenter - vInterval, yCenter + vInterval), Range(xCenter - diff - 5, xCenter - diff + 5)) = 0; 
			complexImage = swapQuadrants(complexImage); 
		}
		else if (manual == true) { // Applying the filter to frequencies specified manually
			for (int i = 0; i < fList.size(); i++) {
				int diff = xCenter - fList[i];
				complexImage(Range(yCenter - vInterval, yCenter + vInterval), Range(xCenter + diff - halfFilterWidth, xCenter + diff + halfFilterWidth)) = 0; // target oscilation is set to 0   
				complexImage(Range(yCenter - vInterval, yCenter + vInterval), Range(xCenter - diff - halfFilterWidth, xCenter - diff + halfFilterWidth)) = 0; 
				mag(Range(yCenter - vInterval, yCenter + vInterval), Range(xCenter + diff - halfFilterWidth, xCenter + diff + halfFilterWidth)) = 0; // target oscilation is set to 0    
				mag(Range(yCenter - vInterval, yCenter + vInterval), Range(xCenter - diff - halfFilterWidth, xCenter - diff + halfFilterWidth)) = 0;  
			}
			complexImage = swapQuadrants(complexImage); 
		}
		imshow("filterd fourier image", mag);
		waitKey();
	}
	else if (mode == 1) {
		float halfGaussFilter[halfFilterWidth]; 
		for (int i = 0;i < halfFilterWidth; i++) { // creating the gaussian filter
			halfGaussFilter[i] = 1 - exp(-pow(i + 1, 2) / (2 * pow(cutoff, 2)));
			// cout << halfGaussFilter[i] << "\n";
		}
		split(complexImage, channels);
		if (manual == false) {  // Applying the filter to frequencies detected automatically
			int diff = xCenter - maxLoc.x;
			for (int i = 0; i <= vInterval; i++) {

				
				// changing the target frequencies in magnitude images
				mag.at<float>(yCenter + i, xCenter + diff) = 0;
				mag.at<float>(yCenter - i, xCenter + diff) = 0;
				mag.at<float>(yCenter + i, xCenter - diff) = 0;
				mag.at<float>(yCenter - i, xCenter - diff) = 0;

				// changing the target frequencies in complex images
				channels[0].at<float>(yCenter + i, xCenter + diff) = 0;
				channels[0].at<float>(yCenter - i, xCenter + diff) = 0;
				channels[0].at<float>(yCenter + i, xCenter - diff) = 0;
				channels[0].at<float>(yCenter - i, xCenter - diff) = 0;
				channels[1].at<float>(yCenter + i, xCenter + diff) = 0;
				channels[1].at<float>(yCenter - i, xCenter + diff) = 0;
				channels[1].at<float>(yCenter + i, xCenter - diff) = 0;
				channels[1].at<float>(yCenter - i, xCenter - diff) = 0;
		
				for (int j = 0; j < halfFilterWidth; j++) {

					//Applying Gaussian filter to the complex image
					// first quadrant
					channels[0].at<float>(yCenter - i, xCenter + diff + (j + 1)) *= halfGaussFilter[j];
					channels[0].at<float>(yCenter - i, xCenter + diff - (j + 1)) *= halfGaussFilter[j];
		
					//second quadrant
					channels[0].at<float>(yCenter - i, xCenter - diff + (j + 1)) *= halfGaussFilter[j];
					channels[0].at<float>(yCenter - i, xCenter - diff - (j + 1)) *= halfGaussFilter[j];
		
					// third quadrant
					channels[0].at<float>(yCenter + i, xCenter - diff + (j + 1)) *= halfGaussFilter[j];
					channels[0].at<float>(yCenter + i, xCenter - diff - (j + 1)) *= halfGaussFilter[j];
		
					// fourth  quadrant
					channels[0].at<float>(yCenter + i, xCenter + diff + (j + 1)) *= halfGaussFilter[j];
					channels[0].at<float>(yCenter + i, xCenter + diff - (j + 1)) *= halfGaussFilter[j];
		
					// first quadrant
					channels[1].at<float>(yCenter - i, xCenter + diff + (j + 1)) *= halfGaussFilter[j];
					channels[1].at<float>(yCenter - i, xCenter + diff - (j + 1)) *= halfGaussFilter[j];
		
					//second quadrant
					channels[1].at<float>(yCenter - i, xCenter - diff + (j + 1)) *= halfGaussFilter[j];
					channels[1].at<float>(yCenter - i, xCenter - diff - (j + 1)) *= halfGaussFilter[j];
		
					// third quadrant
					channels[1].at<float>(yCenter + i, xCenter - diff + (j + 1)) *= halfGaussFilter[j];
					channels[1].at<float>(yCenter + i, xCenter - diff - (j + 1)) *= halfGaussFilter[j];
		
					// fourth  quadrant
					channels[1].at<float>(yCenter + i, xCenter + diff + (j + 1)) *= halfGaussFilter[j];
					channels[1].at<float>(yCenter + i, xCenter + diff - (j + 1)) *= halfGaussFilter[j];

					// drawing changes in the magnitude plot
					// first quadrant
					mag.at<float>(yCenter - i, xCenter + diff + (j + 1)) = mag.at<float>(yCenter - i, xCenter + diff + (j + 1)) * halfGaussFilter[j];
					mag.at<float>(yCenter - i, xCenter + diff - (j + 1)) = mag.at<float>(yCenter - i, xCenter + diff - (j + 1)) * halfGaussFilter[j];
		
					//second quadrant
					mag.at<float>(yCenter - i, xCenter - diff + (j + 1)) = mag.at<float>(yCenter - i, xCenter - diff + (j + 1)) * halfGaussFilter[j];
					mag.at<float>(yCenter - i, xCenter - diff - (j + 1)) = mag.at<float>(yCenter - i, xCenter - diff - (j + 1)) * halfGaussFilter[j];
		
					// third quadrant
					mag.at<float>(yCenter + i, xCenter - diff + (j + 1)) = mag.at<float>(yCenter + i, xCenter - diff + (j + 1)) * halfGaussFilter[j];
					mag.at<float>(yCenter + i, xCenter - diff - (j + 1)) = mag.at<float>(yCenter + i, xCenter - diff - (j + 1)) * halfGaussFilter[j];
		
					// fourth  quadrant
					mag.at<float>(yCenter + i, xCenter + diff + (j + 1)) = mag.at<float>(yCenter + i, xCenter + diff + (j + 1)) * halfGaussFilter[j];
					mag.at<float>(yCenter + i, xCenter + diff - (j + 1)) = mag.at<float>(yCenter + i, xCenter + diff - (j + 1)) * halfGaussFilter[j];
				}
			}
		}
		else if (manual == true) { // Applying the filter to frequencies specified manually
			for (int k = 0; k < fList.size(); k++) {
				int diff = xCenter - fList[k];
		
				for (int i = 0; i <= vInterval; i++) {
					mag.at<float>(yCenter + i, xCenter + diff) = 0;
					mag.at<float>(yCenter - i, xCenter + diff) = 0;
					mag.at<float>(yCenter + i, xCenter - diff) = 0;
					mag.at<float>(yCenter - i, xCenter - diff) = 0;
		
					channels[0].at<float>(yCenter + i, xCenter + diff) = 0;
					channels[0].at<float>(yCenter - i, xCenter + diff) = 0;
					channels[0].at<float>(yCenter + i, xCenter - diff) = 0;
					channels[0].at<float>(yCenter - i, xCenter - diff) = 0;
					channels[1].at<float>(yCenter + i, xCenter + diff) = 0;
					channels[1].at<float>(yCenter - i, xCenter + diff) = 0;
					channels[1].at<float>(yCenter + i, xCenter - diff) = 0;
					channels[1].at<float>(yCenter - i, xCenter - diff) = 0;
		
					for (int j = 0; j < halfFilterWidth; j++) {
						
						// first quadrant
						channels[0].at<float>(yCenter - i, xCenter + diff + (j + 1)) *= halfGaussFilter[j];
						channels[0].at<float>(yCenter - i, xCenter + diff - (j + 1)) *= halfGaussFilter[j];
		
						//second quadrant
						channels[0].at<float>(yCenter - i, xCenter - diff + (j + 1)) *= halfGaussFilter[j];
						channels[0].at<float>(yCenter - i, xCenter - diff - (j + 1)) *= halfGaussFilter[j];
		
						// third quadrant
						channels[0].at<float>(yCenter + i, xCenter - diff + (j + 1)) *= halfGaussFilter[j];
						channels[0].at<float>(yCenter + i, xCenter - diff - (j + 1)) *= halfGaussFilter[j];
		
						// fourth  quadrant
						channels[0].at<float>(yCenter + i, xCenter + diff + (j + 1)) *= halfGaussFilter[j];
						channels[0].at<float>(yCenter + i, xCenter + diff - (j + 1)) *= halfGaussFilter[j];
		
						// first quadrant
						channels[1].at<float>(yCenter - i, xCenter + diff + (j + 1)) *= halfGaussFilter[j];
						channels[1].at<float>(yCenter - i, xCenter + diff - (j + 1)) *= halfGaussFilter[j];
		
						//second quadrant
						channels[1].at<float>(yCenter - i, xCenter - diff + (j + 1)) *= halfGaussFilter[j];
						channels[1].at<float>(yCenter - i, xCenter - diff - (j + 1)) *= halfGaussFilter[j];
		
						// third quadrant
						channels[1].at<float>(yCenter + i, xCenter - diff + (j + 1)) *= halfGaussFilter[j];
						channels[1].at<float>(yCenter + i, xCenter - diff - (j + 1)) *= halfGaussFilter[j];
		
						// fourth  quadrant
						channels[1].at<float>(yCenter + i, xCenter + diff + (j + 1)) *= halfGaussFilter[j];
						channels[1].at<float>(yCenter + i, xCenter + diff - (j + 1)) *= halfGaussFilter[j];
		
		
		
						// first quadrant
						mag.at<float>(yCenter - i, xCenter + diff + (j + 1)) = mag.at<float>(yCenter - i, xCenter + diff + (j + 1)) * halfGaussFilter[j];
						mag.at<float>(yCenter - i, xCenter + diff - (j + 1)) = mag.at<float>(yCenter - i, xCenter + diff - (j + 1)) * halfGaussFilter[j];
		
						//second quadrant
						mag.at<float>(yCenter - i, xCenter - diff + (j + 1)) = mag.at<float>(yCenter - i, xCenter - diff + (j + 1)) * halfGaussFilter[j];
						mag.at<float>(yCenter - i, xCenter - diff - (j + 1)) = mag.at<float>(yCenter - i, xCenter - diff - (j + 1)) * halfGaussFilter[j];
		
						// third quadrant
						mag.at<float>(yCenter + i, xCenter - diff + (j + 1)) = mag.at<float>(yCenter + i, xCenter - diff + (j + 1)) * halfGaussFilter[j];
						mag.at<float>(yCenter + i, xCenter - diff - (j + 1)) = mag.at<float>(yCenter + i, xCenter - diff - (j + 1)) * halfGaussFilter[j];
		
						// fourth  quadrant
						mag.at<float>(yCenter + i, xCenter + diff + (j + 1)) = mag.at<float>(yCenter + i, xCenter + diff + (j + 1)) * halfGaussFilter[j];
						mag.at<float>(yCenter + i, xCenter + diff - (j + 1)) = mag.at<float>(yCenter + i, xCenter + diff - (j + 1)) * halfGaussFilter[j];
					}
				}
			}
		}
		merge(channels, 2, fcomplexImage);  
		complexImage = swapQuadrants(fcomplexImage); 
		
		if (analysis == 1){
			imshow("filterd fourier image", mag); 
			waitKey();
		}	
	}

	// at last we can convert back to the reconstructed image:

	
	idft(complexImage, complexImage, DFT_SCALE | DFT_REAL_OUTPUT); 
	
	minMaxLoc(complexImage, &minVal, &maxVal, &minLoc, &maxLoc);
	cout << "min val is: " << minVal << "  max val is: " << maxVal << "\n";
	
	Mat recImage;

	int myCount0 = 0;
	int myCount255 = 0;

	// Counting the number of pixels being clipped
	for (int i = 0; i < (cols + 2 * hBorder); i++) {
		for (int j = 0; j < (rows + 2 * vBorder); j++) {
			if (complexImage.at<float>(j, i) < 0) {
				complexImage.at<float>(j, i) = 0;
				myCount0++;
			}
			else if (complexImage.at<float>(j, i) > 255) {
				complexImage.at<float>(j, i) = 255;
				myCount255++;
			}
		}
	}

	
	cout << "myCount0: " << myCount0 << "  myCount255: " << myCount255 << "\n";
	complexImage.copyTo(recImage);

	// cutting off the patting
	if (padding == "old version") {
		recImage = recImage.colRange(0, cols);
		recImage = recImage.rowRange(0, rows);
	}
	else {
		recImage = recImage.colRange(hBorder, cols + hBorder);
		recImage = recImage.rowRange(vBorder, rows + vBorder);
	}

	// calculating the residuals
	Mat Residuals;
	Mat OGImage2;
	OGImage.copyTo(OGImage2);
	OGImage2.convertTo(OGImage2, CV_32FC1);
	
	Residuals = OGImage2 - recImage;
	//Residuals = Residuals(Range(15, 400), Range(15, 400));
	
	minMaxLoc(Residuals, &minVal, &maxVal, &minLoc, &maxLoc);
	cout << "min val is: " << minVal << "  max val is: " << maxVal << "\n";
	Residuals = Residuals * (255 / maxVal);

	// saving the images
	imwrite("Residuals_2f.png", Residuals);
	
	//imshow("Reconstructed Image (IDFT)", recImage);
	imwrite("recImage_2f.png", recImage);
	
	
	imshow("Residuals", abs(Residuals / 255.0));
	return recImage / 255.0;
}
