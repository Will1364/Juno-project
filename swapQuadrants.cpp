Mat swapQuadrants(Mat image) {
  // this function is necessary for the DFT function to work
	int xCenter = image.cols / 2;
	int yCenter = image.rows / 2;

	Mat firstQuard(image, Rect(xCenter, 0, xCenter, yCenter));
	Mat secondQuard(image, Rect(0, 0, xCenter, yCenter));
	Mat thirdQuard(image, Rect(0, yCenter, xCenter, yCenter));
	Mat fourthQuard(image, Rect(xCenter, yCenter, xCenter, yCenter));

	Mat Temp;

	//swap first and third quadrant
	firstQuard.copyTo(Temp);
	thirdQuard.copyTo(firstQuard);
	Temp.copyTo(thirdQuard);

	//swap second and third quadrant
	secondQuard.copyTo(Temp);
	fourthQuard.copyTo(secondQuard);
	Temp.copyTo(fourthQuard);
	return image;
}
