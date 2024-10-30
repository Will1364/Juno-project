// Erosion and dilation of binary image
// input: binary image
// output: binary image of same size

Mat erodeLR(Mat img) {
  // Erriosin in a OXO format (only left and right)
  
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
  // Dilation in a cross (updown, leftright
  
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
