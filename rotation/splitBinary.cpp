// Similar to splitImage except it keeps the lines in between valued at zero to preserve coordinates 

pair<Mat, Mat> splitBinary(Mat img) {

	int width = img.cols;		// Parameter setup
	int height = img.rows;

	Mat img_bn1, img_bn2;
	img_bn1 = img > 230;
	img_bn2 = img > 230;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (i % 2 == 0) {

				img_bn2.at<uchar>(i, j) = 0;

			}

			else {

				img_bn1.at<uchar>(i, j) = 0;

			}


		}
	}

	return make_pair(img_bn1, img_bn2);

}
