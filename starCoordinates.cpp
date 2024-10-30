//
// Spits out coordinates of blobs (stars) in console
// Input is binary image of the format produced by splitBinary: Binary image of every other line (with 0 in between)
// Change errosion and dilation to fit when using other format
// 

void blobCoordinates(Mat img) {
	// Input: Image with every other line removed
  // Same for


	Mat img_tmp;
	img.copyTo(img_tmp);		// Temporary image for analysis

	int width = img.cols;		// Parameter setup
	int height = img.rows;

	// Remove speckles
	img_tmp = errodeLR(img_tmp);
	img_tmp = errodeLR(img_tmp);
	img_tmp = errodeLR(img_tmp);

	// Remove left side (due to noise
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (j < 20) {

				img_tmp.at<uchar>(i, j) = 0;
			}

		}
	}

	// dilate back to make stars uniform 
	img_tmp = dilateCross(img_tmp);
	img_tmp = dilateCross(img_tmp);
	img_tmp = dilateCross(img_tmp);


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


	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {

			if (img_tmp.at<uchar>(i, j) == 255) {


				x_min = j;
				x_max = j;

				y_min = i;
				y_max = i;


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

					// Find next pixel, goes through 
					// Directions: 0 -> up, 1 -> right, 2 -> down, 3 -> left
					//cout << "step ";

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

						direction = (direction + 1) % 4;

					}


				}

				// Remove blob from image
				for (int k = y_min; k < y_max + 1; k++) {
					for (int l = x_min; l < x_max + 1; l++) {

						val = img_tmp.at<uchar>(k, l);
						if (val == 255) {
							// Add area
							area += 1;
							x_sum += l;
							y_sum += k;

							img_tmp.at<uchar>(k, l) = 0;
						}




					}

				}


				cout << "Blob at: " << x_sum / area << ", " << y_sum / area << ", area: " << area << endl;

				x_sum = 0;
				y_sum = 0;
				area = 0;
				border = 0;

			}

			

		}
	}

}
