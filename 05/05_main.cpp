#include <opencv2\opencv.hpp>
#include <iostream>
#include <cmath>
#include "05.h"

using namespace cv;
using namespace std;


int main() {
	string imageName("C:\\Users\\samsung\\Desktop\\BEP\\BEP5\\test_lable1.bmp");
	Mat img = imread(imageName.c_str(), IMREAD_COLOR), gray;
	cvtColor(img, gray, COLOR_RGB2GRAY);

	imshow("Original", gray); waitKey();
	imshow("Flood fill", Flood_fill(gray)); waitKey();
	//imshow("Two pass", Two_pass(gray)); waitKey();

	return 0;
}