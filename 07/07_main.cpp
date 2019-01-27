#include <opencv2\opencv.hpp>
#include <iostream>
#include <cmath>
#include "03.h"
#include "04.h"
#include "07.h"

using namespace cv;
using namespace std;

const float PI = 3.14159265358979323846f;

int main() {
	string imageName("C:\\Users\\samsung\\Desktop\\BEP\\BEP7\\image.jpg");
	Mat img = imread(imageName.c_str(), IMREAD_COLOR), gray;
	cvtColor(img, gray, COLOR_RGB2GRAY);
	imshow("Original", img); waitKey();
	//imshow("Color", For_color(img)); waitKey();
	//imshow("Intensity", For_intensity(img)); waitKey();
	//imshow("orienation", For_orientation(img)); waitKey();
	imshow("Saliency", Saliency(img)); waitKey();
	return 0;
}