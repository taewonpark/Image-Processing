#include <opencv2\opencv.hpp>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;
using namespace cv;

void RGB2HSV(double b, double g, double r);
void RGB2CMYK(double b, double g, double r);
void RGB2RGBY(double b, double g, double r);

int main() {
	string imageName("C:\\Users\\samsung\\Desktop\\BEP\\lena.jpg");
	Mat img = imread(imageName.c_str(), IMREAD_COLOR), hsv;

	namedWindow("image", WINDOW_AUTOSIZE);
	
	if (img.empty()){
		std::cout << "Could not open or find the image" << std::endl;
		return -1;
	}
	
	imshow("lena", img);
	waitKey(0);
	srand(3112);
	int i = rand()%img.size().height;
	int j = rand()%img.size().width;
	double b = img.at<Vec3b>(i, j)[0];
	double g = img.at<Vec3b>(i, j)[1];
	double r = img.at<Vec3b>(i, j)[2];
	cvtColor(img, hsv, COLOR_RGB2HSV);;

	std::cout << "For RGB" << endl;
	std::cout << "R:" << r << "  G:" << g << "  B:" << b << "\n" << endl;
	RGB2HSV(b/255, g/255, r/255);
	RGB2CMYK(b, g, r);
	RGB2RGBY(b, g, r);

	double h = hsv.at<Vec3b>(i, j)[0];
	double s = hsv.at<Vec3b>(i, j)[1];
	double v = hsv.at<Vec3b>(i, j)[2];
	std::cout << "For HSV using cvtColor" << endl;
	std::cout << "H:" << h << "  S:" << s << "  V:" << v << "\n" << endl;
	
	return 0;
}


void RGB2HSV(double b, double g, double r) {
	
	double v = max(max(b, g), r);
	double s, h;
	double v_ = v - min(min(b, g), r);

	if (v != 0)
		s = v_ / v;
	else
		s = 0;

	if (s == 0)
		h = 0;
	else {
		if (v == r)
			h = 60 * (g - b) / v_;
		else if (v == g)
			h = 120 + 60 * (b - r) / v_;
		else
			h = 240 + 60 * (r - g) / v_;

		if (h < 0)
			h = h + 360;
	}
	std::cout << "For HSV" << endl;
	std::cout << "H:" << h << "  S:" << s << "  V:" << v << "\n" << endl;
}


void RGB2CMYK(double b, double g, double r) {

	double k = 1 - max(max(r / 255, g / 255), b / 255);
	double c = (1 - k - r / 255) / (1 - k);
	double m = (1 - k - g / 255) / (1 - k);
	double y = (1 - k - b / 255) / (1 - k);

	std::cout << "For CMYK" << endl;
	std::cout << "C:" << c << "  M:" << m << \
		"  Y:" << y << "  K:" << k << "\n" << endl;
}


void RGB2RGBY(double b, double g, double r) {

	double r_ = r - (g + b) / 2;
	double g_ = g - (r + b) / 2;
	double b_ = b - (r + g) / 2;
	double y = (r + g) / 2 - abs(r - g) / 2 - b;

	std::cout << "For RGBY" << endl;
	std::cout << "R:" << r_ << "  G:" << g_ \
		<< "  B:" << b_ << "  Y:" << y << "\n" <<endl;
}


/*
Mat RGB2HSV(Mat pre_img, int height, int width) {
Mat img(height, width, CV_8UC3);
int i, j;
for (i = 0; i < height; i++) {
for (j = 0; j < width; j++) {
double b = pre_img.at<Vec3b>(i, j)[0];
double g = pre_img.at<Vec3b>(i, j)[1];
double r = pre_img.at<Vec3b>(i, j)[2];

double v = max(max(b, g), r);
double v_ = v - min(min(b, g), r);
double s, h;

if (v != 0)
s = v_ / v;
else
s = 0;

if (s == 0)
h = 0;
else{
if (v == r)
h = 60 * (g - b) / v_;
else if (v == g)
h = 120 + 60 * (b - r) / v_;
else
h = 240 + 60 * (r - g) / v_;

if(h < 0)
h = h + 360;
}

img.at<Vec3b>(i, j)[0] = h;
img.at<Vec3b>(i, j)[1] = s;
img.at<Vec3b>(i, j)[2] = v;
}
}
return img;
}*/