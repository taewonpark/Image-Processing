#include <opencv2\opencv.hpp>
#include <iostream>
#include <cmath>
#include "04.h"

using namespace cv;
using namespace std;


int main() {
	string imageName("C:\\Users\\samsung\\Desktop\\BEP\\BEP1\\lena.jpg");
	//string imageName("C:\\Users\\samsung\\Desktop\\BEP\\BEP5\\small.bmp");
	//string imageName("C:\\Users\\samsung\\Desktop\\BEP\\BEP5\\Rotation.jpg");
	Mat img = imread(imageName.c_str(), IMREAD_COLOR), gray, histogram;
	cvtColor(img, gray, COLOR_RGB2GRAY);

	int mask_size = 3;

	imshow("Original", gray); waitKey();

	
	/*
	//Dilation, Erosion, Opening, Closing
	imshow("Dilation", Dilation(gray, mask_size)); waitKey();
	imshow("Erosion", Erosion(gray, mask_size)); waitKey();
	imshow("Opening", Opening(gray, mask_size)); waitKey();
	imshow("Closing", Closing(gray, mask_size)); waitKey();
	*/


	/*
	//Hit or miss
	int C[] = { 1, 1, 0, 1, 1, 0, 0, 0, 0 };
	int D[] = { 0, 0, 0, 0, 0, 1, 0, 1, 1 };
	for (int i = 0; i < mask_size*mask_size; i++) {
		printf("%d ", C[i]);
		if (i%mask_size == mask_size-1) cout << endl;
	}
	cout << endl;
	for (int i = 0; i < mask_size*mask_size; i++) {
		printf("%d ", D[i]);
		if (i%mask_size == mask_size-1) cout << endl;
	}
	imshow("C", Erosion_mask(gray, C, mask_size)); waitKey();
	imshow("D", Erosion_mask(Complement(gray), D, mask_size)); waitKey();
	imshow("Hit or miss", Hit_or_miss(gray, mask_size, C, D)); waitKey();
	*/


	/*
	//Top hat
	imshow("Opening", Opening(gray, mask_size)); waitKey();
	imshow("Top hat", Top_hat(gray, mask_size)); waitKey();
	imshow("Closing", Closing(gray, mask_size)); waitKey();
	imshow("Bottom hat", Bottom_hat(gray, mask_size)); waitKey();
	*/


	/*
	//Grayscale
	//int mask[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	int mask[] = { -10, 1, 2, 10, -3, -7, -2, 9, 1 };
	for (int i = 0; i < mask_size*mask_size; i++) {
		printf("%3d ", mask[i]);
		if (i%mask_size == mask_size - 1) cout << endl << endl;
	}
	imshow("Dilation_gray", Dilation_gray(gray, mask, mask_size)); waitKey();
	imshow("Erosion_gray", Erosion_gray(gray, mask, mask_size)); waitKey();
	imshow("Opening_gray", Opening_gray(gray, mask, mask_size)); waitKey();
	imshow("Closing_gray", Closing_gray(gray, mask, mask_size)); waitKey();
	*/


	/*
	//Interpolation
	Mat cv_nearest, cv_linear, cv_cubic;
	resize(gray, cv_nearest, cv_nearest.size(), 10, 10, INTER_NEAREST);
	resize(gray, cv_linear, cv_linear.size(), 10, 10, INTER_LINEAR);
	resize(gray, cv_cubic, cv_cubic.size(), 10, 10, INTER_CUBIC);

	imshow("Original", gray); waitKey();

	imshow("cv nearest", cv_nearest); waitKey();
	imshow("Nearest interpolation", Nearest_interpolation(gray, 10)); waitKey();

	imshow("cv linear", cv_linear); waitKey();
	imshow("Bilinear interpolation", Bilinear_interpolation(gray, 10)); waitKey();

	imshow("cv cubic", cv_cubic); waitKey();
	imshow("Bicubic interpolation", Bicubic_interpolation(gray, 10, -0.75)); waitKey();
	imshow("Bicubic interpolation x next y", Bicubic_interpolation_x_next_y(gray, 10, -0.75)); waitKey();

	imshow("Bispline interpolation", Bispline_interpolation(gray, 10)); waitKey();
	imshow("Bispline interpolation x next y", Bispline_interpolation_x_next_y(gray, 30)); waitKey();
	*/


	/*
	//Rotation
	imshow("Rotation", Rotation_with_interpolation(gray, PI / 4)); waitKey();
	imshow("Rotation with cut", Rotation_with_cut(gray, PI / 4)); waitKey();
	imshow("Rotation with interpolation", Rotation_with_interpolation(gray, PI / 4)); waitKey();
	*/
	
	return 0;
}