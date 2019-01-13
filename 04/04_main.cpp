#include <opencv2\opencv.hpp>
#include <iostream>
#include <cmath>
#include "04.h"

using namespace cv;
using namespace std;


int main() {
	string imageName("C:\\Users\\samsung\\Desktop\\BEP\\BEP1\\lena.jpg");
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
	
	return 0;
}