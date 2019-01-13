#include <opencv2\opencv.hpp>
#include <iostream>
#include <cmath>
#include "03.h"

using namespace cv;
using namespace std;


int main() {
	string imageName("C:\\Users\\samsung\\Desktop\\BEP\\BEP3\\woman.bmp");
	//string imageName("C:\\Users\\samsung\\Desktop\\BEP\\BEP3\\cat.jpg");
	Mat img = imread(imageName.c_str(), IMREAD_COLOR), gray, histogram;
	cvtColor(img, gray, COLOR_RGB2GRAY);

	int padding_size = 1;
	int k_size = 2 * padding_size + 1;

	


	/*
	//sort
	int a[] = { 4, 2, 1, 3 };
	Seletionsort(a, sizeof(a) / sizeof(*a));
	Quicksort(a, 0, sizeof(a) / sizeof(*a) - 1);
	*/


	/*
	//padding
	Mat padded1 = zero_padding(gray, padding_size);
	Mat padded2 = edge_padding(gray, padding_size);
	Mat padded3 = reflect_padding(gray, padding_size);
	imshow("s", gray); waitKey();
	imshow("a", padded1); waitKey();
	imshow("b", padded2); waitKey();
	imshow("c", padded3); waitKey();
	*/


	/*
	//filter
	Mat padded_1 = zero_padding(gray, 1);
	Mat padded_2 = zero_padding(gray, 5);
	imshow("a", Average_Filter(padded_1, gray.rows, gray.cols));
	waitKey();
	imshow("b", Average_Filter(padded_2, gray.rows, gray.cols));
	waitKey();
	*/


	/*
	//need two int type mask (error;; just use float)
	float *x_mask = new float[k_size * k_size], *y_mask = new float[k_size * k_size];
	Get_laplacian(x_mask, y_mask, k_size);
	for (int i = 0; i < k_size*k_size; i++) {
	printf("%3.2f ", x_mask[i]);
	if (i%k_size == k_size - 1) cout << endl;
	}
	cout << endl;
	for (int i = 0; i < k_size*k_size; i++) {
	printf("%3.2f ", y_mask[i]);
	if (i%k_size == k_size - 1) cout << endl;
	}

	Mat padded = zero_padding(gray, padding_size);
	Mat x_img = Convolution_F(padded, x_mask, gray.rows, gray.cols);
	Mat y_img = Convolution_F(padded, y_mask, gray.rows, gray.cols);
	imshow("x", x_img); waitKey();
	imshow("y", y_img); waitKey();
	//imshow("S", RSS(x_img, y_img)); waitKey();
	*/
	
	
	/*
	//LoG
	float *mask1 = new float[k_size * k_size];
	float *mask2 = new float[k_size * k_size];
	float *mask3 = new float[k_size * k_size];

	Mat padded = zero_padding(gray, padding_size);
	Get_gaussian(mask1, k_size, 1.6);
	Mat gauss = zero_padding(Convolution_F(padded, mask1, gray.rows, gray.cols), padding_size);
	imshow("g", gauss); waitKey();
	Get_laplacian(mask2, mask3, k_size);
	imshow("a", Convolution_F(gauss, mask2, gray.rows, gray.cols)); waitKey();
	imshow("b", Convolution_F(gauss, mask3, gray.rows, gray.cols)); waitKey();
	*/


	/*
	//DoG
	float *mask = new float[k_size * k_size];
	Mat padded = zero_padding(gray, padding_size);
	Get_DoG(mask, k_size, 0.4, 0.8);
	for (int i = 0; i < k_size*k_size; i++) {
		printf("%3.2f ", mask[i]);
		if (i%k_size == k_size - 1) cout << endl << endl;
	}
	imshow("b", Convolution_F(padded, mask, gray.rows, gray.cols)); waitKey();
	*/

	/*
	//highboost
	int *x_mask = new int[k_size * k_size], *y_mask = new int[k_size * k_size];
	Get_laplacian(x_mask, y_mask, k_size);
	imshow("g", gray); waitKey();
	Mat padded = zero_padding(gray, padding_size);
	Mat x_img = Convolution(padded, x_mask, gray.rows, gray.cols);
	Mat y_img = Convolution(padded, y_mask, gray.rows, gray.cols);
	imshow("x", adding(gray, x_img, 0.2)); waitKey();
	imshow("y", adding(gray, y_img, 0.2)); waitKey();
	*/
	

	return 0;
}