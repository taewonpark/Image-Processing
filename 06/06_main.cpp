#include <opencv2\opencv.hpp>
#include <iostream>
#include "06.h"

using namespace cv;
using namespace std;


int main() {
	string imageName("C:\\Users\\samsung\\Desktop\\BEP\\BEP6\\small.jpg");
	Mat img = imread(imageName.c_str(), IMREAD_COLOR), gray;
	cvtColor(img, gray, COLOR_RGB2GRAY);

	double **realImg = new double*[img.rows], **imagImg = new double*[img.rows];

	for (int i = 0; i < img.rows; i++) {
		realImg[i] = new double[img.cols];
		imagImg[i] = new double[img.cols];
		memset(realImg[i], 0, sizeof(double)*img.cols);
		memset(imagImg[i], 0, sizeof(double)*img.cols);
	}

	imshow("Original", gray); waitKey();
	imshow("DFT", DFT(gray, realImg, imagImg)); waitKey();
	imshow("IDFT", IDFT(realImg, imagImg, gray.rows, gray.cols)); waitKey();

	/*
	//add Noise and remove it
	Scalar Black(0, 0, 0);
	int freq = 40, interval = gray.cols / freq, bound = gray.cols / 2 / interval + 1;
	int center_row = (img.rows - 1) / 2, center_col = (img.cols - 1) / 2;
	for (int i = 0; i < bound; i++) {
		line(gray, Point(center_col + interval * i, 0), Point(center_col + interval * i, img.rows - 1), Black, 1);
		line(gray, Point(center_col - interval * i, 0), Point(center_col - interval * i, img.rows - 1), Black, 1);
	}
	imshow("Noise", gray); waitKey();
	imshow("DFT", DFT(gray, realImg, imagImg)); waitKey();
	for (int i = 1; i < (interval + 1) / 2; i++)
		for (int m = -6; m < 7; m++)
			for (int n = -5; n < 6; n++){
					realImg[center_row + m][center_col + i*freq - n] = 0;
					imagImg[center_row + m][center_col + i*freq - n] = 0;
					realImg[center_row + m][center_col - i*freq + n] = 0;
					imagImg[center_row + m][center_col - i*freq + n] = 0;
				}
	imshow("remove Noise", Spectrum(realImg, imagImg, gray.rows, gray.cols)); waitKey();
	imshow("IDFT", IDFT(realImg, imagImg, gray.rows, gray.cols)); waitKey();
	*/

	return 0;
}