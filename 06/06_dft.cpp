#include <opencv2\opencv.hpp>
#include <iostream>
#include <math.h>
#include "06.h"

using namespace cv;
using namespace std;

const float PI = 3.14159265358979323846f;


Mat DFT(Mat img, double **realImg, double **imagImg) {
	int H = img.rows, W = img.cols, max_val = 0;
	Mat result(H, W, CV_8UC1, Scalar(0, 0, 0));
	double real, imagine, temp;
	double **arr = new double*[img.rows];
	int center_row = (img.rows - 1) / 2, center_col = (img.cols - 1) / 2;

	for (int i = 0; i < img.rows; i++) {
		arr[i] = new double[img.cols];
		memset(arr[i], 0, sizeof(double)*img.cols);
	}

	for (int i = 0; i < H; i++) {
		cout << i << endl << endl;
		for (int j = 0; j < W; j++) {
			real = 0; imagine = 0;
			for (int m = 0; m < H; m++)
				for (int n = 0; n < W; n++) {
					real += (double)img.at<uchar>(m, n) * cos(2 * PI*((double)m*i / H + (double)n*j / W));
					imagine -= (double)img.at<uchar>(m, n) * sin(2 * PI*((double)m*i / H + (double)n*j / W));
				}
			realImg[(i + center_row) % H][(j + center_col) % W] = real;
			imagImg[(i + center_row) % H][(j + center_col) % W] = imagine;
			temp = log2(sqrt(real*real + imagine*imagine) + 1);
			if (temp < 0) temp = 0.0;
			else if (temp > 255) temp = 255.0;
			if (temp > max_val) max_val = (int)temp;
			arr[(i + center_row) % H][(j + center_col) % W] = temp;
		}
	}
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++) {
			result.at<uchar>(i, j) = (int)arr[i][j] * 255 / max_val;
		}
	return result;
}



Mat IDFT(double **realImg, double **imagImg, int H, int W) {
	Mat result(H, W, CV_8UC1, Scalar(0, 0, 0));
	double real, imagine, angle;
	int temp, center_row = (H - 1) / 2, center_col = (W - 1) / 2;

	for (int i = 0; i < H; i++) {
		cout << i << endl << endl;
		for (int j = 0; j < W; j++) {
			real = 0; imagine = 0;
			for (int m = 0; m < H; m++)
				for (int n = 0; n < W; n++) {
					angle = 2 * PI*((double)m*i / H + (double)n*j / W);
					real += (realImg[m][n] * cos(angle) - imagImg[m][n] * sin(angle));
					imagine += (realImg[m][n] * sin(angle) + imagImg[m][n] * cos(angle));
				}
			temp = (int)sqrt(real*real + imagine*imagine) / (H*W);
			if (temp < 0) temp = 0;
			else if (temp > 255) temp = 255;
			result.at<uchar>(i, j) = temp;
		}
	}
	return result;
}



Mat Spectrum(double **realImg, double **imagImg, int H, int W) {
	Mat result(H, W, CV_8UC1, Scalar(0, 0, 0));
	int temp;

	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++) {
			temp = sqrt(realImg[i][j] * realImg[i][j] + imagImg[i][j] * imagImg[i][j]);
			if (temp > 255) temp = 255;
			result.at<uchar>(i, j) = temp;
		}
	return result;
}