#include <opencv2\opencv.hpp>
#include <iostream>
#include <cmath>
#include "02.h"

using namespace cv;
using namespace std;

int *check_freq(Mat img, int histSize) {

	int *result = new int[histSize];
	for (int i = 0; i < histSize; i++) result[i] = 0;

	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			result[(int)img.at<uchar>(i, j)] += 1;
		}
	}
	return result;
}



int check_height(int* freq, int histSize) {
	int result = 0;
	for (int i = 0; i < histSize; i++) {
		result = max(freq[i], result);
	}
	return result;
}



void PlotPDF(int* freq, int h, int bin_w, double bin_h, int histSize) {
	Scalar black(0, 0, 0);
	Mat histImg(h, (int)bin_w*histSize, CV_8UC1, Scalar(255, 255, 255));

	for (int i = 0; i < histSize; i++) {
		for (int j = 0; j < bin_w; j++) {
			line(histImg, Point(i*bin_w + j, h), Point(i*bin_w + j, h - cvRound((double)bin_h*freq[i])), black);
		}
	}
	imshow("Histogram", histImg);
}



void PlotCDF(int* freq, int h, int bin_w, double bin_h, int histSize) {
	Scalar black(0, 0, 0);
	Mat histImg(h, (int)bin_w*histSize, CV_8UC1, Scalar(255, 255, 255));
	float summation = 0.0;

	for (int i = 0; i < histSize; i++) {
		summation += freq[i];
		for (int j = 0; j < bin_w; j++) {
			line(histImg, Point(i*bin_w + j, h), Point(i*bin_w + j, h - cvFloor(bin_h*summation)), black);
		}
	}
	imshow("Histogram", histImg);
}



void Sliding(Mat img, Mat& result, int offset) {
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			uchar* pixel = result.ptr<uchar>(i) + j;
			int temp = (int)img.at<uchar>(i, j) + offset;
			if (temp < 0) {
				*pixel = 0;
			}
			else if (temp > 255) {
				*pixel = 255;
			}
			else {
				*pixel = temp;
			}
		}
	}
}



void Stretching(Mat img, Mat& result, int max_val, int min_val) {
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			uchar* pixel = result.ptr<uchar>(i) + j;
			*pixel = cvRound(((int)img.at<uchar>(i, j) - min_val) * 255 / (max_val - min_val));
		}
	}
}



void Shrink(Mat img, Mat& result, int max_val, int min_val, int s_max, int s_min) {
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			uchar* pixel = result.ptr<uchar>(i) + j;
			*pixel = cvRound(((int)img.at<uchar>(i, j) - min_val) * (s_max - s_min) / (max_val - min_val) + s_min);
		}
	}
}



void Equalize(Mat img, Mat& result, int total, int* freq, int histSize, int max_val, int min_val) {
	int summation = freq[min_val];
	int *value = new int[histSize];
	double slope = (double)total / histSize;

	for (int i = 0; i < histSize; i++) (i <= min_val) ? value[i] = 0 : value[i] = histSize - 1;

	for (int i = 0, j = min_val; i < histSize; i++) {
		if (summation < (int)(slope*i)) {
			while (summation < (int)(slope*i)) {
				value[j] = i;
				j += 1;
				summation += freq[j];
			}
		}
	}

	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			uchar* pixel = result.ptr<uchar>(i) + j;
			*pixel = value[(int)img.at<uchar>(i, j)];
		}
	}
}



void Quantization(Mat img, Mat& result, int L, int histSize) {

	int* value = new int[L];
	for (int i = 0; i < L; i++) value[i] = (histSize - 1) * i / (L - 1);

	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			uchar* pixel = result.ptr<uchar>(i) + j;
			*pixel = value[(int)img.at<uchar>(i, j) / L];
		}
	}
}



void Static_bin(Mat img, Mat& result, int T, int* check) {
	*check = T;
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			uchar* pixel = result.ptr<uchar>(i) + j;
			if ((int)img.at<uchar>(i, j) > T) {
				*pixel = 255;
			}else {
				*pixel = 0;
			}
		}
	}
}



void P_tile_bin(Mat img, Mat& result, int P, int* freq, int total, int histSize, int* check) {
	int threshold;
	int ckpt = total / P, summation = 0;

	for (int i = 0; i < histSize; i++) {
		summation += freq[i];
		if (summation >= ckpt) {
			threshold = i;
			break;
		}
	}
	*check = threshold;
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			uchar* pixel = result.ptr<uchar>(i) + j;
			if ((int)img.at<uchar>(i, j) > threshold) {
				*pixel = 255;
			}
			else {
				*pixel = 0;
			}
		}
	}
}



void Iter_bin(Mat img, Mat& result, int threshold, int eps, int* freq, int histSize, int* check) {
	int T = threshold, T_ = 0;
	double summation, mean1, mean2, local_total;
	
	while (T - T_ > eps) {
		summation = 0, local_total = 0;
		for (int i = 0; i < T; i++) {
			summation += freq[i]*i;
			local_total += freq[i];
		}
		(local_total == 0) ? mean1 = 53489389 : mean1 = (double)summation / local_total;

		summation = 0, local_total = 0;
		for (int i = histSize - 1; i >= T; i--) {
			summation += freq[i]*i;
			local_total += freq[i];
		}
		(local_total == 0) ? mean2 = 53489389 : mean2 = (double)summation / local_total;

		T_ = T;
		T = (int)(mean1 + mean2) / 2;
	}
	*check = T;
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			uchar* pixel = result.ptr<uchar>(i) + j;
			if ((int)img.at<uchar>(i, j) > T) {
				*pixel = 255;
			}
			else {
				*pixel = 0;
			}
		}
	}
}



void Otsu(Mat img, Mat& result, int total, int* freq, int histSize, int max_val, int min_val, int* check) {
	double* value = new double[histSize];
	double check_val=0, val, local_total, summation, mean0, mean1, dev0, dev1, w_0, w_1;
	int threshold;

	for (int i = min_val + 1 ; i < max_val-1; i++) {
		mean0 = 0, mean1 = 0, dev0 = 0, dev1 = 0;
		summation = 0, local_total = 0;
		for (int j = 0; j < i; j++) {
			summation += freq[j]*j;
			local_total += freq[j];
		}
		mean0 = (double)summation / local_total;

		for (int j = 0; j < i; j++) dev0 += (freq[j] - mean0)*(freq[j] - mean0);
		dev0 /= local_total;

		w_0 = local_total / total;
		w_1 = 1 - w_0;

		summation = 0, local_total = 0;
		for (int j = histSize - 1; j >= i; j--) {
			summation += freq[j] * j;
			local_total += freq[j];
		}
		mean1 = (double)summation / local_total;

		for (int j = histSize - 1; j >= i; j--) dev1 += (freq[j] - mean1)*(freq[j] - mean1);
		dev1 /= local_total;

		val = w_0*w_1*(mean0 - mean1)*(mean0 - mean1);

		if (val >= check_val) {
			check_val = val;
			threshold = i;
		}
	}
	*check = threshold;
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			uchar* pixel = result.ptr<uchar>(i) + j;
			if ((int)img.at<uchar>(i, j) > threshold) {
				*pixel = 255;
			}
			else {
				*pixel = 0;
			}
		}
	}
	cout << threshold << endl;
}



void Valley(Mat img, Mat& result, int total, int* freq, int histSize, int max_val, int min_val, int* check) {
	double* value = new double[histSize];
	double check_val = 0, val, local_total, summation, mean0, mean1, dev0, dev1, w_0, w_1;
	int threshold;

	for (int i = min_val + 1; i < max_val - 1; i++) {
		mean0 = 0, mean1 = 0, dev0 = 0, dev1 = 0;
		summation = 0, local_total = 0;
		for (int j = 0; j < i; j++) {
			summation += freq[j] * j;
			local_total += freq[j];
		}
		mean0 = (double)summation / local_total;

		for (int j = 0; j < i; j++) dev0 += (freq[j] - mean0)*(freq[j] - mean0);
		dev0 /= local_total;

		w_0 = local_total / total;
		w_1 = 1 - w_0;

		summation = 0, local_total = 0;
		for (int j = histSize - 1; j >= i; j--) {
			summation += freq[j] * j;
			local_total += freq[j];
		}
		mean1 = (double)summation / local_total;

		for (int j = histSize - 1; j >= i; j--) dev1 += (freq[j] - mean1)*(freq[j] - mean1);
		dev1 /= local_total;

		val = w_0*w_1*(mean0 - mean1)*(mean0 - mean1) * (1 - (double)freq[i] / total);

		if (val >= check_val) {
			check_val = val;
			threshold = i;
		}
	}
	*check = threshold;
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			uchar* pixel = result.ptr<uchar>(i) + j;
			if ((int)img.at<uchar>(i, j) > threshold) {
				*pixel = 255;
			}
			else {
				*pixel = 0;
			}
		}
	}
	cout << threshold << endl;
}



void PlotPDF_threshold(int* freq, int h, int bin_w, double bin_h, int histSize, int threshold) {
	Scalar some(127, 127, 127), black(0, 0, 0);
	Mat histImg(h, (int)bin_w*histSize, CV_8UC1, Scalar(255, 255, 255));

	for (int i = 0; i < histSize; i++) {
		for (int j = 0; j < bin_w; j++) {
			line(histImg, Point(i*bin_w + j, h), Point(i*bin_w + j, h - cvRound((double)bin_h*freq[i])), black);
		}
	}
	for (int i = 0; i < bin_w; i++) line(histImg, Point(threshold*bin_w + i, h), Point(threshold*bin_w + i, 0), some);
	imshow("Histogram", histImg);
}




void PlotCDF_threshold(int* freq, int h, int bin_w, double bin_h, int histSize, int threshold) {
	Scalar black(0, 0, 0);
	Mat histImg(h, (int)bin_w*histSize, CV_8UC1, Scalar(255, 255, 255));
	float summation = 0.0;

	for (int i = 0; i < histSize; i++) {
		summation += freq[i];
		for (int j = 0; j < bin_w; j++) {
			line(histImg, Point(i*bin_w + j, h), Point(i*bin_w + j, h - cvFloor(bin_h*summation)), black);
		}
	}
	for (int i = 0; i < bin_w; i++) line(histImg, Point(threshold*bin_w + i, h), Point(threshold*bin_w + i, 0), black);
	imshow("Histogram", histImg);
}