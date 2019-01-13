#include <opencv2\opencv.hpp>
#include <iostream>
#include <cmath>
#include "03.h"

using namespace cv;
using namespace std;

const float PI = 3.14159265358979323846f;


Mat zero_padding(Mat img, int pad) {
	Mat padded(img.rows + 2 * pad, img.cols + 2 * pad, CV_8UC1, Scalar(0, 0, 0));
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			padded.at<uchar>(i + pad, j + pad) = img.at<uchar>(i, j);
		}
	}
	return padded;
}



Mat edge_padding(Mat img, int pad) {
	Mat padded = zero_padding(img, pad);
	for (int i = 0; i < pad; i++) {
		for (int j = 0; j < img.rows; j++) {
			padded.at<uchar>(j + pad, i) = img.at<uchar>(j, 1);
			padded.at<uchar>(j + pad, padded.cols - i - 1) = img.at<uchar>(j, img.cols - 2);
		}
		for (int j = 0; j < img.cols; j++) {
			padded.at<uchar>(i, j + pad) = img.at<uchar>(1, j);
			padded.at<uchar>(padded.rows - i - 1, j + pad) = img.at<uchar>(img.rows - 2, j);
		}
		for (int j = 0; j < pad; j++) {
			padded.at<uchar>(i, j) = img.at<uchar>(1, 1);
			padded.at<uchar>(padded.rows - i - 1, j) = img.at<uchar>(img.rows - 2, 1);
			padded.at<uchar>(i, padded.cols - i) = img.at<uchar>(1, img.cols - 2);
			padded.at<uchar>(padded.rows - i - 1, padded.cols - i - 1) = img.at<uchar>(img.rows - 2, img.cols - 2);
		}
	}
	return padded;
}



Mat reflect_padding(Mat img, int pad) {
	Mat padded = zero_padding(img, pad);
	for (int i = 0; i < pad; i++) {
		for (int j = 0; j < img.rows; j++) {
			padded.at<uchar>(j + pad, i) = img.at<uchar>(j, pad - i);
			padded.at<uchar>(j + pad, padded.cols - i - 1) = img.at<uchar>(j, img.cols - pad + i - 1);
		}
		for (int j = 0; j < img.cols; j++) {
			padded.at<uchar>(i, j + pad) = img.at<uchar>(pad - i, j);
			padded.at<uchar>(padded.rows - i - 1, j + pad) = img.at<uchar>(img.rows - pad + i - 1, j);
		}
		for (int j = 0; j < pad; j++) {
			padded.at<uchar>(i, j) = img.at<uchar>(pad - i, pad - j);
			padded.at<uchar>(padded.rows - i - 1, j) = img.at<uchar>(img.rows - pad + i - 1, pad - j);
			padded.at<uchar>(i, padded.cols - i) = img.at<uchar>(pad - i, img.rows - pad + i - 1);
			padded.at<uchar>(padded.rows - i - 1, padded.cols - i - 1) = img.at<uchar>(img.rows - pad + i - 1, img.rows - pad + i - 1);
		}
	}
	return padded;
}



void Seletionsort(int* arr, int size) {
	int min_val, check, temp;

	for (int i = 0; i < size; i++) {
		min_val = arr[i]; check = i;
		for (int j = i + 1; j < size; j++)
			if (min_val > arr[j]) {
				min_val = arr[j];
				check = j;
			}
		temp = arr[i];
		arr[i] = arr[check];
		arr[check] = temp;
	}
}



void Quicksort(int* arr, int left, int right) {
	int i = left, j = right;
	int pivot = arr[(left + right) / 2];
	int temp;

	while (1) {
		while (arr[i] < pivot) i++;
		while (arr[j] > pivot) j--;

		if (i > j) break;

		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
		i++; j--;
	}

	if (left < j) Quicksort(arr, left, j);
	if (i < right) Quicksort(arr, i, right);
}



Mat Average_Filter(Mat img, int height, int width) {
	Mat result(height, width, CV_8UC1, Scalar(0, 0, 0));
	int num = (img.rows - result.rows) + 1;
	int total_num = num*num;
	int mean;

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			mean = 0;

			for (int m = 0; m < num; m++)
				for (int n = 0; n < num; n++) {
					mean += img.at<uchar>(i + m, j + n)/total_num;
				}

			result.at<uchar>(i, j) = mean;
		}
	return result;
}



Mat Median_Filter(Mat img, int height, int width) {
	Mat result(height, width, CV_8UC1, Scalar(0, 0, 0));
	int num = (img.rows - result.rows) + 1;
	int total_num = num*num;
	int* arr = new int[total_num];
	int check;

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			check = 0;
			for (int m = 0; m < num; m++)
				for (int n = 0; n < num; n++) {
					 arr[check] = img.at<uchar>(i + m, j + n);
					 check += 1;
				}
			Quicksort(arr, 0, total_num - 1);
			result.at<uchar>(i, j) = arr[(total_num + 1)/2];
		}
	delete[] arr;

	return result;
}



Mat Maximum_Filter(Mat img, int height, int width) {
	Mat result(height, width, CV_8UC1, Scalar(0, 0, 0));
	int num = (img.rows - result.rows) + 1;
	int total_num = num*num;
	int* arr = new int[total_num];
	int check;

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			check = 0;
			for (int m = 0; m < num; m++)
				for (int n = 0; n < num; n++) {
					arr[check] = img.at<uchar>(i + m, j + n);
					check += 1;
				}
			Quicksort(arr, 0, total_num - 1);
			result.at<uchar>(i, j) = arr[total_num - 1];
		}
	delete[] arr;

	return result;
}



Mat Minimum_Filter(Mat img, int height, int width) {
	Mat result(height, width, CV_8UC1, Scalar(0, 0, 0));
	int num = (img.rows - result.rows) + 1;
	int total_num = num*num;
	int* arr = new int[total_num];
	int check;

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			check = 0;
			for (int m = 0; m < num; m++)
				for (int n = 0; n < num; n++) {
					arr[check] = img.at<uchar>(i + m, j + n);
					check += 1;
				}
			Quicksort(arr, 0, total_num - 1);
			result.at<uchar>(i, j) = arr[0];
		}
	delete[] arr;

	return result;
}



Mat Midpoint_Filter(Mat img, int height, int width) {
	Mat result(height, width, CV_8UC1, Scalar(0, 0, 0));
	int num = (img.rows - result.rows) + 1;
	int total_num = num*num;
	int* arr = new int[total_num];
	int check;

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			check = 0;
			for (int m = 0; m < num; m++)
				for (int n = 0; n < num; n++) {
					arr[check] = img.at<uchar>(i + m, j + n);
					check += 1;
				}
			Quicksort(arr, 0, total_num - 1);
			result.at<uchar>(i, j) = (arr[0] + arr[total_num - 1])/2;
		}
	delete[] arr;

	return result;
}



Mat Alpha_Filter(Mat img, int height, int width) {
	Mat result(height, width, CV_8UC1, Scalar(0, 0, 0));
	int num = (img.rows - result.rows) + 1;
	int total_num = num*num;
	int* arr = new int[total_num];
	int check, summation;

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			check = 0; summation = 0;
			for (int m = 0; m < num; m++)
				for (int n = 0; n < num; n++) {
					arr[check] = img.at<uchar>(i + m, j + n);
					check += 1;
				}
			Quicksort(arr, 0, total_num - 1);

			for (int k = 1; k < total_num - 1; k++) summation += arr[k];
			result.at<uchar>(i, j) = summation / (total_num - 2);
		}
	delete[] arr;

	return result;
}



Mat Convolution(Mat img, int* mask, int height, int width) {
	// mask [[ 1st row vec ], [ 2nd row vec ], ... ]
	// but one dimension array
	Mat result(height, width, CV_8UC1, Scalar(0, 0, 0));
	int num = (img.rows - result.rows) + 1;
	int total_num = num*num;
	int summation;

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			summation = 0;
			for (int m = 0; m < num; m++)
				for (int n = 0; n < num; n++) {
					summation += img.at<uchar>(i + m, j + n) * mask[num*m + n];
				}
			summation = (int)abs(summation);
			if (summation < 0) result.at<uchar>(i, j) = 0;
			else if (summation > 255) result.at<uchar>(i, j) = 255;
			else result.at<uchar>(i, j) = summation;
		}
	
	return result;
}



Mat Convolution_F(Mat img, float* mask, int height, int width) {
	// mask [[ 1st row vec ], [ 2nd row vec ], ... ]
	// but one dimension array
	Mat result(height, width, CV_8UC1, Scalar(0, 0, 0));
	int num = (img.rows - result.rows) + 1;
	int total_num = num*num, fin;
	float summation;

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			summation = 0;
			for (int m = 0; m < num; m++)
				for (int n = 0; n < num; n++) {
					summation += (float)img.at<uchar>(i + m, j + n) * mask[num*m + n];
				}
			fin = (int)abs(summation);
			if (fin < 0) result.at<uchar>(i, j) = 0;
			else if (fin > 255) result.at<uchar>(i, j) = 255;
			else result.at<uchar>(i, j) = (int)abs(summation);
		}

	return result;
}



Mat RSS(Mat x_img, Mat y_img) {
	Mat result = x_img.clone();
	int x, y;

	for (int i = 0; i<result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			x = x_img.at<uchar>(i, j);
			y = y_img.at<uchar>(i, j);
			result.at<uchar>(i, j) = sqrt(x*x + y*y);
		}
	return result;
}



void Get_sobel(float* x_mask, float* y_mask, int size) {
	int half = size / 2;
	int total = 0;
	for (int i = 0; i < half; i++) total += 2 * (half + 1) *(i + 1);
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size + 1; j++) {
			((i - half) == 0) ? x_mask[i*size + j] = (float)2 * (j - half) / total : x_mask[i*size + j] = (float)(j - half) / total;
			((j - half) == 0) ? y_mask[i*size + j] = (float)2 * (i - half) / total : y_mask[i*size + j] = (float)(i - half) / total;
		}
}



void Get_prewitt(float* x_mask, float* y_mask, int size) {
	int half = size / 2;
	int total = 0;
	for (int i = 0; i < half; i++) total += (2 * half + 1) *(i + 1);
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size + 1; j++) {
			x_mask[i*size + j] = (float)(j - half) / total;
			y_mask[i*size + j] = (float)(i - half) / total;
		}
}



void Get_roberts(float* x_mask, float* y_mask, int size) {
	int half = size / 2;

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size + 1; j++) {
			x_mask[i*size + j] = 0;
			y_mask[i*size + j] = 0;
		}

	x_mask[half*size + half] = 1;
	x_mask[0*size + 0] = -1;
	y_mask[half*size + half] = 1;
	y_mask[0*size + size - 1] = -1;
}



void Get_laplacian(float *without, float* with, int size) {
	int half = size / 2;

	for (int i = 0; i < size*size; i++) {
		with[i] = (float)1 / 8;
		without[i] = 0;
	}
	with[size*half + half] = -1;

	for (int i = 1; i < size*size; i += 2) without[i] = (float)1/4;
	without[size*half + half] = -1;
}



void Get_gaussian(float *mask, int size, float stddev) {
	int half = size / 2;
	float variance = stddev * stddev;
	float summation = 0;
	int x, y;

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size + 1; j++) {
			x = j - half; y = i - half;
			mask[size*i + j] = exp(-(x*x + y*y) / (2 * variance)) / sqrt(2 * PI*variance);
			summation += mask[size*i + j];
		}
	for (int i = 0; i < size*size; i++) mask[i] = mask[i] / summation;
}



void Get_LoG(float *mask, int size, float stddev) {
	int half = size / 2;
	float variance = stddev * stddev;
	float summation = 0;
	int x, y;

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size + 1; j++) {
			x = j - half; y = i - half;
			mask[size*i + j] = -exp(-(x*x + y*y) / (2 * variance)) * (1 - (x*x + y*y) / (2 * variance) / (PI*variance*variance));
			summation += mask[size*i + j];
		}
	for (int i = 0; i < size*size; i++) mask[i] = mask[i] / summation;
}



void Get_DoG(float *mask, int size, float stddev1, float stddev2) {
	int half = size / 2;
	float *mask1 = new float[size*size], *mask2 = new float[size*size];

	Get_gaussian(mask1, size, stddev1);
	Get_gaussian(mask2, size, stddev2);

	for (int i = 0; i < size*size; i++) {
		mask[i] = mask1[i] - mask2[i];
	}
	delete[] mask1, mask2;
	for (int i = 0; i < size*size; i++) mask[i] = mask[i];
}



void Get_Highboost(float* x_mask, float* y_mask, int size, double c) {
	int half = size*size / 2;
	
	Get_laplacian(x_mask, y_mask, size);
	for (int i = 0; i < size*size; i++) {
		x_mask[i] *= -1;
		y_mask[i] *= -1;
	}
	x_mask[half] += c;
	y_mask[half] += c;
}



Mat adding(Mat x_img, Mat y_img, float alpha) {
	Mat result = x_img.clone();
	float x, y, temp;

	for (int i = 0; i<result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			x = x_img.at<uchar>(i, j);
			y = y_img.at<uchar>(i, j);
			temp = x + alpha*y;
			(temp > 255) ? result.at<uchar>(i, j) = 255 : \
				result.at<uchar>(i, j) = temp;
		}
	return result;
}

