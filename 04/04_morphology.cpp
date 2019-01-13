#include <opencv2\opencv.hpp>
#include <iostream>
#include <cmath>
#include "..\03\03.h"
#include "04.h"

using namespace cv;
using namespace std;


Mat Dilation(Mat img, int mask_size) {
	int pad = mask_size / 2, check;
	Mat padded = zero_padding(img, pad);
	Mat result(img.rows, img.cols, CV_8UC1, Scalar(0, 0, 0));

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			check = 0;
			for (int m = 0; m < mask_size; m++)
				for (int n = 0; n < mask_size; n++) {
					check = check || padded.at<uchar>(i + m, j + n);
				}
			result.at<uchar>(i, j) = check * 255;
		}
	return result;
}



Mat Erosion(Mat img, int mask_size) {
	int pad = mask_size / 2, check;
	Mat padded = zero_padding(img, pad);
	Mat result(img.rows, img.cols, CV_8UC1, Scalar(0, 0, 0));

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			check = 1;
			for (int m = 0; m < mask_size; m++)
				for (int n = 0; n < mask_size; n++) {
					check = check && padded.at<uchar>(i + m, j + n);
				}
			result.at<uchar>(i, j) = check * 255;
		}
	return result;
}



Mat Opening(Mat img, int mask_size) {
	Mat result = Erosion(img, mask_size);
	result = Dilation(result, mask_size);
	return result;
}



Mat Closing(Mat img, int mask_size) {
	Mat result = Dilation(img, mask_size);
	result = Erosion(result, mask_size);
	return result;
}



Mat Complement(Mat img) {
	Mat result(img.rows, img.cols, CV_8UC1, Scalar(0, 0, 0));

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			result.at<uchar>(i, j) = 255 - img.at<uchar>(i, j);
		}
	return result;
}



Mat Image_and(Mat img1, Mat img2) {
	Mat result(img1.rows, img1.cols, CV_8UC1, Scalar(0, 0, 0));

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			result.at<uchar>(i, j) = 255 * (img1.at<uchar>(i, j) && img2.at<uchar>(i, j));
		}
	return result;
}



Mat Erosion_mask(Mat img, int *mask, int mask_size) {
	int pad = mask_size / 2, check;
	Mat padded = zero_padding(img, pad);
	Mat result(img.rows, img.cols, CV_8UC1, Scalar(0, 0, 0));

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			check = 1;
			for (int m = 0; m < mask_size; m++)
				for (int n = 0; n < mask_size; n++) {
					if (mask[mask_size*m + n])
						check = check && padded.at<uchar>(i + m, j + n);
				}
			result.at<uchar>(i, j) = check * 255;
		}
	return result;
}



Mat Hit_or_miss(Mat img, int mask_size, int *C, int *D) {
	Mat temp = Complement(img);
	Mat result = Image_and(Erosion_mask(img, C, mask_size), Erosion_mask(temp, D, mask_size));
	return result;
}



Mat Top_hat(Mat img, int mask_size) {
	Mat result(img.rows, img.cols, CV_8UC1, Scalar(0, 0, 0));
	Mat temp = Opening(img, mask_size);

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			result.at<uchar>(i, j) = img.at<uchar>(i, j) - temp.at<uchar>(i, j);
		}
	return result;
}



Mat Bottom_hat(Mat img, int mask_size) {
	Mat result(img.rows, img.cols, CV_8UC1, Scalar(0, 0, 0));
	Mat temp = Closing(img, mask_size);

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			result.at<uchar>(i, j) =  temp.at<uchar>(i, j) - img.at<uchar>(i, j);
		}
	return result;
}



Mat Dilation_gray(Mat img, int *mask, int mask_size) {
	int pad = mask_size / 2, check, max_val;
	Mat padded = zero_padding(img, pad);
	Mat result(img.rows, img.cols, CV_8UC1, Scalar(0, 0, 0));

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			max_val = 0;
			for (int m = 0; m < mask_size; m++)
				for (int n = 0; n < mask_size; n++) {
					check = mask[mask_size*m + n] + padded.at<uchar>(i + m, j + n);
					if (check > max_val) max_val = check;
				}
			if (max_val > 255) max_val = 255;
			result.at<uchar>(i, j) = max_val;
		}
	return result;
}



Mat Erosion_gray(Mat img, int *mask, int mask_size) {
	int pad = mask_size / 2, check, min_val;
	Mat padded = zero_padding(img, pad);
	Mat result(img.rows, img.cols, CV_8UC1, Scalar(0, 0, 0));

	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			min_val = 255;
			for (int m = 0; m < mask_size; m++)
				for (int n = 0; n < mask_size; n++) {
					check = padded.at<uchar>(i + m, j + n) - mask[mask_size*m + n];
					if (check < min_val) min_val = check;
				}
			if (min_val < 0) min_val = 0;
			result.at<uchar>(i, j) = min_val;
		}
	return result;
}



Mat Opening_gray(Mat img, int *mask, int mask_size) {
	Mat result = Erosion_gray(img, mask, mask_size);
	result = Dilation_gray(result, mask, mask_size);
	return result;
}



Mat Closing_gray(Mat img, int *mask, int mask_size) {
	Mat result = Dilation_gray(img, mask, mask_size);
	result = Erosion_gray(result, mask, mask_size);
	return result;
}