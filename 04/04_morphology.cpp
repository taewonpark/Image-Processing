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



Mat Forward_mapping(Mat img, int k) {
	Mat result(k*img.rows, k*img.cols, CV_8UC1, Scalar(0, 0, 0));

	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols; j++)
			result.at<uchar>(i*k, j*k) = img.at<uchar>(i, j);
	return result;
}



Mat Nearest_interpolation(Mat img, int k) {
	Mat result(k*img.rows, k*img.cols, CV_8UC1, Scalar(0, 0, 0));
	
	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols; j++)
			for (int m = 0; m < k; m++)
				for (int n = 0; n < k; n++)
					result.at<uchar>(i*k + m, j*k + n) = img.at<uchar>(i, j);
	return result;
}



Mat Bilinear_interpolation(Mat img, int k) {
	Mat result(k*img.rows, k*img.cols, CV_8UC1, Scalar(0, 0, 0));
	double slope, interval = (double)1 / (k + 1);
	int val1, val2, temp_k = k + 1;
	//calculate x-axis pixel value
	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols - 1; j++) {
			val1 = img.at<uchar>(i, j);
			val2 = img.at<uchar>(i, j + 1);
			slope = (double)(val2 - val1) * interval;
			for (int m = 0; m < k + 2; m++)
				result.at<uchar>(i*temp_k, j*temp_k + m) = (int)(slope * m + val1);
		}
	//calculate y-axis pixel value
	for (int i = 0; i < img.rows - 1; i++)
		for (int j = 0; j < result.cols; j++) {
			val1 = result.at<uchar>(i*temp_k, j);
			val2 = result.at<uchar>((i + 1)*temp_k, j);
			slope = (double)((val2 - val1) * interval);
			for (int m = 1; m < temp_k; m++)
				result.at<uchar>(i*temp_k + m, j) = (int)(slope * m + val1);
		}
	return result;
}



double Bicubic_weight(double x, double a) {
	double temp = abs(x);

	if (temp < 1)
		return (a + 2.0)*temp*temp*temp - (a + 3.0)*temp*temp + 1.0;
	else if (temp >= 1 && temp < 2)
		return a*temp*temp*temp - 5.0 * a*temp*temp + 8.0 * a*temp - 4.0 * a;
	else
		return 0;
}



//consider 16 element at same moment
Mat Bicubic_interpolation(Mat img, int k, double a) {
	Mat result(k*img.rows, k*img.cols, CV_8UC1, Scalar(0, 0, 0));
	Mat padded = zero_padding(img, 1);
	double x, y, interval = (double)1 / (k + 1);
	double temp, val, distance;
	int temp_k = k + 1;

	for (int i = 0; i < padded.rows - 3; i++)
		for (int j = 0; j < padded.cols - 3; j++)
			for (int m = 0; m < k + 2; m++)
				for (int n = 0; n < k + 2; n++) {
					temp = 0;
					for (int o = 0; o < 4; o++)
						for (int p = 0; p < 4; p++) {
							x = (o - 1.5) - (m*interval - 0.5);
							y = (p - 1.5) - (n*interval - 0.5);
							distance = sqrt(x*x + y*y);
							val = padded.at<uchar>(i + o, j + p);
							temp += Bicubic_weight(distance, a) * val;
						}
					if (temp > 255) temp = 255;
					else if (temp < 0) temp = 0;
					result.at<uchar>(i*temp_k + m, j*temp_k + n) = temp;
				}
	return result;
}



//compute x-axis interpolation and then consider y-axis
Mat Bicubic_interpolation_x_next_y(Mat img, int k, double a) {
	Mat result(k*img.rows, k*img.cols, CV_8UC1, Scalar(0, 0, 0));
	double *left_two = new double[k + 2], *left_one = new double[k + 2];
	double *right_two = new double[k + 2], *right_one = new double[k + 2];
	double x, interval = (double)1 / (k + 1);
	double val1, val2, val3, val4;
	int temp, temp_k = k + 1;
	//calculate weight of each position
	for (int i = 0; i < k + 2; i++) {
		x = interval*i;
		left_two[i] = Bicubic_weight(x + 1.0, a);
		left_one[i] = Bicubic_weight(x, a);
		right_one[i] = Bicubic_weight(1.0 - x, a);
		right_two[i] = Bicubic_weight(2.0 - x, a);
	}
	//calculate x-axis pixel value
	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols - 1; j++) {
			(j == 0) ? val1 = 0 : val1 = img.at<uchar>(i, j - 1);
			val2 = img.at<uchar>(i, j);
			val3 = img.at<uchar>(i, j + 1);
			(j == img.cols - 2) ? val4 = 0 : val4 = img.at<uchar>(i, j + 2);
			for (int m = 0; m < k + 2; m++) {
				temp = (int)(val1*left_two[m] + val2*left_one[m] + val3*right_one[m] + val4*right_two[m]);
				if (temp > 255) temp = 255;
				else if (temp < 0) temp = 0;
				result.at<uchar>(i*temp_k, j*temp_k + m) = temp;
			}
		}
	//calculate y-axis pixel value
	for (int i = 0; i < img.rows - 1; i++)
		for (int j = 0; j < result.cols; j++) {
			(i == 0) ? val1 = 0 : val1 = result.at<uchar>((i - 1)*temp_k, j);
			val2 = result.at<uchar>(i*temp_k, j);
			val3 = result.at<uchar>((i + 1)*temp_k, j);
			(i == img.rows - 2) ? val4 = 0 : val4 = result.at<uchar>((i + 2)*temp_k, j);
			for (int m = 1; m < temp_k; m++) {
				temp = (int)(val1*left_two[m] + val2*left_one[m] + val3*right_one[m] + val4*right_two[m]);
				if (temp > 255) temp = 255;
				else if (temp < 0) temp = 0;
				result.at<uchar>(i*temp_k + m, j) = temp;
			}
		}
	return result;
}



double Bspline_weight(double x) {
	double temp = abs(x);

	if (temp < 1)
		return temp*temp*temp / 2.0 - temp*temp + 2.0 / 3.0;
	else if (temp >= 1 && temp < 2)
		return -temp*temp*temp / 6.0 + temp*temp - 2.0* temp + 4.0 / 3.0;
	else
		return 0;
}


//consider 16 element at same moment
Mat Bispline_interpolation(Mat img, int k) {
	Mat result(k*img.rows, k*img.cols, CV_8UC1, Scalar(0, 0, 0));
	Mat padded = zero_padding(img, 1);
	double x, y, interval = (double)1 / (k + 1);
	double temp, val, distance;
	int temp_k = k + 1;

	for (int i = 0; i < padded.rows - 3; i++)
		for (int j = 0; j < padded.cols - 3; j++)
			for (int m = 0; m < k + 2; m++)
				for (int n = 0; n < k + 2; n++) {
					temp = 0;
					for (int o = 0; o < 4; o++)
						for (int p = 0; p < 4; p++) {
							x = (o - 1.5) - (m*interval - 0.5);
							y = (p - 1.5) - (n*interval - 0.5);
							distance = sqrt(x*x + y*y);
							val = padded.at<uchar>(i + o, j + p);
							temp += Bspline_weight(distance) * val;
						}
					if (temp > 255) temp = 255;
					else if (temp < 0) temp = 0;
					result.at<uchar>(i*temp_k + m, j*temp_k + n) = temp;
				}
	return result;
}




Mat Bispline_interpolation_x_next_y(Mat img, int k) {
	Mat result(k*img.rows, k*img.cols, CV_8UC1, Scalar(0, 0, 0));
	double *left_two = new double[k + 2], *left_one = new double[k + 2];
	double *right_two = new double[k + 2], *right_one = new double[k + 2];
	double x, interval = (double)1 / (k + 1);
	double val1, val2, val3, val4;
	int temp, temp_k = k + 1;
	//calculate weight of each position
	for (int i = 0; i < k + 2; i++) {
		x = interval*i;
		left_two[i] = Bspline_weight(x + 1.0);
		left_one[i] = Bspline_weight(x);
		right_one[i] = Bspline_weight(1.0 - x);
		right_two[i] = Bspline_weight(2.0 - x);
	}
	//calculate x-axis pixel value
	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols - 1; j++) {
			(j == 0) ? val1 = 0 : val1 = img.at<uchar>(i, j - 1);
			val2 = img.at<uchar>(i, j);
			val3 = img.at<uchar>(i, j + 1);
			(j == img.cols - 2) ? val4 = 0 : val4 = img.at<uchar>(i, j + 2);
			for (int m = 0; m < k + 2; m++) {
				temp = (int)(val1*left_two[m] + val2*left_one[m] + val3*right_one[m] + val4*right_two[m]);
				if (temp > 255) temp = 255;
				else if (temp < 0) temp = 0;
				result.at<uchar>(i*temp_k, j*temp_k + m) = temp;
			}
		}
	//calculate y-axis pixel value
	for (int i = 0; i < img.rows - 1; i++)
		for (int j = 0; j < result.cols; j++) {
			(i == 0) ? val1 = 0 : val1 = result.at<uchar>((i - 1)*temp_k, j);
			val2 = result.at<uchar>(i*temp_k, j);
			val3 = result.at<uchar>((i + 1)*temp_k, j);
			(i == img.rows - 2) ? val4 = 0 : val4 = result.at<uchar>((i + 2)*temp_k, j);
			for (int m = 0; m < temp_k+1; m++) {
				temp = (int)(val1*left_two[m] + val2*left_one[m] + val3*right_one[m] + val4*right_two[m]);
				if (temp > 255) temp = 255;
				else if (temp < 0) temp = 0;
				result.at<uchar>(i*temp_k + m, j) = temp;
			}
		}
	return result;
}



Mat Rotation(Mat img, double theta) {
	double cos_ = cos(theta), sin_ = sin(theta);
	int img_center_row = (img.rows - 1) / 2, img_center_col = (img.cols - 1) / 2;
	int total_cols = img.rows*sin_ + img.cols*cos_;
	int total_rows = img.rows*cos_ + img.cols*sin_;
	int total_center_row = (total_rows - 1) / 2, total_center_col = (total_cols - 1) / 2;
	int m, n;

	Mat result(total_rows, total_cols, CV_8UC1, Scalar(0, 0, 0));

	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols; j++) {
			m = (i - img_center_row)*cos_ + (j - img_center_col)*sin_;
			n = (i - img_center_row)*(-sin_) + (j - img_center_col)*cos_;
			result.at<uchar>(m + total_center_row, n + total_center_col) = img.at<uchar>(i, j);
		}
	return result;
}



Mat Rotation_with_cut(Mat img, double theta) {
	double cos_ = cos(theta), sin_ = sin(theta);
	int img_center_row = (img.rows - 1) / 2, img_center_col = (img.cols - 1) / 2;
	int m, n;

	Mat result(img.rows, img.cols, CV_8UC1, Scalar(0, 0, 0));

	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols; j++) {
			m = (i - img_center_row)*cos_ + (j - img_center_col)*sin_;
			n = (i - img_center_row)*(-sin_) + (j - img_center_col)*cos_;
			m += img_center_row; n += img_center_col;
			if (m >= 0 && m < img.rows && n >= 0 && n < img.cols)
				result.at<uchar>(m , n) = img.at<uchar>(i, j);
		}
	return result;
}



Mat Rotation_with_interpolation(Mat img, double theta) {
	double cos_ = cos(theta), sin_ = sin(theta);
	int img_center_row = (img.rows - 1) / 2, img_center_col = (img.cols - 1) / 2;
	int total_cols = img.rows*sin_ + img.cols*cos_;
	int total_rows = img.rows*cos_ + img.cols*sin_;
	int total_center_row = (total_rows - 1) / 2, total_center_col = (total_cols - 1) / 2;
	int m, n, count, summation, temp;

	Mat result(total_rows, total_cols, CV_8UC1, Scalar(0, 0, 0));
	Mat check(total_rows, total_cols, CV_8UC1, Scalar(0, 0, 0));

	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols; j++) {
			m = (i - img_center_row)*cos_ + (j - img_center_col)*sin_;
			n = (i - img_center_row)*(-sin_) + (j - img_center_col)*cos_;
			result.at<uchar>(m + total_center_row, n + total_center_col) = img.at<uchar>(i, j);
			check.at<uchar>(m + total_center_row, n + total_center_col) = 1;
		}

	Mat padded_result = zero_padding(result, 1);
	Mat padded_check = zero_padding(check, 1);
	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			if (padded_check.at<uchar>(i + 1, j + 1) == 0) {
				count = 0, summation = 0;
				for (int m = 0; m < 3; m++)
					for (int n = 0; n < 3; n++) {
						temp = padded_check.at<uchar>(i + m, j + n);
						count += temp;
						summation += temp*padded_result.at<uchar>(i + m, j + n);
					}
				if (count > 0) {
					summation /= count;
					result.at<uchar>(i, j) = summation;
				}
			}
		}
	return result;
}