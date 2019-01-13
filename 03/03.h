#pragma once

cv::Mat zero_padding(cv::Mat img, int pad);
cv::Mat edge_padding(cv::Mat img, int pad);
cv::Mat reflect_padding(cv::Mat img, int pad);

void Seletionsort(int* arr, int size);
void Quicksort(int* arr, int left, int right);

cv::Mat Average_Filter(cv::Mat img, int height, int width);
cv::Mat Median_Filter(cv::Mat img, int height, int width);
cv::Mat Maximum_Filter(cv::Mat img, int height, int width);
cv::Mat Minimum_Filter(cv::Mat img, int height, int width);
cv::Mat Midpoint_Filter(cv::Mat img, int height, int width);
cv::Mat Alpha_Filter(cv::Mat img, int height, int width);

cv::Mat Convolution(cv::Mat img, int* mask, int height, int width);
cv::Mat Convolution_F(cv::Mat img, float* mask, int height, int width);
cv::Mat RSS(cv::Mat x_img, cv::Mat y_img);

void Get_sobel(float* x_mask, float* y_mask, int size);
void Get_prewitt(float* x_mask, float* y_mask, int size);
void Get_roberts(float* x_mask, float* y_mask, int size);
void Get_laplacian(float *without, float* with, int size);

void Get_gaussian(float *mask, int size, float stddev);
//void Get_LoG(float *mask, int size, float stddev);
//void Get_DoG(float *mask, int size, float stddev1, float stddev2);
//void Get_Highboost(int* x_mask, int* y_mask, int size, int c);

cv::Mat adding(cv::Mat x_img, cv::Mat y_img, float alpha);