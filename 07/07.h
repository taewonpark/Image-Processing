#pragma once

cv::Mat Saliency(cv::Mat img);
cv::Mat For_color(cv::Mat img);
cv::Mat For_intensity(cv::Mat img);
cv::Mat For_orientation(cv::Mat img);
int*** color(cv::Mat img, int first, int second, int third);
int** intensity(cv::Mat img, int first, int second);
int*** orientation(cv::Mat img, int first, int second, int third);
cv::Mat Arr2Mat(int** arr, int H, int W);
void Mat2Arr(cv::Mat img, int** arr);
cv::Mat Gaussian(cv::Mat img, float sigma);
cv::Mat Octave(cv::Mat img);
cv::Mat Down_sampling(cv::Mat img);
int*** Allocation_3D(int first, int second, int third);
int** Allocation_2D(int first, int second);
void Deallocation_3D(int ***arr, int first, int second);
void Deallocation_2D(int **arr, int first);
cv::Mat Img_absSubtract(cv::Mat img1, cv::Mat img2);
cv::Mat Img_Subtract(cv:: Mat img1, cv::Mat img2);
cv::Mat Normalization(cv::Mat img);
void Img_Plus(cv::Mat img, cv::Mat& dest);
void color_threshold(cv::Mat img, int* threshold, int* i_max);
cv::Mat Last_summation(cv::Mat img1, cv::Mat img2, cv::Mat img3);
cv::Mat gaborKernel(int ks, double sig, double th, double lm, double ps);