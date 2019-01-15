#pragma once

cv::Mat Dilation(cv::Mat img, int mask_size);
cv::Mat Erosion(cv::Mat img, int mask_size);
cv::Mat Opening(cv::Mat img, int mask_size);
cv::Mat Closing(cv::Mat img, int mask_size);

cv::Mat Complement(cv::Mat img);
cv::Mat Image_and(cv::Mat img1, cv::Mat img2);
cv::Mat Erosion_mask(cv::Mat img, int *mask, int mask_size);
cv::Mat Hit_or_miss(cv::Mat img, int mask_size, int *C, int *D);

cv::Mat Top_hat(cv::Mat img, int mask_size);
cv::Mat Bottom_hat(cv::Mat img, int mask_size);

cv::Mat Dilation_gray(cv::Mat img, int *mask, int mask_size);
cv::Mat Erosion_gray(cv::Mat img, int *mask, int mask_size);
cv::Mat Opening_gray(cv::Mat img, int *mask, int mask_size);
cv::Mat Closing_gray(cv::Mat img, int *mask, int mask_size);

cv::Mat Forward_mapping(cv::Mat img, int k);
cv::Mat Nearest_interpolation(cv::Mat img, int k);
cv::Mat Bilinear_interpolation(cv::Mat img, int k);
double Bicubic_weight(double x, double a);
cv::Mat Bicubic_interpolation(cv::Mat img, int k, double a);
cv::Mat Bicubic_interpolation_x_next_y(cv::Mat img, int k, double a);
double Bspline_weight(double x);
cv::Mat Bispline_interpolation(cv::Mat img, int k);
cv::Mat Bispline_interpolation_x_next_y(cv::Mat img, int k);

cv::Mat Rotation(cv::Mat img, double theta);
cv::Mat Rotation_with_cut(cv::Mat img, double theta);
cv::Mat Rotation_with_interpolation(cv::Mat img, double theta);