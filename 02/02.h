#pragma once

int *check_freq(cv::Mat img, int histSize);
int check_height(int* freq, int histSize);
void PlotPDF(int* freq, int h, int bin_w, double bin_h, int histSize);
void PlotCDF(int* freq, int h, int bin_w, double bin_h, int histSize);
void PlotPDF_threshold(int* freq, int h, int bin_w, double bin_h, int histSize, int threshold);
void PlotCDF_threshold(int* freq, int h, int bin_w, double bin_h, int histSize, int threshold);

// Histogram
void Sliding(cv::Mat img, cv::Mat& result, int offset);
void Stretching(cv::Mat img, cv::Mat& result, int max_val, int min_val);
void Shrink(cv::Mat img, cv::Mat& result, int max_val, int min_val, int s_max, int s_min);
void Equalize(cv::Mat img, cv::Mat& result, int total, int* freq, int histSize, int max_val, int min_val);
void Quantization(cv::Mat img, cv::Mat& result, int L, int histSize);

// Threshold
void Static_bin(cv::Mat img, cv::Mat& result, int T, int* check);
void P_tile_bin(cv::Mat img, cv::Mat& result, int P, int* freq, int total, int histSize, int* check);
void Iter_bin(cv::Mat img, cv::Mat& result, int threshold, int eps, int* freq, int histSize, int* check);
void Otsu(cv::Mat img, cv::Mat& result, int total, int* freq, int histSize, int max_val, int min_val, int* check);
void Valley(cv::Mat img, cv::Mat& result, int total, int* freq, int histSize, int max_val, int min_val, int* check);