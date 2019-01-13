#include <opencv2\opencv.hpp>
#include <iostream>
#include <cmath>
#include "02.h"

using namespace cv;
using namespace std;


int main() {
	int min_val, max_val;
	int histSize = 256, height = 400, width = 1200, h, threshold;
	int bin_w = cvRound((double)width / histSize);

	string imageName("C:\\Users\\samsung\\Desktop\\BEP\\BEP2\\girl.jpg");
	Mat img = imread(imageName.c_str(), IMREAD_COLOR), gray, histogram;

	// Convert GBR to gray
	cvtColor(img, gray, COLOR_RGB2GRAY);

	// find each grayscale frequency
	int *freq = check_freq(gray, histSize);

	// find total # of pixel
	double total = gray.cols*gray.rows;

	// find height to plot
	h = check_height(freq, histSize);

	// declare temporary variable for plot
	Mat temp_img = gray.clone();
	int *temp_freq;

	// find min and max grayscale that pixel have
	for (int i = 0; i < histSize; i++)
		if (freq[i] != 0) { min_val = i; break; }

	for (int i = histSize-1; i >= 0; i--) 
		if (freq[i] != 0) { max_val = i; break; }
	
	// plot original image's pdf and cdf
	imshow("gray", gray);
	waitKey();
	PlotPDF(freq, height, bin_w, (double)height/h, histSize);
	waitKey();
	PlotCDF(freq, height, bin_w, (double)height / total, histSize);
	waitKey();

	Iter_bin(gray, temp_img, 150, 2, freq, histSize, &threshold);
	PlotPDF_threshold(freq, height, bin_w, (double)height / h, histSize, threshold);
	waitKey();

	/*
	Sliding(gray, temp_img, 100);
	Stretching(gray, temp_img, max_val, min_val);
	Shrink(gray, temp_img, max_val, min_val, 20, 100);
	Equalize(gray, temp_img, (int)total, freq, histSize, max_val, min_val);
	Quantization(gray, temp_img, 16, histSize);

	Static_bin(gray, temp_img, 150, &threshold);
	P_tile_bin(gray, temp_img, 4, freq, total, histSize, &threshold);
	Iter_bin(gray, temp_img, 150, 2, freq, histSize, &threshold);
	Otsu(gray, temp_img, total, freq, histSize, max_val, min_val, &threshold);
	Valley(gray, temp_img, total, freq, histSize, max_val, min_val, &threshold);

	imshow("S_b", temp_img);
	waitKey();
	temp_freq = check_freq(temp_img, histSize);
	*/
	return 0;
}