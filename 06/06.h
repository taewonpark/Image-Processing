#pragma once

cv::Mat DFT(cv:: Mat img, double **realImg, double **imagImg);
cv::Mat IDFT(double **realImg, double **imagImg, int H, int W);
cv::Mat Spectrum(double **realImg, double **imagImg, int H, int W);