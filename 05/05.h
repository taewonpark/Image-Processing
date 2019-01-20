#pragma once

cv::Mat Flood_fill(cv::Mat img);
void Ff_recursive(int **arr, int **temp, int m, int n, int label);
cv::Mat Two_pass(cv::Mat img);