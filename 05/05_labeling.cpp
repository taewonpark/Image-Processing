#include <opencv2\opencv.hpp>
#include <iostream>
#include <cmath>
#include <string.h>
#include "05.h"

using namespace cv;
using namespace std;


Mat Flood_fill(Mat img) {
	Mat result(img.rows, img.cols, CV_8UC1, Scalar(0, 0, 0));
	int **arr = new int*[img.rows], **temp = new int*[img.rows];
	int interval, label = 1;

	for (int i = 0; i < img.rows; i++) {
		arr[i] = new int[img.cols];
		temp[i] = new int[img.cols];
		memset(temp[i], 0, sizeof(int)*img.cols);
	}
	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			arr[i][j] = img.at<uchar>(i, j);
		}
	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++)
			if (arr[i][j] != 0 && temp[i][j] == 0){
				Ff_recursive(arr, temp, i, j, label);
				label++;
			}
	interval = 255 / label;
	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			result.at<uchar>(i, j) = temp[i][j] * interval;
		}
	for (int i = 0; i < img.rows; i++) {
		delete[] arr[i];
		delete[] temp[i];
	}
	return result;
}



void Ff_recursive(int **arr, int **temp, int m, int n, int label) {
	if (arr[m][n] != 0 && temp[m][n] == 0){
		temp[m][n] = label;
		/*
		//4-neighbourhood
		Ff_recursive(arr, temp, m, n + 1, label);
		Ff_recursive(arr, temp, m - 1, n, label);
		Ff_recursive(arr, temp, m, n - 1, label);
		Ff_recursive(arr, temp, m + 1, n, label);
		*/
		
		//8-neighbourhood
		Ff_recursive(arr, temp, m, n + 1, label);
		Ff_recursive(arr, temp, m - 1, n + 1, label);
		Ff_recursive(arr, temp, m - 1, n, label);
		Ff_recursive(arr, temp, m - 1, n - 1, label);
		Ff_recursive(arr, temp, m, n - 1, label);
		Ff_recursive(arr, temp, m + 1, n - 1, label);
		Ff_recursive(arr, temp, m + 1, n, label);
		Ff_recursive(arr, temp, m + 1, n + 1, label);
	}
}



Mat Two_pass(Mat img) {
	Mat result(img.rows, img.cols, CV_8UC1, Scalar(0, 0, 0));
	int **temp = new int*[img.rows];
	int a, b, interval, label_size = 500, label = 0;
	int* table = new int[label_size];
	int* equivalent_table = new int[label_size];
	//initialization
	memset(table, 0, sizeof(int)*label_size);
	memset(equivalent_table, 0, sizeof(int)*label_size);
	for (int i = 0; i < img.rows; i++) {
		temp[i] = new int[img.cols];
		memset(temp[i], 0, sizeof(int)*img.cols);
	}
	//1-pass
	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++)
			if (img.at<uchar>(i, j) != 0) {
				a = table[temp[i - 1][j]];//above pixel's label
				b = table[temp[i][j - 1]];//below pixel's label
				if ((a == 0) && (b == 0)) {
					//allocate new label value
					temp[i][j] = ++label;
					table[label] = label;
				}else if ((a != 0) && (b != 0)) {
					//consider priority
					if (a > b) {
						temp[i][j] = b;
						table[temp[i - 1][j]] = b;
					}else {
						temp[i][j] = a;
						table[temp[i][j - 1]] = a;
					}
				}else if (a != 0) temp[i][j] = a;
				else if (b != 0) temp[i][j] = b;
			}
	//take care of unhandled table part
	for (int i = 0; i < label_size; i++) table[i] = table[table[i]];
	//make equivalent table
	label = 0;
	for (int i = 0; i < label_size; i++)
		if (table[i] != 0)
			if (equivalent_table[table[i]] == 0)
				equivalent_table[table[i]] = ++label;
	interval = 255 / label;
	//2-pass
	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			result.at<uchar>(i, j) = equivalent_table[table[temp[i][j]]] * interval;
		}
	//remove dynamic allocation
	for (int i = 0; i < img.rows; i++) {
		delete[] temp[i];
	}
	return result;
}