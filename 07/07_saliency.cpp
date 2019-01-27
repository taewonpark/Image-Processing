#include <opencv2\opencv.hpp>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "..\\02.h"
#include "..\\03.h"
#include "..\\04.h"
#include "..\\07.h"

using namespace cv;
using namespace std;

const float PI = 3.14159265358979323846f;


//pyramid Node
class pNode
{
public:
	pNode* next;
	Mat data;
};

//pyramid list
class PyramidList
{
public:
	int length;
	pNode* head;
	PyramidList* next;

	PyramidList();
	~PyramidList();
	void add(Mat data);
};

PyramidList::PyramidList() {
	this->length = 0;
	this->head = NULL;
}

PyramidList::~PyramidList() {
	pNode* next = this->head;
	pNode* cur = NULL;
	while (next != NULL) {
		cur = next;
		next = next->next;
		delete cur;
	}
}


void PyramidList::add(Mat data) {
	pNode* node = new pNode();
	node->data = data;
	node->next = this->head;
	this->head = node;
	this->length++;
}

//type list
class TypeList
{
public:
	int length;
	PyramidList* head;

	TypeList();
	~TypeList();
	void add();
};

TypeList::TypeList() {
	this->length = 0;
	this->head = NULL;
}

TypeList::~TypeList() {
	PyramidList *next = this->head;
	PyramidList *cur = NULL;
	while (next != NULL) {
		cur = next;
		next = next->next;
		delete cur;
	}
}

void TypeList::add() {
	PyramidList* type = new PyramidList();
	type->next = this->head;
	this->head = type;
	this->length++;
}



Mat Saliency(Mat img) {
	//return Last_summation(Normalization(For_color(img)), Normalization(For_intensity(img)), Normalization(For_orientation(img)));
	//return Last_summation(For_color(img), For_intensity(img), For_orientation(img));

	//Strenching after summation
	double min, max;
	Mat temp = Last_summation(Normalization(For_color(img)), Normalization(For_intensity(img)), Normalization(For_orientation(img)));
	Mat result = temp.clone();
	minMaxLoc(temp, &min, &max);
	Stretching(temp, result, max, min);
	return result;
}



Mat For_color(Mat img) {
	//for RGBY
	TypeList* RGBY = new TypeList();
	for (int i = 0; i < 4; i++) RGBY->add();
	//[ R(Head), G, B, Y ]
	int ***color_arr = color(img, 4, img.rows, img.cols);
	//for pyramid ( sigma : [ 8(Head), 7, 6, 5, 4, 3, 2] )
	Mat temp_img;
	PyramidList* temp_Pyram;
	for (int i = 0; i < 4; i++) {
		temp_Pyram = RGBY->head;
		for (int j = 0; j < i; j++) temp_Pyram = temp_Pyram->next;
		temp_img = Octave(Gaussian(Arr2Mat(color_arr[i], img.rows, img.cols), 1.0)); // 1
		for (int j = 2; j <= 8; j++) {
			temp_img = Octave(temp_img);
			temp_Pyram->add(temp_img.clone());
		}
	}
	Deallocation_3D(color_arr, 4, img.rows);
	//CSD and ASC&N
	int H = img.rows, W = img.cols;
	for (int i = 0; i < 4; i++) {
		H = (H + 1) / 2; W = (W + 1) / 2;
	}
	Mat result(H, W, CV_8UC1, Scalar(0, 0, 0)), firstC, firstS, firstStemp, secondC, secondS, secondStemp;
	PyramidList *first, *second;
	pNode *firstNode, *secondNode;
	for (int c = 2; c <= 4; c++)
		for (int sigma = 3; sigma <= 4; sigma++) {
			for (int i = 0; i < 2; i++) {
				first = RGBY->head; second = RGBY->head;
				//find each PyramidList
				for (int j = 0; j < 2 * i; j++) first = first->next;
				for (int j = 0; j < 2 * i + 1; j++) second = second->next;
				//find each image
				firstNode = first->head; secondNode = second->head;
				if (c + sigma == 8) { firstStemp = firstNode->data; secondStemp = secondNode->data; }
				for (int j = 7; j >= 2; j--) {
					firstNode = firstNode->next; secondNode = secondNode->next;
					if(j == c) { firstC = firstNode->data; secondC = secondNode->data; }
					if (j == c + sigma) { firstStemp = firstNode->data; secondStemp = secondNode->data; }
				}
				//interpolation
				firstS = firstC.clone(); secondS = secondC.clone();
				resize(firstStemp, firstS, firstS.size(), 0, 0, INTER_CUBIC);
				resize(secondStemp, secondS, secondS.size(), 0, 0, INTER_CUBIC);
				//Center-surround
				if (i == 0)
					temp_img = Normalization(Img_absSubtract(Img_Subtract(firstC, secondC), Img_Subtract(secondS, firstS)));
				else
					Img_Plus(Normalization(Img_absSubtract(Img_Subtract(firstC, secondC), Img_Subtract(secondS, firstS))), temp_img);
			}
			for (int i = c; i < 4; i++) temp_img = Down_sampling(temp_img);
			Img_Plus(temp_img, result);
		}
	//deallocate every List
	delete RGBY;
	return result;
}



Mat For_intensity(Mat img) {
	//for intensity
	PyramidList* inten = new PyramidList();
	//intensity
	int **intensity_arr = intensity(img, img.rows, img.cols);
	//for pyramid ( sigma : [ 8(Head), 7, 6, 5, 4, 3, 2] )
	Mat temp_img;
	temp_img = Octave(Gaussian(Arr2Mat(intensity_arr, img.rows, img.cols), 1.0)); // 1
	for (int i = 2; i <= 8; i++) {
		temp_img = Octave(temp_img);
		inten->add(temp_img.clone());
	}
	Deallocation_2D(intensity_arr, img.rows);
	//CSD and ASC&N
	int H = img.rows, W = img.cols;
	for (int i = 0; i < 4; i++) {
		H = (H + 1) / 2; W = (W + 1) / 2;
	}
	Mat result(H, W, CV_8UC1, Scalar(0, 0, 0)), I_C, I_S, I_Stemp;
	pNode *I_Node;
	for (int c = 2; c <= 4; c++)
		for (int sigma = 3; sigma <= 4; sigma++) {
			//find each image
			I_Node = inten->head;
			if (c + sigma == 8) { I_Stemp = I_Node->data; }
			for (int j = 7; j >= 2; j--) {
				I_Node = I_Node->next;
				if (j == c) { I_C = I_Node->data; }
				if (j == c + sigma) { I_Stemp = I_Node->data; }
			}
			//interpolation
			I_S = I_C.clone();
			resize(I_Stemp, I_S, I_S.size(), 0, 0, INTER_CUBIC);
			//Center-surround
			temp_img = Normalization(Img_absSubtract(I_C, I_S));
			for (int i = c; i < 4; i++) temp_img = Down_sampling(temp_img);
			Img_Plus(temp_img, result);
		}
	//deallocate every List
	delete inten;
	return result;
}



Mat For_orientation(Mat img) {
	//for orientation
	TypeList* Orient = new TypeList();
	for (int i = 0; i < 4; i++) Orient->add();
	//[ 0(Head), 45, 90, 135 ]
	int ***orientation_arr = orientation(img, 4, img.rows, img.cols);
	//for pyramid ( sigma : [ 8(Head), 7, 6, 5, 4, 3, 2] )
	Mat temp_img;
	PyramidList* temp_Pyram;
	for (int i = 0; i < 4; i++) {
		temp_Pyram = Orient->head;
		for (int j = 0; j < i; j++) temp_Pyram = temp_Pyram->next;
		temp_img = Octave(Gaussian(Arr2Mat(orientation_arr[i], img.rows, img.cols), 1.0)); // 1
		for (int j = 2; j <= 8; j++) {
			temp_img = Octave(temp_img);
			temp_Pyram->add(temp_img.clone());
		}
	}
	Deallocation_3D(orientation_arr, 4, img.rows);
	//CSD and ASC&N
	int H = img.rows, W = img.cols;
	for (int i = 0; i < 4; i++) {
		H = (H + 1) / 2; W = (W + 1) / 2;
	}
	Mat result(H, W, CV_8UC1, Scalar(0, 0, 0)), O_C, O_S, O_Stemp, orient_img;
	PyramidList *O_Pyram;
	pNode *O_Node;
	for (int i = 0; i < 4; i++) {
		O_Pyram = Orient->head;
		//find each PyramidList
		for (int j = 0; j < i; j++) O_Pyram = O_Pyram->next;
		//c - s
		orient_img = Mat(H, W, CV_8UC1, Scalar(0, 0, 0));
		for (int c = 2; c <= 4; c++)
			for (int sigma = 3; sigma <= 4; sigma++) {
				//find each image
				O_Node = O_Pyram->head;
				if (c + sigma == 8) { O_Stemp = O_Node->data; }
				for (int j = 7; j >= 2; j--) {
					O_Node = O_Node->next;
					if (j == c) { O_C = O_Node->data; }
					if (j == c + sigma) { O_Stemp = O_Node->data; }
				}
				/*imshow("a", O_C); waitKey();
				imshow("b", O_Stemp); waitKey();*/
				//interpolation
				O_S = O_C.clone();
				resize(O_Stemp, O_S, O_S.size(), 0, 0, INTER_CUBIC);
				//Center-surround
				temp_img = Normalization(Img_absSubtract(O_C, O_S));
				//imshow("a", temp_img); waitKey();
				for (int i = c; i < 4; i++) temp_img = Down_sampling(temp_img);
				Img_Plus(temp_img, orient_img);
			}
		Img_Plus(Normalization(orient_img), result);
	}
	//deallocate every List
	delete Orient;
	return result;
}



int*** color(Mat img, int first, int second, int third) {
	//RGBY
	//memory allocation for each color value
	//[ R, G, B, Y ]
	// R, G, B, Y : 2D array
	int r, g, b, rgb[3] = { 0,0,0 }, ***arr = Allocation_3D(first, second, third);
	//normalization
	int temp, i_max, threshold, check_max[3] = { 0,0,0 }, check_min[3] = { 255, 255, 255 };
	color_threshold(img, &threshold, &i_max);
	for (int i = 0; i<img.rows; i++)
		for (int j = 0; j < img.cols; j++) {
			temp = 0;
			for (int ch = 0; ch < 3; ch++ ) {
				rgb[ch] = img.at<Vec3b>(i, j)[2 - ch];
				temp += rgb[ch];
			}
			if (temp / 3 >= threshold)
				for (int ch = 0; ch < 3; ch++) {
					if (check_max[ch] < rgb[ch]) check_max[ch] = rgb[ch];
					if (check_min[ch] > rgb[ch]) check_min[ch] = rgb[ch];
					arr[ch][i][j] = rgb[ch];
				}
		}
	for (int ch = 0; ch < 3; ch++)
		for (int i = 0; i < img.rows; i++)
			for (int j = 0; j < img.cols; j++) {
				arr[ch][i][j] = arr[ch][i][j] * i_max / (check_max[ch] - check_min[ch]);
			}
	//calculate RGBY
	for (int i = 0; i<img.rows; i++)
		for (int j = 0; j < img.cols; j++) {
			r = arr[0][i][j];
			g = arr[1][i][j];
			b = arr[2][i][j];
			if ((r + g + b) / 3 >= threshold) {
				arr[0][i][j] = r - (g + b) / 2;// R
				arr[1][i][j] = g - (b + r) / 2;// G
				arr[2][i][j] = b - (r + g) / 2;// B
				arr[3][i][j] = (r + g) / 2 - abs(r - g) / 2 - b;// Y
				for (int k = 0; k < 4; k++) {
					if (arr[k][i][j] < 0) arr[k][i][j] = 0;
					else if (arr[k][i][j] > 255) arr[k][i][j] = 255;
				}
			}
		}
	return arr;
}



int** intensity(Mat img, int first, int second) {
	//memory allocation for intensity
	// intensity : 2D array
	int **arr = Allocation_2D(first, second);
	//calculate intensity
	int temp;
	for (int i = 0; i < first; i++)
		for (int j = 0; j < second; j++) {
			temp = 0;
			for (int k = 0; k < 3; k++) temp += img.at<Vec3b>(i, j)[k];
			temp /= 3;
			if (temp > 255) temp = 255;
			arr[i][j] = temp;
		}
	return arr;
}



int*** orientation(Mat img, int first, int second, int third) {
	//memory allocation for each color value
	//[ 0, 45, 90, 135 ]
	// 0, 45, 90, 135 : 2D array
	int ***arr = Allocation_3D(first, second, third);
	//calculate orientation
	Mat gray, temp, kernel, dest, mag;
	cvtColor(img, temp, COLOR_RGB2GRAY);
	temp.convertTo(gray, CV_32F, 10.0 / 255, 0);
	int k_size = 21;
	double theta[4] = { 0.0, 45.0, 90.0, 135.0 };
	double sig = 5, lm = 1, gm = 1, ps = 90;
	for (int i = 0; i < 4; i++) {
		kernel = gaborKernel(k_size, sig, theta[i], lm, ps);
		filter2D(gray, dest, CV_32F, kernel);
		pow(dest, 2.0, mag);
		mag.convertTo(temp, CV_8UC1);
		for (int m = 0; m < img.rows; m++)
			for (int n = 0; n < img.cols; n++)
				arr[i][m][n] = (int)mag.at<float>(m, n);
	}
	return arr;
}



Mat Arr2Mat(int** arr, int H, int W) {
	Mat result(H, W, CV_8UC1, Scalar(0, 0, 0));
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			result.at<uchar>(i, j) = arr[i][j];
	return result;
}



void Mat2Arr(Mat img, int** arr) {
	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols; j++)
			arr[i][j] = img.at<uchar>(i, j);
}



Mat Gaussian(Mat img, float sigma) {
	int padding_size = 1, k_size = 2 * padding_size + 1;
	float *mask = new float[k_size*k_size];

	Get_gaussian(mask, k_size, sigma);
	return Convolution_F(zero_padding(img, padding_size), mask, img.rows, img.cols);
}



Mat Octave(Mat img) {
	int padding_size = 1, k_size = 2 * padding_size + 1;
	float *mask = new float[k_size*k_size];
	
	Get_gaussian(mask, k_size, 1.0);
	//make scale 2
	Mat temp1 = Convolution_F(zero_padding(img, padding_size), mask, img.rows, img.cols);
	Mat temp2 = Convolution_F(zero_padding(temp1, padding_size), mask, img.rows, img.cols);
	Mat temp3 = Convolution_F(zero_padding(temp2, padding_size), mask, img.rows, img.cols);
	//down sampling
	return Down_sampling(temp3);
}



Mat Down_sampling(Mat img) {
	Mat result((img.rows + 1) / 2, (img.cols + 1) / 2, CV_8UC1, Scalar(0, 0, 0));
	for (int i = 0; i < result.rows; i ++)
		for (int j = 0; j < result.cols; j ++) {
			result.at<uchar>(i, j) = img.at<uchar>(2 * i, 2 * j);
		}
	return result;
}



int*** Allocation_3D(int first, int second, int third) {
	int ***arr = new int**[first];
	for (int i = 0; i < first; i++) {
		arr[i] = new int*[second];
		for (int j = 0; j < second; j++) {
			arr[i][j] = new int[third];
			memset(arr[i][j], 0, sizeof(int)*third);
		}
	}
	return arr;
}



int** Allocation_2D(int first, int second) {
	int **arr = new int*[first];
	for (int i = 0; i < first; i++) {
		arr[i] = new int[second];
		memset(arr[i], 0, sizeof(int)*second);
	}
	return arr;
}



void Deallocation_3D(int ***arr, int first, int second) {
	for (int i = 0; i < first; i++)
		for (int j = 0; j < second; j++)
			delete[] arr[i][j];
	for (int i = 0; i < first; i++)
		delete[] arr[i];
	delete[] arr;
}



void Deallocation_2D(int **arr, int first) {
	for (int i = 0; i < first; i++)
		delete[] arr[i];
	delete[] arr;
}



Mat Img_absSubtract(Mat img1, Mat img2) {
	int temp;
	Mat result(img1.rows, img1.cols, CV_8UC1, Scalar(0, 0, 0));
	for (int i = 0; i < img1.rows; i++)
		for (int j = 0; j < img1.cols; j++) {
			temp = abs(img1.at<uchar>(i, j) - img2.at<uchar>(i, j));
			if (temp > 255) temp = 255;
			result.at<uchar>(i, j) = temp;
		}
	return result;
}



Mat Img_Subtract(Mat img1, Mat img2) {
	int temp;
	Mat result(img1.rows, img1.cols, CV_8UC1, Scalar(0, 0, 0));
	for (int i = 0; i < img1.rows; i++)
		for (int j = 0; j < img1.cols; j++) {
			temp = img1.at<uchar>(i, j) - img2.at<uchar>(i, j);
			if (temp < 0) temp = 0;
			else if (temp > 255) temp = 255;
			result.at<uchar>(i, j) = temp;
		}
	return result;
}



Mat Normalization(Mat img) {
	Mat temp, result(img.rows, img.cols, CV_8UC1, Scalar(0, 0, 0)), preliminary=img.clone();
	int center, count = 0, M = 4;
	double max, min, avg, temp_max, temp_min, summation=0;
	minMaxLoc(img, &min, &max);
	if (max != 0) {
		for (int i = 0; i < img.rows - 2; i++)
			for (int j = 0; j < img.cols - 2; j++) {
				temp = preliminary(Rect(j, i, 3, 3));
				center = temp.at<uchar>(1, 1);
				temp.at<uchar>(1, 1) = 0;
				minMaxLoc(temp, &temp_min, &temp_max);
				if (center > temp_max)
					if (center < max) {
						summation += center;
						count++;
					}
			}
		avg = (double)summation*M / (max - min) / count;
		for (int i = 0; i < img.rows; i++)
			for (int j = 0; j < img.cols; j++) {
				result.at<uchar>(i, j) = img.at<uchar>(i, j)*M*(M - avg)*(M - avg) / (max - min);
			}
		return result;
	}else if (max == 0) return img;
}



void Img_Plus(Mat img, Mat& dest) {
	int temp;
	for (int i = 0; i < img.rows; i++)
		for (int j = 0; j < img.cols; j++) {
			uchar* pixel = dest.ptr<uchar>(i) + j;
			temp = (int)*pixel + img.at<uchar>(i, j);
			if (temp < 0) temp = 0;
			else if (temp > 255) temp = 255;
			*pixel = temp;
		}
}



void color_threshold(Mat img, int* threshold, int* i_max) {
	double max, min;
	Mat gray;
	cvtColor(img, gray, COLOR_RGB2GRAY);
	minMaxLoc(gray, &min, &max);
	*threshold = (int)max / 10;
	*i_max = max;
}



Mat Last_summation(Mat img1, Mat img2, Mat img3) {
	Mat result(img1.rows, img1.cols, CV_8UC1, Scalar(0, 0, 0));
	int temp;
	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++) {
			temp = ((int)img1.at<uchar>(i, j) + (int)img2.at<uchar>(i, j) + (int)img3.at<uchar>(i, j)) / 3;
			if (temp < 0) temp = 0;
			else if (temp > 255) temp = 255;
			result.at<uchar>(i, j) = temp;
		}
	return result;
}



Mat gaborKernel(int ks, double sig, double th, double lm, double ps){
	int hks = (ks - 1) / 2;
	double theta = th*CV_PI / 180;
	double psi = ps*CV_PI / 180;
	double del = 2.0 / (ks - 1);
	double lmbd = lm;
	double sigma = sig / ks; 
	double x_theta;
	double y_theta;
	cv::Mat kernel(ks, ks, CV_32F);
	for (int y = -hks; y <= hks; y++)
	{
		for (int x = -hks; x <= hks; x++)
		{
			x_theta = x*del*cos(theta) + y*del*sin(theta);
			y_theta = -x*del*sin(theta) + y*del*cos(theta);
			kernel.at<float>(hks + y, hks + x) = (float)exp(-0.5*(pow(x_theta, 2) + pow(y_theta, 2)) / pow(sigma, 2))* cos(2 * CV_PI*x_theta / lmbd + psi);
		}
	}
	return kernel;
}