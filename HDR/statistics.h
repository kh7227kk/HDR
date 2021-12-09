#pragma once
#include <opencv2/opencv.hpp>

//#include "mydefine.h"

void Histogram(cv::Mat dst, int x, int y, int col_size, int row_size, std::vector<int>& val);

void mean_Histogram(std::vector<int> val, double& mean);

void SD_Histogram(std::vector<int> val, double mean, double& SD);

//void printHistogram(Filedata img);

void Regression(std::string path, const int p, std::vector<double>& coeff);

void Yxy_2_RGB(cv::Mat& img);

void RGB_2_XYZ(cv::Mat& img);


void sRGB_2_RGB(cv::Mat& img);
void RGB_2_sRGB(cv::Mat& img);
void RGB_2_xyY(cv::Mat& img);

double norm2(std::vector<double> vec1, std::vector<double> vec2);