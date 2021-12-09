#pragma once
#include <opencv2/opencv.hpp>
//
//cv::Mat tonereproduction(cv::Mat img);
//
//void imgprocessing(cv::Mat img, cv::Mat& out_img);
//cv::Mat reproduction_Lab(cv::Mat img);
//cv::Mat tonereproduction_Lab(cv::Mat img);
//cv::Mat tonemapping_Lab(cv::Mat img_L);

void edgedetcion(cv::Mat img, cv::Mat& edge_img);
void macro(cv::Mat img, cv::Mat& out_img);
void macroEqualizationhis(cv::Mat img, cv::Mat& out_img);
void macroEqualizationhis_Yuv(cv::Mat img, cv::Mat& out_img);
void macroEqualizationhis_EPD(cv::Mat img, cv::Mat& out_img);
void macroEqualizationhis_RGB(cv::Mat img, cv::Mat& out_img);
