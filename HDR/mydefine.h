#pragma once
#include <stdlib.h>
#include <string>
#include <opencv2/opencv.hpp>
////====================================================================
////========================struct======================================
////====================================================================
struct Filedata
{
	std::string img_name;
	cv::Mat str;
};
struct FilePath
{
	std::string input_images_path;
	std::string output_images_path;
	std::string temp_images_path;

};
