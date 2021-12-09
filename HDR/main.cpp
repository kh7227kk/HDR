/**********************************************************
Name :
Date : 2020/01/29
By   : Sean Chen
Final: 2020/02/
**********************************************************/
#include <iostream>
#include <stdio.h> 
//#include <stdlib.h>
//#include <opencv2/opencv.hpp>
#include <algorithm>
#include <vector>
#include <fstream>
//#include <string>
#include "mydefine.h"
#include "tonereproduction.h"
#include "imageprocessing.h"
//#include <io.h>

////====================================================================
////========================read========================================
////====================================================================
void read(std::string input_images_path, std::vector<Filedata>& img_data)
{
	std::vector<cv::String> files_fomat;
	std::vector<cv::String> all_files;
	std::vector<cv::String> files;
	files_fomat.push_back("/*.bmp");
	files_fomat.push_back("/*.png");
	files_fomat.push_back("/*.jpeg");
	files_fomat.push_back("/*.jpg");
	files_fomat.push_back("/*.tiff");
	files_fomat.push_back("/*.TIFF");
	files_fomat.push_back("/*.TIF");

	Filedata data;
	std::string filepath;
	int files_number = 0;
	for (int fomat_number = 0; fomat_number < files_fomat.size(); fomat_number++)
	{
		filepath = input_images_path + files_fomat[fomat_number];
		std::cout << filepath << std::endl;
		glob(filepath, files, false);
		files_number += files.size();
		std::cout << files_number << std::endl;
		for (int i = 0; i < files.size(); i++)
		{
			all_files.push_back(files[i]);
		}
		files.clear();
	}


	std::vector<cv::Mat> img;
	for (int img_number = 0; img_number < files_number; img_number++)
	{

		int length_sub = strlen(all_files[img_number].c_str()) - strlen(input_images_path.c_str());
		int length = strlen(input_images_path.c_str());
		std::string filesname = all_files[img_number].substr(length + 1, length_sub);
		img.push_back(imread(all_files[img_number]));
		std::cout << all_files[img_number] << std::endl;
		std::cout << filesname << std::endl;

		data.img_name = filesname;
		data.str = img[img_number];
		img_data.push_back(data);

	}
}

void read2(const FilePath& path, std::vector<Filedata>& img_data)
{
	std::vector<cv::Mat> img;
	std::vector<cv::String> fomat, all_files_fomat, all_files, files;
	fomat.push_back("/*.bmp");
	fomat.push_back("/*.png");
	fomat.push_back("/*.jpeg");
	fomat.push_back("/*.jpg");
	fomat.push_back("/*.tiff");
	Filedata data;
	std::string filepath;
	int files_number = 0;
	for (int fomat_number = 0; fomat_number < fomat.size(); fomat_number++)
	{
		filepath = path.input_images_path + fomat[fomat_number];
		std::cout << filepath << std::endl;
		glob(filepath, files, false);
		files_number += files.size();
		for (int i = 0; i < files.size(); i++)
		{
			all_files_fomat.push_back(fomat[fomat_number].substr(2));
			all_files.push_back(files[i]);
		}
		files.clear();
	}

	for (int img_number = 0; img_number < files_number; img_number++)
	{
		int length_sub = strlen(all_files[img_number].c_str()) - strlen(path.input_images_path.c_str()) - strlen(all_files_fomat[img_number].c_str()) - 1;
		int length = strlen(path.input_images_path.c_str());
		std::string filesname = all_files[img_number].substr(length + 1, length_sub);
		img.push_back(imread(all_files[img_number]));
		data.img_name = filesname;
		data.str = img[img_number];
		img_data.push_back(data);

	}
}


void data(std::string input_path, std::vector< Filedata>img_data)
{
	std::vector<cv::Mat> out;
	int level_size = 16;
	int ed_mode = 1;
	std::string output_path = "./merge/";
	for (int number = 0; number < img_data.size(); number++)
	{
		cv::Mat buff(1072, 1448, CV_8UC3);
		cv::Mat buff1;
		cv::Mat buff2;
		cv::Mat buff3;
		//imgprocessing(img_data[number].str, buff1);
		imgprocessing(img_data[number], buff1, level_size, ed_mode);
		imgprocessing_v3(img_data[number].str, buff2, level_size, ed_mode);

		//macroEqualizationhis(img_data[number].str, buff2);
		macroEqualizationhis(buff2, buff3);

		//macroEqualizationhis_Yuv(img_data[number].str, buff3);

		for (int y = 0; y < buff1.rows; y++)
		{
			for (int x = 0; x < buff1.cols; x++)
			{
				for (int c = 0; c < 3; c++)
				{
					buff.at<cv::Vec3b>(y, x)[c] = img_data[number].str.at<cv::Vec3b>(y, x)[c];
				}
			}
		}

		for (int y = 0; y < buff1.rows; y++)
		{
			for (int x = buff1.cols; x < buff.cols; x++)
			{
				for (int c = 0; c < 3; c++)
				{
					if (x - buff1.cols < buff1.cols)
					{
						buff.at<cv::Vec3b>(y, x)[c] = buff1.at<cv::Vec3b>(y, x - buff1.cols)[c];
					}
				}
			}
		}
		
		for (int y = buff1.rows; y < buff.rows; y++)
		{
			for (int x = 0; x < buff1.cols; x++)
			{
				for (int c = 0; c < 3; c++)
				{
					if (y - buff1.rows < buff1.rows)
					{
						buff.at<cv::Vec3b>(y, x)[c] = buff2.at<cv::Vec3b>(y - buff1.rows, x )[c];
					}
				}
			}
		}
				
		for (int y = buff1.rows; y < buff.rows; y++)
		{
			for (int x = buff1.cols; x < buff.cols; x++)
			{
				for (int c = 0; c < 3; c++)
				{
					if (x - buff1.cols < buff1.cols && y - buff1.rows < buff1.rows)
					{
						buff.at<cv::Vec3b>(y, x)[c] = buff3.at<cv::Vec3b>(y - buff1.rows, x - buff1.cols)[c];
					}
				}
			}
		}

		out.push_back(buff);
		//cv::imshow("a", buff);
		//cv::waitKey(0);
		std::string output_file_path = output_path + img_data[number].img_name;
		cv::imwrite(output_file_path, out[number]);

	}



}
void bmp(std::string input_path, std::vector< Filedata>img_data)
{
	std::string output_path = "./bmp/";

	for (int number = 0; number < img_data.size(); number++)
	{
		std::string output_file_path = output_path + img_data[number].img_name ;
		cv::imwrite(output_file_path, img_data[number].str);
	}
}

int main()
{
	std::vector< Filedata> img_data;
	std::string input_path = "./data";
	std::string output_path = "./tone/";
	std::string output_path2 = "./macro/";
	std::string output_path3 = "./macro2/";
	//read(input_path, img_data);

	FilePath path;
	//std::vector< Filedata> img_data;
	path.input_images_path = "./img";///kodim
	path.output_images_path = "./out/";
	path.temp_images_path = "./temp/";

	read(input_path, img_data);
	//bmp(input_path, img_data);
	//data(input_path, img_data);

	//system("pause");
	std::vector<cv::Mat> out;
	std::vector<cv::Mat> out2;
	std::vector<cv::Mat> out3;
	if (img_data.size() == 0)
	{
		std::cout << "there are no files or path is wrong";
	}
	else
	{
		for (int number = 0; number < img_data.size(); number++)
		{
			cv::Mat buff;
			cv::Mat buff2;
			cv::Mat buff3;
			//imgprocessing(img_data[number].str, buff);
			imgprocessing(img_data[number], buff, 16, 1);
			macroEqualizationhis(img_data[number].str, buff2);
			//macroEqualizationhis_Yuv(img_data[number].str, buff3);
			macroEqualizationhis_RGB(img_data[number].str, buff3);

			out.push_back(buff);
			out2.push_back(buff2);
			out3.push_back(buff3);
		}
		std::system("del /q tone\\");
		std::system("del /q macro\\");
		for (int number = 0; number < img_data.size(); number++)
		{
			cv::namedWindow("input"+std::to_string(number), 0);
			cv::imshow("input" + std::to_string(number), img_data[number].str);
			cv::resizeWindow("input" + std::to_string(number), 640, 480);
			cv::moveWindow("input" + std::to_string(number), 0, 0);

			cv::namedWindow("tone" + std::to_string(number), 0);
			cv::imshow("tone" + std::to_string(number), out[number]);
			cv::resizeWindow("tone" + std::to_string(number), 640, 480);
			cv::moveWindow("tone" + std::to_string(number), 640, 0);

			cv::namedWindow("macro" + std::to_string(number), 0);
			cv::imshow("macro" + std::to_string(number), out2[number]);
			cv::resizeWindow("macro" + std::to_string(number), 640, 480);
			cv::moveWindow("macro" + std::to_string(number), 1280, 0);

			cv::namedWindow("macro2_" + std::to_string(number), 0);
			cv::imshow("macro2_" + std::to_string(number), out3[number]);
			cv::resizeWindow("macro2_" + std::to_string(number), 640, 480);
			cv::moveWindow("macro2_" + std::to_string(number), 0, 480);

			//cv::moveWindow("output", 640, 0);
			
			cv::waitKey(0);
			std::string output_file_path = output_path + img_data[number].img_name;
			std::string output_file_path2 = output_path2 + img_data[number].img_name;
			std::string output_file_path3 = output_path3 + img_data[number].img_name;
			cv::imwrite(output_file_path, out[number]);
			cv::imwrite(output_file_path2, out2[number]);
			cv::imwrite(output_file_path3, out3[number]);
		}
	}
	return 0;
}
