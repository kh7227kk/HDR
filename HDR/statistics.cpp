
#include "statistics.h"

void Histogram(cv::Mat dst, int x, int y, int col_size, int row_size, std::vector<int>& val)
{
	for (int col = 0 + x; col < x + col_size; col++)
	{
		for (int row = 0 + y; row < y + row_size; row++)
		{
			int temp = dst.at<uchar>(row, col);
			val[temp] += 1;
		}
	}
}

void mean_Histogram(std::vector<int> val, double& mean)
{
	double total = 0;
	double temp = 0;
	for (int n = 0; n < val.size(); n++)
	{
		total += val[n];
		temp += val[n] * n;
	}
	mean = temp / total;
}

void SD_Histogram(std::vector<int> val, double mean, double& SD)
{
	double total = 0;
	double temp = 0;
	for (int n = 0; n < val.size(); n++)
	{
		total += val[n];
		temp += val[n] * (n - mean) * (n - mean);
	}
	SD = sqrt(temp / total);

}
//
//void printHistogram(Filedata img)
//{
//
//	std::vector<std::vector<int>> val(3, std::vector<int>(256, 0));
//	const int size = 256;
//	std::vector<cv::Mat> channel;
//	std::vector<double> mean(3);
//	std::vector<double> SD_val(3);
//	cv::split(img.str, channel);
//	for (int i = 0; i < 3; i++)
//	{
//		Histogram(channel[i], 0, 0, img.str.cols, img.str.rows, val[i]);
//		mean_Histogram(val[i], mean[i]);
//		SD_Histogram(val[i], mean[i], SD_val[i]);
//	}
//	cv::Mat dst(size, size, CV_8UC3, cv::Scalar(0, 0, 0));
//	std::vector<cv::Mat> dstImage;
//	cv::split(dst, dstImage);
//	std::vector<double> LTminValue(3), LTmaxValue(3);
//	std::vector<cv::Point> LTminLoc(3), LTmaxLoc(3);
//	for (int n = 0; n < 3; n++)
//	{
//		int val_max = *max_element(val[n].begin(), val[n].end());
//		int val_min = *min_element(val[n].begin(), val[n].end());
//		for (int i = 0; i < 256; i++)
//		{
//			int temp = (val[n][i] - val_min) * 255.0 / (val_max - val_min);
//			cv::line(dstImage[n], cv::Point(i, size - 1), cv::Point(i, 255 - temp), cv::Scalar(255));
//		}
//
//		cv::minMaxLoc(channel[n], &LTminValue[n], &LTmaxValue[n], &LTminLoc[n], &LTmaxLoc[n]);
//	}
//
//	cv::Mat textbox(cv::Size(1100, 256), CV_8UC3);
//	//設定背景
//	textbox.setTo(cv::Scalar(0, 0, 0));
//	//設定繪製文字的相關引數
//	std::string text_R = "R mean:" + std::to_string(mean[2]) + " SD:" + std::to_string(SD_val[2]) + " max:" + std::to_string(LTmaxValue[2]) + " min:" + std::to_string(LTminValue[2]);
//	std::string text_G = "G mean:" + std::to_string(mean[1]) + " SD:" + std::to_string(SD_val[1]) + " max:" + std::to_string(LTmaxValue[1]) + " min:" + std::to_string(LTminValue[1]);
//	std::string text_B = "B mean:" + std::to_string(mean[0]) + " SD:" + std::to_string(SD_val[0]) + " max:" + std::to_string(LTmaxValue[0]) + " min:" + std::to_string(LTminValue[0]);
//
//	int font_face = cv::FONT_HERSHEY_SIMPLEX;
//	double font_scale = 1;
//	int thickness = 2;
//	int baseline;
//	//獲取文字框的長寬
//	cv::Size text_size = cv::getTextSize(text_R, font_face, font_scale, thickness, &baseline);
//	//將文字框居中繪製
//	cv::Point origin;
//	origin.x = textbox.cols / 2 - text_size.width / 2;
//	origin.y = textbox.rows / 2 - text_size.height / 2;
//	cv::putText(textbox, text_R, origin, font_face, font_scale, cv::Scalar(255, 255, 255), thickness, 1, 0);
//	origin.y = origin.y - 2 * text_size.height;
//	cv::putText(textbox, text_G, origin, font_face, font_scale, cv::Scalar(255, 255, 255), thickness, 1, 0);
//	origin.y = origin.y + 4 * text_size.height;
//	cv::putText(textbox, text_B, origin, font_face, font_scale, cv::Scalar(255, 255, 255), thickness, 1, 0);
//	//顯示繪製解果
//	cv::namedWindow("textbox_" + img.img_name, 0);
//	cv::resizeWindow("textbox_" + img.img_name, textbox.cols, textbox.rows);
//	cv::imshow("textbox_" + img.img_name, textbox);
//	cv::namedWindow("B_" + img.img_name, 0);
//	cv::namedWindow("G_" + img.img_name, 0);
//	cv::namedWindow("R_" + img.img_name, 0);
//	cv::namedWindow("image_" + img.img_name, 0);
//	cv::resizeWindow("R_" + img.img_name, dstImage[0].cols, dstImage[0].rows);
//	cv::resizeWindow("G_" + img.img_name, dstImage[0].cols, dstImage[0].rows);
//	cv::resizeWindow("B_" + img.img_name, dstImage[0].cols, dstImage[0].rows);
//	cv::imshow("R_" + img.img_name, dstImage[2]);
//	cv::imshow("G_" + img.img_name, dstImage[1]);
//	cv::imshow("B_" + img.img_name, dstImage[0]);
//
//	int limit_rows;
//	int limit_cols;
//	if (img.str.cols < 120 || img.str.rows < 100)
//	{
//
//		if (img.str.cols / 120.0 > img.str.rows / 100.0)
//		{
//			limit_rows = 100;
//			limit_cols = img.str.cols * (100.0 / img.str.rows);
//		}
//		else
//		{
//			limit_cols = 120;
//			limit_rows = img.str.rows * (120.0 / img.str.cols);
//
//		}
//
//		cv::resizeWindow("image_" + img.img_name, limit_cols, img.str.rows);
//		cv::imshow("image_" + img.img_name, img.str);
//		cv::moveWindow("textbox_" + img.img_name, 0, 256);
//		cv::moveWindow("B_" + img.img_name, 125 + 2 * dstImage[0].cols, 0);
//		cv::moveWindow("G_" + img.img_name, 125 + dstImage[0].cols, 0);
//		cv::moveWindow("R_" + img.img_name, 125, 0);
//		cv::moveWindow("image_" + img.img_name, 0, 0);
//	}
//	else
//	{
//		cv::resizeWindow("image_" + img.img_name, img.str.cols, img.str.rows);
//		cv::imshow("image_" + img.img_name, img.str);
//		cv::moveWindow("textbox_" + img.img_name, img.str.cols, 256);
//		cv::moveWindow("B_" + img.img_name, 125 + 2 * dstImage[0].cols, 0);
//		cv::moveWindow("G_" + img.img_name, 125 + dstImage[0].cols, 0);
//		cv::moveWindow("R_" + img.img_name, 125, 0);
//		cv::moveWindow("image_" + img.img_name, 0, 0);
//	}
//
//
//	cv::waitKey(0);
//	//cv::destroyAllWindows();
//	cv::imwrite("./printHistogram/" + img.img_name + "_str.jpg", img.str);
//	cv::imwrite("./printHistogram/" + img.img_name + "_B.jpg", dstImage[0]);
//	cv::imwrite("./printHistogram/" + img.img_name + "_G.jpg", dstImage[1]);
//	cv::imwrite("./printHistogram/" + img.img_name + "_R.jpg", dstImage[2]);
//
//}


void Regression(std::string path, const int p, std::vector<double>& coeff)
{
	std::fstream fin;
	fin.open(path, std::ios::in);
	std::vector<double> vec_y;
	std::vector<double> vec_x;
	double reg_x;
	double reg_y;
	while (!fin.eof())
	{
		fin >> reg_x;
		fin >> reg_y;
		vec_x.push_back(reg_x);
		vec_y.push_back(reg_y);
		std::cout << reg_x << "\t" << reg_y << std::endl;
	}
	fin.close();
	const int vec_size = vec_y.size();
	cv::Mat matrix(vec_size, p, CV_64FC1);
	cv::Mat inv_matrix(vec_size, vec_size, CV_64FC1);
	cv::Mat transpose_matrix(p, vec_size, CV_64FC1);
	cv::Mat coeffvector(p, 1, CV_64FC1);
	cv::Mat temp(vec_size, 1, CV_64FC1);
	cv::Mat varible(vec_size, 1, CV_64FC1);
	for (int y = 0; y < vec_size; y++)
	{
		varible.at<double>(y, 0) = vec_y[y];
	}


	for (int y = 0; y < vec_size; y++)
	{
		for (int x = 0; x < p; x++)
		{
			matrix.at<double>(y, x) = pow(vec_x[y], x);
			std::cout << matrix.at<double>(y, x) << "\t";
		}
		std::cout << std::endl;
	}

	cv::transpose(matrix, transpose_matrix);
	inv_matrix = (transpose_matrix * matrix).inv() * transpose_matrix;
	coeffvector = inv_matrix.clone() * varible.clone();
	for (int y = 0; y < p; y++)
	{
		std::cout << coeffvector.at<double>(y, 0) << "\t";
		coeff.push_back(coeffvector.at<double>(y, 0));
	}
	std::cout << std::endl;

}

////====================================================================
////========================Color space transform===================================
////====================================================================
void sRGB_2_RGB(cv::Mat& img)
{
	double ThrLin2Gamma = 0.04045;
	for (int r = 0; r < img.rows; r++)
	{
		for (int c = 0; c < img.cols; c++)
		{
			for (int channel = 0; channel < 3; channel++)
			{
				double temp = (int)img.at<cv::Vec3f>(r, c)[channel] / 255.0;
				if (temp > ThrLin2Gamma)
				{
					img.at<cv::Vec3f>(r, c)[channel] = pow(((temp + 0.055) / 1.055), 2.4) * 255;
				}
				else
				{
					img.at<cv::Vec3f>(r, c)[channel] = temp / 12.92;
				}
			}
		}
	}
}

void RGB_2_sRGB(cv::Mat& img)
{
	double ThrLin2Gamma = 0.00304;
	for (int r = 0; r < img.rows; r++)
	{
		for (int c = 0; c < img.cols; c++)
		{
			for (int channel = 0; channel < 3; channel++)
			{
				double temp = (int)img.at<cv::Vec3b>(r, c)[channel] / 255.0;
				if (temp > ThrLin2Gamma)
				{
					img.at<cv::Vec3b>(r, c)[channel] = ((1 + 0.055) * pow(temp, 1 / 2.44) - 0.055) * 255;
				}
				else
				{
					img.at<cv::Vec3b>(r, c)[channel] = 12.92 * temp;
				}
			}
		}
	}
}

void RGB_2_xyY(cv::Mat& img) {
	cv::Mat data_XYZ;
	cv::cvtColor(img, data_XYZ, cv::COLOR_BGR2XYZ);
	for (int r = 0; r < img.rows; r++)
	{
		for (int c = 0; c < img.cols; c++)
		{
			double reg = data_XYZ.at<cv::Vec3f>(r, c)[0] + data_XYZ.at<cv::Vec3f>(r, c)[1] + data_XYZ.at<cv::Vec3f>(r, c)[2];
			img.at<cv::Vec3f>(r, c)[0] = data_XYZ.at<cv::Vec3f>(r, c)[0] / reg; //x
			img.at<cv::Vec3f>(r, c)[1] = data_XYZ.at<cv::Vec3f>(r, c)[1] / reg; //y
			img.at<cv::Vec3f>(r, c)[2] = data_XYZ.at<cv::Vec3f>(r, c)[1] * 100.0; //Y
		}
	}
}

////====================================================================
////========================YxytoRGB====================================
////====================================================================
void Yxy_2_RGB(cv::Mat& img) {
	cv::Mat data_XYZ(img.rows, img.cols, CV_32FC3);
	for (int r = 0; r < img.rows; r++)
	{
		for (int c = 0; c < img.cols; c++)
		{
			data_XYZ.at<cv::Vec3f>(r, c)[1] = img.at<cv::Vec3f>(r, c)[2];//Y
			data_XYZ.at<cv::Vec3f>(r, c)[0] = img.at<cv::Vec3f>(r, c)[2] *
				img.at<cv::Vec3f>(r, c)[0] / img.at<cv::Vec3f>(r, c)[1];//X
			data_XYZ.at<cv::Vec3f>(r, c)[2] = img.at<cv::Vec3f>(r, c)[2] *
				(1 - img.at<cv::Vec3f>(r, c)[0] - img.at<cv::Vec3f>(r, c)[1])
				/ img.at<cv::Vec3f>(r, c)[1];//y
		}
	}
	cv::cvtColor(data_XYZ, img, cv::COLOR_XYZ2BGR);
	for (int r = 0; r < img.rows; r++)
	{
		for (int c = 0; c < img.cols; c++)
		{
			for (int channel = 0; channel < 3; channel++)
			{
				if (img.at<cv::Vec3f>(r, c)[channel] > 255)
					img.at<cv::Vec3f>(r, c)[channel] = 255;
				else if (img.at<cv::Vec3f>(r, c)[channel] < 0)
					img.at<cv::Vec3f>(r, c)[channel] = 0;

			}
		}
	}
}

void RGB_2_XYZ(cv::Mat& img)
{
	cv::Mat data_XYZ = img.clone();
	cv::Mat RGB(3, 1, CV_32FC1);
	cv::Mat XYZ(3, 1, CV_32FC1);
	cv::Mat matrix(3, 3, CV_32FC1);
	matrix.at<float>(0, 0) = 0.49;
	matrix.at<float>(0, 1) = 0.31;
	matrix.at<float>(0, 2) = 0.20;
	matrix.at<float>(1, 0) = 0.17697;
	matrix.at<float>(1, 1) = 0.81240;
	matrix.at<float>(1, 2) = 0.01063;
	matrix.at<float>(2, 0) = 0.00;
	matrix.at<float>(2, 1) = 0.01;
	matrix.at<float>(2, 2) = 0.99;


	for (int r = 0; r < img.rows; r++)
	{
		for (int c = 0; c < img.cols; c++)
		{
			RGB.at<float>(2, 0) = img.at<cv::Vec3f>(r, c)[0];
			RGB.at<float>(1, 0) = img.at<cv::Vec3f>(r, c)[1];
			RGB.at<float>(0, 0) = img.at<cv::Vec3f>(r, c)[2];
			XYZ = (1 / 0.17697) * matrix * RGB;
			img.at<cv::Vec3f>(r, c)[0] = XYZ.at<float>(0, 0);
			img.at<cv::Vec3f>(r, c)[1] = XYZ.at<float>(1, 0);
			img.at<cv::Vec3f>(r, c)[2] = XYZ.at<float>(2, 0);
		}
	}
}

////====================================================================
////========================norm2=======================================
////====================================================================
double norm2(std::vector<double> vec1, std::vector<double> vec2)
{
	int vec_size = vec1.size();
	double val = 0;

	for (int i = 0; i < vec_size; i++)
	{
		val = val + pow(vec1[i] - vec2[i], 2);
	}

	return sqrt(val);

}
