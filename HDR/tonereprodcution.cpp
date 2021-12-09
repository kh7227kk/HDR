#include "tonereproduction.h"
#include "imageprocessing.h"
//cv::Mat tonereproduction(cv::Mat img)
//{
//	cv::Mat result,img_YCbCr;
//	std::vector<cv::Mat> channel;
//	cv::cvtColor(img, img_YCbCr, cv::COLOR_BGR2YCrCb);
//	cv::split(img_YCbCr, channel);
//	double Y_avg = 0, Y_W, Y_W_pow;
//	cv::Mat Y_S(img.rows, img.cols, CV_64FC1);
//	cv::Mat Y_T(img.rows, img.cols, CV_64FC1);
//	cv::Mat	Y_N(img.rows, img.cols, CV_8UC1);
//	double delta = 0.1;
//	double alpha, alpha_power;
//	double minValue, maxValue;
//	cv::Point minLoc, maxLoc;
//	for (int col = 0; col < img.cols; col++)
//	{
//		for (int row = 0; row < img.rows; row++)
//		{
//			Y_avg += log(img_YCbCr.at<cv::Vec3b>(row, col)[0] + delta);
//		}
//	}
//	Y_avg = exp(Y_avg / (img.rows * img.cols));
//
//	cv::minMaxLoc(channel[0], &minValue, &maxValue, &minLoc, &maxLoc);
//	if (maxValue > 252)	maxValue = 252;
//	if (minValue < 3)	minValue = 3;
//
//	alpha_power = 2 * (2 * log2(Y_avg) - log2(minValue) - log2(maxValue)) / (log2(maxValue) - log2(minValue));
//	alpha = 0.18 * pow(2, alpha_power);
//	Y_W_pow = log2(maxValue) - log2(minValue) - 8 / 2;
//	Y_W = 1.5 * pow(2, Y_W_pow);
//
//	for (int col = 0; col < img.cols; col++)
//	{
//		for (int row = 0; row < img.rows; row++)
//		{
//			Y_S.at<double>(row, col) = alpha * img_YCbCr.at<cv::Vec3b>(row, col)[0] / Y_avg;
//			Y_T.at<double>(row, col) = Y_S.at<double>(row, col) * (1 + (Y_S.at<double>(row, col) / pow(Y_W, 2))) / (1 + Y_S.at<double>(row, col));
//		}
//	}
//	double YTminValue, YTmaxValue;
//	cv::Point YTminLoc, YTmaxLoc;
//	cv::minMaxLoc(Y_T, &YTminValue, &YTmaxValue, &YTminLoc, &YTmaxLoc);
//	for (int col = 0; col < img.cols; col++)
//	{
//		for (int row = 0; row < img.rows; row++)
//		{
//			//std::cout << (Y_T.at<double>(row, col) - YTminValue) * 256 / (YTmaxValue - YTminValue) <<std::endl;
//			Y_N.at<uchar>(row, col) = int((Y_T.at<double>(row, col) - YTminValue) * 255 / (YTmaxValue - YTminValue));
//		}
//	}
//	cv::minMaxLoc(Y_N, &YTminValue, &YTmaxValue, &YTminLoc, &YTmaxLoc);
//	channel[0] = Y_N;
//	cv::merge(channel, img_YCbCr);
//	cv::cvtColor(img_YCbCr, result, cv::COLOR_YCrCb2BGR);
//	return result;
//}
//
//void imgprocessing(cv::Mat img, cv::Mat& out_img)
//{
//	cv::Mat str(img.rows, img.cols, CV_32FC3);
//	for (int y = 0; y < str.rows; y++)
//	{
//		for (int x = 0; x < str.cols; x++)
//		{
//			for (int c = 0; c < 3; c++)
//			{
//				str.at<cv::Vec3f>(y, x)[c] = (float)img.at<cv::Vec3b>(y, x)[c];
//			}
//		}
//	}
//	//tone reproduction and color enhancement
//	cv::Mat img_RL;
//	img_RL = reproduction_Lab(str);
//	out_img = img_RL.clone();
//}
//
//cv::Mat reproduction_Lab(cv::Mat img)
//{
//	cv::Mat result, img_Lab, img_L, img_LH, img_LL;
//	std::vector<cv::Mat> channel;
//	for (int y = 0; y < img.rows; y++)
//	{
//		for (int x = 0; x < img.cols; x++)
//		{
//			for (int c = 0; c < 3; c++)
//			{
//				img.at<cv::Vec3f>(y, x)[c] = img.at<cv::Vec3f>(y, x)[c] / 255.0;
//
//			}
//		}
//	}
//	cv::cvtColor(img, img_Lab, cv::COLOR_BGR2Lab);
//	cv::split(img_Lab, channel);
//	cv::Mat	L_N(img.rows, img.cols, CV_8UC1);
//	img_L = channel[0].clone();
//
//	cv::GaussianBlur(img_L, img_LL, cv::Size(13, 13), 0, 0);
//	img_LH = img_L - img_LL;
//	cv::Mat img_LL_reg = tonemapping_Lab(img_LL);
//	cv::Mat temp(img.rows, img.cols, CV_32FC1);
//	for (int y = 0; y < img.rows; y++)
//	{
//		for (int x = 0; x < img.cols; x++)
//		{
//			temp.at<float>(y, x) = img_LH.at<float>(y, x) + img_LL_reg.at<float>(y, x);
//		}
//	}
//	L_N = tonemapping_Lab(temp);
//
//	channel[0] = L_N;
//	channel[1] = channel[1] * 1.5;
//	channel[2] = channel[2] * 1.5;
//
//	cv::merge(channel, img_Lab);
//	cv::cvtColor(img_Lab, result, cv::COLOR_Lab2BGR);
//	result.convertTo(result, CV_8UC3, 255);
//	return result;
//}
//
//cv::Mat tonereproduction_Lab(cv::Mat img)
//{
//	cv::Mat result, img_Lab, img_L, img_LH, img_LL;
//	std::vector<cv::Mat> channel;
//	cv::cvtColor(img, img_Lab, cv::COLOR_BGR2Lab);
//	cv::split(img_Lab, channel);
//	cv::Mat	Y_N(img.rows, img.cols, CV_8UC1);
//	img_L = channel[0].clone();
//	cv::GaussianBlur(img_L, img_LL, cv::Size(3, 3), 0, 0);
//	img_LH = img_L - img_LL;
//
//	img_LL = tonemapping_Lab(img_LL);
//	Y_N = tonemapping_Lab(img_LH + img_LL);
//
//	channel[0] = Y_N;
//	cv::merge(channel, img_Lab);
//	cv::cvtColor(img_Lab, result, cv::COLOR_Lab2BGR);
//	return result;
//}
//
//
//cv::Mat tonemapping_Lab(cv::Mat img_L)
//{
//	int rowsize = img_L.rows;
//	int colsize = img_L.cols;
//	double L_avg = 0, L_W, L_W_pow;
//	cv::Mat L_S(rowsize, colsize, CV_32FC1);
//	cv::Mat L_T(rowsize, colsize, CV_32FC1);
//	cv::Mat	L_N(rowsize, colsize, CV_32FC1);
//	double delta = 1;
//	double alpha, alpha_power;
//	double minValue = 100, maxValue = 0;
//	cv::Point minLoc, maxLoc;
//	double Lmax = 100 - 0.01;
//	double Lmin = 0.01;
//
//	for (int col = 0; col < colsize; col++)
//	{
//		for (int row = 0; row < rowsize; row++)
//		{
//			if (img_L.at<float>(row, col) <= 0)
//			{
//				L_avg += 0;
//			}
//			else
//			{
//				L_avg += log(img_L.at<float>(row, col) + delta);
//			}
//		}
//	}
//	L_avg = exp(L_avg / (rowsize * colsize));
//	for (int col = 0; col < img_L.cols; col++)
//	{
//		for (int row = 0; row < img_L.rows; row++)
//		{
//			if ((img_L.at<float>(row, col) > Lmin && img_L.at<float>(row, col) < Lmax))
//			{
//				if (img_L.at<float>(row, col) > maxValue)
//				{
//					maxValue = img_L.at<float>(row, col);
//				}
//
//				if (img_L.at<float>(row, col) < minValue)
//				{
//					if (img_L.at<float>(row, col) == 0)
//					{
//						minValue = 1;
//					}
//					else
//					{
//						minValue = img_L.at<float>(row, col);
//					}
//				}
//			}
//		}
//	}
//
//	alpha_power = 2 * (2 * log2(L_avg) - log2(minValue) - log2(maxValue)) / (log2(maxValue) - log2(minValue));
//	alpha = 0.07 * pow(2, alpha_power);
//
//	L_W_pow = log2(maxValue) - log2(minValue) - log2(100) / 2;
//	L_W = 1.5 * pow(2, L_W_pow);
//	for (int col = 0; col < colsize; col++)
//	{
//		for (int row = 0; row < rowsize; row++)
//		{
//
//			L_S.at<float>(row, col) = alpha * img_L.at<float>(row, col) / L_avg;
//			L_T.at<float>(row, col) = L_S.at<float>(row, col) * (1 + (L_S.at<float>(row, col) / pow(L_W, 2))) / (1 + L_S.at<float>(row, col));
//		}
//	}
//	double LTminValue, LTmaxValue;
//	cv::Point LTminLoc, LTmaxLoc;
//	cv::minMaxLoc(L_T, &LTminValue, &LTmaxValue, &LTminLoc, &LTmaxLoc);
//	for (int col = 0; col < colsize; col++)
//	{
//		for (int row = 0; row < rowsize; row++)
//		{
//			L_N.at<float>(row, col) = (L_T.at<float>(row, col) - LTminValue) * 100 / (LTmaxValue - LTminValue);
//		}
//	}
//	//cv::minMaxLoc(L_N, &LTminValue, &LTmaxValue, &LTminLoc, &LTmaxLoc);
//
//	return L_N;
//
//}

////====================================================================
////========================edgedetcion=================================
////====================================================================
void edgedetcion(cv::Mat img, cv::Mat& edge_img)
{
	cv::Mat gray;
	cv::Mat src;
	int scale = 7;
	int delta = 0;
	int ddepth = CV_16S;

	cv::GaussianBlur(img, src, cv::Size(3, 3), 0, 0, cv::BORDER_DEFAULT);
	cv::cvtColor(src, gray, cv::COLOR_BGR2GRAY);
	cv::Mat grad, grad_x, grad_y;
	cv::Mat abs_grad_x, abs_grad_y;
	Sobel(gray, grad_x, ddepth, 1, 0, 3, scale, delta, cv::BORDER_DEFAULT);
	convertScaleAbs(grad_x, abs_grad_x);
	Sobel(gray, grad_y, ddepth, 0, 1, 3, scale, delta, cv::BORDER_DEFAULT);
	convertScaleAbs(grad_y, abs_grad_y);
	addWeighted(abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad);
	//cv::threshold(grad, grad, 127, 255, 0);
	edge_img = grad.clone();
}

////====================================================================
////========================macro_edge==================================
////====================================================================
void macro(cv::Mat img, cv::Mat& out_img)
{
	cv::Mat edge_img;
	cv::Mat	macro_img(img.rows, img.cols, CV_8UC1);
	edgedetcion(img, edge_img);
	//cv::imshow("edge", edge_img);
	for (int col = 0; col < macro_img.cols; col++)
	{
		for (int row = 0; row < macro_img.rows; row++)
		{
			if ((int)edge_img.at<uchar>(row, col) == 255)
			{
				for (int i = -4; i < 4; i++)
				{
					for (int j = -4; j < 4; j++)
					{
						if (row + i < macro_img.rows && col + j < macro_img.cols && row + i > 0 && col + j > 0)
						{
							macro_img.at<uchar>(row + i, col + j) = 255;
						}
					}
				}
			}
		}
	}

	macro_img.copyTo(out_img);
}

void macroEqualizationhis(cv::Mat img, cv::Mat& out_img)
{
	cv::Mat macro_img;
	macro(img, macro_img);
	std::vector<double> val(256, 0);
	std::vector<cv::Mat> channel_Lab(3);
	cv::Mat img_Lab;


	cv::cvtColor(img, img_Lab, cv::COLOR_BGR2Lab);
	cv::split(img_Lab, channel_Lab);
	double sum = 0;
	for (int x = 0; x < channel_Lab[0].cols; x++)
	{
		for (int y = 0; y < channel_Lab[0].rows; y++)
		{
			if ((int)macro_img.at<uchar>(y, x) == 255)
			{
				int temp = (int)channel_Lab[0].at<uchar>(y, x);
				val[temp] += 1;
				sum += 1;
			}
		}
	}
	std::vector<double>PDF(256);
	for (int P = 0; P < val.size(); P++)
	{
		PDF[P] = val[P] / sum;
	}
	std::vector<double>CDF(256);
	//CDF[0] = PDF[0] * 255.0;
	CDF[0] = PDF[0];

	for (int C = 1; C < 256; C++)
	{
		CDF[C] = CDF[C - 1] + PDF[C];
		//CDF[C] = CDF[C - 1] + PDF[C]*255.0;
	}
	//for (int C = 0; C < 256; C++)
	//{
	//	CDF[C] = round(CDF[C]);
	//}
	for (int i = 1; i < PDF.size(); i++)
	{
		if (CDF[i] / i * 255.0 > 2)
			CDF[i] = 2 * i * 0.8;
		else
			CDF[i] = round(CDF[i] * 255.0 * 0.8);

	}



	for (int y = 0; y < channel_Lab[0].rows; y++)
	{
		for (int x = 0; x < channel_Lab[0].cols; x++)
		{
			int temp = (int)channel_Lab[0].at<uchar>(y, x);
			channel_Lab[0].at<uchar>(y, x) = CDF[temp];
		}
	}


	cv::merge(channel_Lab, img_Lab);
	cv::cvtColor(img_Lab, img_Lab, cv::COLOR_Lab2BGR);
	img_Lab.convertTo(img_Lab, CV_32FC3, 1 / 255.0);
	cv::cvtColor(img_Lab, img_Lab, cv::COLOR_BGR2Lab);
	cv::split(img_Lab, channel_Lab);
	channel_Lab[1] = channel_Lab[1] * 1.2;
	channel_Lab[2] = channel_Lab[2] * 1.2;
	cv::merge(channel_Lab, img_Lab);
	cv::cvtColor(img_Lab, img_Lab, cv::COLOR_Lab2BGR);
	img_Lab.convertTo(img_Lab, CV_32FC3, 255.0);
	img_Lab.convertTo(img_Lab, CV_8UC3);
	out_img = img_Lab;
}



void macroEqualizationhis_Yuv(cv::Mat img, cv::Mat& out_img)
{
	cv::Mat macro_img;
	macro(img, macro_img);
	std::vector<double> val(256, 0);
	std::vector<cv::Mat> channel_Lab(3);
	cv::Mat img_Lab;


	cv::cvtColor(img, img_Lab, cv::COLOR_BGR2YUV);
	cv::split(img_Lab, channel_Lab);
	double sum = 0;
	for (int x = 0; x < channel_Lab[0].cols; x++)
	{
		for (int y = 0; y < channel_Lab[0].rows; y++)
		{
			if ((int)macro_img.at<uchar>(y, x) == 255)
			{
				int temp = (int)channel_Lab[0].at<uchar>(y, x);
				val[temp] += 1;
				sum += 1;
			}
		}
	}
	std::vector<double>PDF(256);
	for (int P = 0; P < val.size(); P++)
	{
		PDF[P] = val[P] / sum;
	}
	std::vector<double>CDF(256);
	CDF[0] = PDF[0] * 255.0;

	for (int C = 1; C < 256; C++)
	{
		CDF[C] = CDF[C - 1] + PDF[C] * 255.0;
	}
	for (int C = 0; C < 256; C++)
	{
		CDF[C] = round(CDF[C]);
	}
	for (int y = 0; y < channel_Lab[0].rows; y++)
	{
		for (int x = 0; x < channel_Lab[0].cols; x++)
		{
			int temp = (int)channel_Lab[0].at<uchar>(y, x);
			channel_Lab[0].at<uchar>(y, x) = CDF[temp];

			if (channel_Lab[1].at<uchar>(y, x) > 128)
			{

				if (channel_Lab[1].at<uchar>(y, x) * 1.05 > 255)
				{
					channel_Lab[1].at<uchar>(y, x) = channel_Lab[1].at<uchar>(y, x);
				}
				else
				{
					channel_Lab[1].at<uchar>(y, x) = channel_Lab[1].at<uchar>(y, x) * 1.05;
				}
			}
			else
			{
				channel_Lab[1].at<uchar>(y, x) = channel_Lab[1].at<uchar>(y, x) / 1.05;
			}

			if (channel_Lab[2].at<uchar>(y, x) > 128)
			{
				if (channel_Lab[2].at<uchar>(y, x) * 1.05 > 255)
				{
					channel_Lab[2].at<uchar>(y, x) = channel_Lab[2].at<uchar>(y, x);
				}
				else
				{
					channel_Lab[2].at<uchar>(y, x) = channel_Lab[2].at<uchar>(y, x) * 1.05;
				}
			}
			else
			{
				channel_Lab[2].at<uchar>(y, x) = channel_Lab[2].at<uchar>(y, x) / 1.05;
			}
		}
	}



	cv::merge(channel_Lab, img_Lab);
	cv::cvtColor(img_Lab, img_Lab, cv::COLOR_YUV2BGR);

	out_img = img_Lab;
}

void macroEqualizationhis_EPD(cv::Mat img, cv::Mat& out_img)
{
	cv::Mat macro_img;
	macro(img, macro_img);
	std::vector<double> val(256, 0);
	std::vector<cv::Mat> channel_Lab(3);
	cv::Mat img_Lab;


	cv::cvtColor(img, img_Lab, cv::COLOR_BGR2Lab);
	cv::split(img_Lab, channel_Lab);
	double sum = 0;
	for (int x = 0; x < channel_Lab[0].cols; x++)
	{
		for (int y = 0; y < channel_Lab[0].rows; y++)
		{
			if ((int)macro_img.at<uchar>(y, x) == 255)
			{
				int temp = (int)channel_Lab[0].at<uchar>(y, x);
				val[temp] += 1;
				sum += 1;
			}
		}
	}
	std::vector<double>PDF(256);
	for (int P = 0; P < val.size(); P++)
	{
		PDF[P] = val[P] / sum;
	}
	std::vector<double>CDF(256);
	//CDF[0] = PDF[0] * 255.0;
	CDF[0] = PDF[0];

	for (int C = 1; C < 256; C++)
	{
		CDF[C] = CDF[C - 1] + PDF[C];
	}
	//for (int C = 0; C < 256; C++)
	//{
	//	CDF[C] = round(CDF[C]);
	//}
	for (int i = 1; i < PDF.size(); i++)
	{
		if (CDF[i] / i * 255.0 > 2)
			CDF[i] = 2 * i * 0.8;
		else
			CDF[i] = round(CDF[i] * 255.0 * 0.8);

	}


	for (int y = 0; y < channel_Lab[0].rows; y++)
	{
		for (int x = 0; x < channel_Lab[0].cols; x++)
		{
			int temp = (int)channel_Lab[0].at<uchar>(y, x);
			channel_Lab[0].at<uchar>(y, x) = CDF[temp];
		}
	}
	cv::merge(channel_Lab, img_Lab);
	cv::cvtColor(img_Lab, img_Lab, cv::COLOR_Lab2BGR);
	img_Lab.convertTo(img_Lab, CV_32FC3, 1 / 255.0);
	cv::cvtColor(img_Lab, img_Lab, cv::COLOR_BGR2Lab);
	cv::split(img_Lab, channel_Lab);
	channel_Lab[1] = channel_Lab[1] * 1.5;
	channel_Lab[2] = channel_Lab[2] * 1.5;
	cv::merge(channel_Lab, img_Lab);
	cv::cvtColor(img_Lab, img_Lab, cv::COLOR_Lab2BGR);
	img_Lab.convertTo(img_Lab, CV_32FC3, 255.0);
	img_Lab.convertTo(img_Lab, CV_8UC3, 240.0 / 255.0);
	out_img = img_Lab;
}


void macroEqualizationhis_RGB(cv::Mat img, cv::Mat& out_img)
{
	std::vector<double> matrix(9);
	std::vector<double> vec(6);
	std::fstream fin;
	fin.open("./pso/CFA_results.txt", std::ios::in);
	for (int i = 0; i < 6; i++)
	{
		fin >> vec[i];
	}
	matrix[0] = vec[0];
	matrix[1] = vec[1];
	matrix[2] = 1 - vec[0] - vec[1];
	matrix[3] = vec[2];
	matrix[4] = vec[3];
	matrix[5] = 1 - vec[2] - vec[3];
	matrix[6] = vec[4];
	matrix[8] = vec[5];
	matrix[7] = 1 - vec[4] - vec[5];
	img.convertTo(out_img, CV_32FC3);
	cv::Mat img_cc = colorcorrection_one(out_img, matrix); /// 0~255
	img_cc.convertTo(out_img, CV_8UC3);

	cv::Mat macro_img;
	macro(img, macro_img);
	std::vector<double> val(256, 0);
	std::vector<cv::Mat> channel_Lab(3);
	cv::Mat img_Lab;

	cv::cvtColor(out_img, img_Lab, cv::COLOR_BGR2Lab);
	cv::split(img_Lab, channel_Lab);
	double sum = 0;
	for (int x = 0; x < channel_Lab[0].cols; x++)
	{
		for (int y = 0; y < channel_Lab[0].rows; y++)
		{
			if ((int)macro_img.at<uchar>(y, x) == 255)
			{
				int temp = (int)channel_Lab[0].at<uchar>(y, x);
				val[temp] += 1;
				sum += 1;
			}
		}
	}
	std::vector<double>PDF(256);
	for (int P = 0; P < val.size(); P++)
	{
		PDF[P] = val[P] / sum;
	}
	std::vector<double>CDF(256);
	//CDF[0] = PDF[0] * 255.0;
	CDF[0] = PDF[0];

	for (int C = 1; C < 256; C++)
	{
		CDF[C] = CDF[C - 1] + PDF[C];
		//CDF[C] = CDF[C - 1] + PDF[C]*255.0;
	}
	//for (int C = 0; C < 256; C++)
	//{
	//	CDF[C] = round(CDF[C]);
	//}
	for (int i = 1; i < PDF.size(); i++)
	{
		if (CDF[i] / i * 255.0 > 2)
			CDF[i] = 2 * i * 0.8;
		else
			CDF[i] = round(CDF[i] * 255.0 * 0.8);

	}



	for (int y = 0; y < channel_Lab[0].rows; y++)
	{
		for (int x = 0; x < channel_Lab[0].cols; x++)
		{
			for (int channel = 0; channel < 3; channel++)
			{
				int temp = (int)out_img.at<cv::Vec3b>(y, x)[channel];
				out_img.at<cv::Vec3b>(y, x)[channel] = CDF[temp];
			}

		}
	}

	cv::Mat out_lab;
	std::vector<cv::Mat> out_lab_ch;
	out_img.convertTo(out_img, CV_32FC3, 1 / 255.0);
	cv::cvtColor(out_img, out_lab, cv::COLOR_BGR2Lab);
	cv::split(out_lab, out_lab_ch);
	out_lab_ch[1] *= 1;
	out_lab_ch[2] *= 1;
	cv::merge(out_lab_ch, out_lab);
	cv::cvtColor(out_lab, out_img, cv::COLOR_Lab2BGR);
	out_img.convertTo(out_img, CV_8UC3, 255);
	//cv::merge(channel_Lab, img_Lab);
	//cv::cvtColor(img_Lab, img_Lab, cv::COLOR_Lab2BGR);
	//img_Lab.convertTo(img_Lab, CV_32FC3, 1 / 255.0);
	//cv::cvtColor(img_Lab, img_Lab, cv::COLOR_BGR2Lab);
	//cv::split(img_Lab, channel_Lab);
	//channel_Lab[1] = channel_Lab[1] * 1.2;
	//channel_Lab[2] = channel_Lab[2] * 1.2;
	//cv::merge(channel_Lab, img_Lab);
	//cv::cvtColor(img_Lab, img_Lab, cv::COLOR_Lab2BGR);
	//img_Lab.convertTo(img_Lab, CV_32FC3, 255.0);
	//img_Lab.convertTo(img_Lab, CV_8UC3);
	//out_img = img_Lab;



}