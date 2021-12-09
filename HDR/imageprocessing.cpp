#include "imageprocessing.h"

////====================================================================
////========================imgprocessing===============================
////====================================================================
void imgprocessing(Filedata img, cv::Mat& out_img, int level_size, int ed_mode)
{
	cv::Mat str;
	cv::Mat img_wb;
	cv::Mat img_cc;
	cv::Mat img_RL;
	//cv::Mat img_ed;
	img.str.convertTo(str, CV_32FC3); ///range is 0 ~ 255
	//white balance
	//img_wb = balance_white(str); ///range R= (240 - 18.395) * x / 240 + 18.395; G = (210 - 2.2535) * x / 240 + 2.2535; B = (233 - 0) * x / 240 + 0;

	//color correction
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

	//img_cc = colorcorrection_one(str, matrix); /// 0~255
	img_cc = str.clone();
	//cv::Mat show;
	//img_cc.convertTo(show, CV_8UC3);
	//cv::imshow("rgb", show);
	//cv::waitKey(0);
	//tone reproduction and color enhancement
	img_RL = reproduction_Lab(img_cc, img);
	//error_diffusion
	cv::Mat img_ed(img.str.rows, img.str.cols, CV_32FC3);
	for (int y = 0; y < img.str.rows; y++)
	{
		for (int x = 0; x < img.str.cols; x++)
		{
			for (int c = 0; c < 3; c++)
			{
				img_ed.at<cv::Vec3f>(y, x)[c] = img_RL.at<cv::Vec3f>(y, x)[c] * 240;
			}
		}
	}
	std::vector<std::vector<int>> level(3, std::vector<int>(level_size));
	for (int i = 0; i < level_size; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (level_size == 2)
			{

				if (i == 0)
				{
					level[j][i] = 0;
				}
				else if (i == 1)
				{
					level[j][i] = 255;
				}
			}
			else
			{
				level[j][i] = i * 16;

				//level[j][i] = i * 256 / level_size;
				if (level[j][i] == 256)
					level[j][i] == 255;
			}


		}

	}
	img_ed.convertTo(out_img, CV_8UC3);
	//Error_diffusion_v2(img_ed, ed_mode, level);
	//img_ed.convertTo(out_img, CV_8UC3);
	//out_img = img_ed.clone();
}


void imgprocessing_v3(cv::Mat img, cv::Mat& out_img, int level_size, int ed_mode)
{
	cv::Mat str;
	cv::Mat img_wb;
	cv::Mat img_cc;
	cv::Mat img_RL;
	//cv::Mat img_ed;
	img.convertTo(str, CV_32FC3); ///range is 0 ~ 255
	//white balance
	//img_wb = balance_white(str); ///range R= (240 - 18.395) * x / 240 + 18.395; G = (210 - 2.2535) * x / 240 + 2.2535; B = (233 - 0) * x / 240 + 0;

	//color correction
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

	img_cc = colorcorrection_one(str, matrix); /// 0~255
	//cv::Mat show;
	//img_cc.convertTo(show, CV_8UC3);
	//cv::imshow("rgb", show);
	//cv::waitKey(0);
	//tone reproduction and color enhancement
	//img_RL = reproduction_Lab(img_cc, img);
	img_RL = img_cc.clone();
	//error_diffusion
	cv::Mat img_ed(img.rows, img.cols, CV_32FC3);
	for (int y = 0; y < img.rows; y++)
	{
		for (int x = 0; x < img.cols; x++)
		{
			for (int c = 0; c < 3; c++)
			{
				//img_ed.at<cv::Vec3f>(y, x)[c] = img_RL.at<cv::Vec3f>(y, x)[c] * 240;
				img_ed.at<cv::Vec3f>(y, x)[c] = img_RL.at<cv::Vec3f>(y, x)[c] * 240 / 250.0;
			}
		}
	}
	std::vector<std::vector<int>> level(3, std::vector<int>(level_size));
	for (int i = 0; i < level_size; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (level_size == 2)
			{

				if (i == 0)
				{
					level[j][i] = 0;
				}
				else if (i == 1)
				{
					level[j][i] = 255;
				}
			}
			else
			{
				//level[j][i] = i * 16;

				level[j][i] = i * 256 / level_size;
				if (level[j][i] == 256)
					level[j][i] == 255;
			}


		}

	}
	img_ed.convertTo(out_img, CV_8UC3);
	//Error_diffusion_v2(img_ed, ed_mode, level);
	//img_ed.convertTo(out_img, CV_8UC3);
	//out_img = img_ed.clone();
}

cv::Mat oversaturationset(cv::Mat img)
{
	cv::Mat out;
	out = img.clone();

	std::vector<double> img_sort;
	for (int y = 0; y < img.rows; y++)
	{
		for (int x = 0; x < img.cols; x++)
		{
			for (int channel = 0; channel < 3; channel++)
			{
				img_sort.push_back(img.at<cv::Vec3f>(y, x)[channel]);
				//std::cout << img.at<cv::Vec3f>(y, x)[channel] <<std::endl;
			}
		}
	}
	std::sort(img_sort.begin(), img_sort.end(), std::greater<double>());
	double th = img_sort[img.rows * img.cols * 3 / 100];
	for (int y = 0; y < img.rows; y++)
	{
		for (int x = 0; x < img.cols; x++)
		{
			for (int channel = 0; channel < 3; channel++)
			{
				if (out.at<cv::Vec3f>(y, x)[channel] > th)
				{
					out.at<cv::Vec3f>(y, x)[channel] = th;
				}
				else if (out.at<cv::Vec3f>(y, x)[channel] < 0)
				{
					out.at<cv::Vec3f>(y, x)[channel] = 0;
				}
			}
		}
	}
	std::vector<cv::Mat>img_channel;
	cv::split(out, img_channel);
	//////normal
	double LTminValue, LTmaxValue;
	cv::Point LTminLoc, LTmaxLoc;
	int max_val = 0;
	for (int i = 0; i < 3; i++)
	{
		cv::minMaxLoc(img_channel[i], &LTminValue, &LTmaxValue, &LTminLoc, &LTmaxLoc);

		if (max_val < LTmaxValue)
		{
			max_val = LTmaxValue;
		}

	}


	for (int col = 0; col < out.cols; col++)
	{
		for (int row = 0; row < out.rows; row++)
		{
			for (int i = 0; i < 3; i++)
			{
				img_channel[i].at<float>(row, col) = (img_channel[i].at<float>(row, col) - LTminValue) * 255.0 / (max_val - LTminValue);
			}
		}
	}

	cv::merge(img_channel, out);
	//cv::imshow("a",out );

	//cv::waitKey(0);
	return out;

}

void imgprocessing_v2(cv::Mat img, cv::Mat& out_img, int level_size, int ed_mode)
{
	cv::Mat img_ed = img.clone();
	std::vector<std::vector<int>> level(3, std::vector<int>(level_size));
	for (int i = 0; i < level_size; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			level[j][i] = i * 256 / level_size;
			if (level[j][i] == 256)
				level[j][i] == 255;

			//if (i == 0)
			//{
			//	level[j][i] = 0;
			//}
			//else if (i == 1)
			//{
			//	level[j][i] = 255;
			//}
		}

	}



	Error_diffusion_v2(img_ed, ed_mode, level);

	out_img = img_ed.clone();
}

////====================================================================
////========================reproduction_Lab============================
////====================================================================
cv::Mat reproduction_Lab(cv::Mat img, Filedata org_img)
{
	cv::Mat result, img_Lab, img_L, img_LH, img_LL;
	std::vector<cv::Mat> channel;
	for (int y = 0; y < img.rows; y++)
	{
		for (int x = 0; x < img.cols; x++)
		{
			for (int c = 0; c < 3; c++)
			{
				img.at<cv::Vec3f>(y, x)[c] = img.at<cv::Vec3f>(y, x)[c] / 255.0;
			}
		}
	}
	cv::cvtColor(img, img_Lab, cv::COLOR_BGR2Lab);
	cv::split(img_Lab, channel);
	//org_img is cv_8uc3

	cv::Mat	L_N(img.rows, img.cols, CV_8UC1);
	img_L = channel[0].clone();

	//cv::Mat show;
	//img_L.convertTo(show, CV_8UC3);
	//std::vector<int> val(100);
	//Histogram(show, 0, 0, show.cols, show.rows, val);

	//cv::GaussianBlur(img_L, img_LL, cv::Size(5, 5), 0, 0);
	cv::blur(img_L, img_LL, cv::Size(5, 5));
	img_LH = img_L - img_LL;

	for (int y = 0; y < img.rows; y++)
	{
		for (int x = 0; x < img.cols; x++)
		{
			img_LH.at<float>(y, x) = img_L.at<float>(y, x) - img_LL.at<float>(y, x);
		}
	}

	//img_LH.convertTo(img_LH, CV_32FC3, 1.0, -50);
	//img_LH.convertTo(img_LH, CV_32FC3, 0.5);



	cv::Mat img_LL_reg = tonemapping_Lab(img_LL, img_LL);
	cv::Mat temp(img.rows, img.cols, CV_32FC1);
	for (int y = 0; y < img.rows; y++)
	{
		for (int x = 0; x < img.cols; x++)
		{
			temp.at<float>(y, x) = img_LH.at<float>(y, x) + img_LL_reg.at<float>(y, x);
		}
	}


	L_N = tonemapping_Lab(temp, temp);

	//cv::Mat org_Lab;
	//org_img.convertTo(org_Lab, CV_32FC3, 1 / 255.0);
	//cv::cvtColor(org_Lab, org_Lab, cv::COLOR_BGR2Lab);
	//std::vector<cv::Mat> org_Lab_channel;
	//cv::split(org_Lab, org_Lab_channel);
	//cv::Mat org_img_L;
	//cv::Mat org_img_H = channel[0].clone();
	//org_img_H.convertTo(org_img_H, CV_32FC1);
	//cv::GaussianBlur(org_Lab_channel[0], org_img_L, cv::Size(5, 5), 0, 0);
	//for (int y = 0; y < org_img.rows; y++)
	//{
	//	for (int x = 0; x < org_img.cols; x++)
	//	{
	//		org_img_H.at<float>(y, x) = org_Lab_channel[0].at<float>(y, x) - org_img_L.at<float>(y, x);
	//	}
	//}
	////cv::Mat org_img_H_show= org_img_H.clone();
	double LTminValue, LTmaxValue;
	cv::Point LTminLoc, LTmaxLoc;
	//cv::minMaxLoc(org_img_H, &LTminValue, &LTmaxValue, &LTminLoc, &LTmaxLoc);
	//cv::Mat org_img_H_show = org_img_H.clone();
	//for (int col = 0; col < org_img.cols; col++)
	//{
	//	for (int row = 0; row < org_img.rows; row++)
	//	{
	//		org_img_H_show.at<float>(row, col) = (org_img_H_show.at<float>(row, col) - LTminValue) * 100 / (LTmaxValue - LTminValue);
	//	}
	//}
	//org_img_H_show.convertTo(org_img_H_show, CV_8UC1, 2.25);

	//cv::imshow("org_img_H_show", org_img_H_show);
	//cv::waitKey(0);

	//for (int y = 0; y < org_img.rows; y++)
	//{
	//	for (int x = 0; x < org_img.cols; x++)
	//	{
	//		L_N.at<float>(y, x) = L_N.at<float>(y, x) +img_LH.at<float>(y, x);
	//	}
	//}

	////////normal
	//cv::minMaxLoc(L_N, &LTminValue, &LTmaxValue, &LTminLoc, &LTmaxLoc);
	//std::cout << "L_N:" << std::endl;
	//std::cout << LTminValue << "\t" << LTmaxValue << std::endl;

	//for (int col = 0; col < L_N.cols; col++)
	//{
	//	for (int row = 0; row < L_N.rows; row++)
	//	{
	//		L_N.at<float>(row, col) = (L_N.at<float>(row, col) - LTminValue) * 100 / (LTmaxValue - LTminValue);
	//	}
	//}
	for (int y = 0; y < img.rows; y++)
	{
		for (int x = 0; x < img.cols; x++)
		{
			L_N.at<float>(y, x) = img_LH.at<float>(y, x) + L_N.at<float>(y, x);
		}
	}
	cv::minMaxLoc(L_N, &LTminValue, &LTmaxValue, &LTminLoc, &LTmaxLoc);


	//for (int col = 0; col < L_N.cols; col++)
	//{
	//	for (int row = 0; row < L_N.rows; row++)
	//	{
	//		L_N.at<float>(row, col) = (L_N.at<float>(row, col) - LTminValue) * 100 / (LTmaxValue - LTminValue);
	//	}
	//}


	channel[0] = L_N;
	cv::Mat show;
	L_N.convertTo(show, CV_8UC3);
	//cv::imshow("a", show);
	cv::imwrite("./test/"+org_img.img_name+".bmp",show);
	//cv::waitKey(0);

	channel[1] = channel[1]*1.2 ;
	channel[2] = channel[2]*1.2 ;
	cv::merge(channel, img_Lab);

	cv::cvtColor(img_Lab, result, cv::COLOR_Lab2BGR);
	return result;
}

////====================================================================
////========================balance_white===============================
////====================================================================
cv::Mat balance_white(cv::Mat img)
{
	//range R= (240 - 18.395) * x / 240 + 18.395; G = (210 - 2.2535) * x / 240 + 2.2535; B = (233 - 0) * x / 240 + 0;
	double R;
	double G;
	double B;
	cv::Mat out = img.clone();
	for (int y = 0; y < out.rows; y++)
	{
		for (int x = 0; x < out.cols; x++)
		{
			R = (240 - 18.395) * out.at<cv::Vec3f>(y, x)[2] / 240 + 18.395;
			G = (210 - 2.2535) * out.at<cv::Vec3f>(y, x)[1] / 240 + 2.2535;
			B = (233 - 0) * out.at<cv::Vec3f>(y, x)[0] / 240 + 0;

			out.at<cv::Vec3f>(y, x)[0] = B;
			out.at<cv::Vec3f>(y, x)[1] = G;
			out.at<cv::Vec3f>(y, x)[2] = R;
		}
	}
	return out;
}



////====================================================================
////========================tone reproduction lab=======================
////====================================================================
cv::Mat tonemapping_Lab(cv::Mat img_L, cv::Mat str)
{
	int rowsize = img_L.rows;
	int colsize = img_L.cols;
	double L_avg = 0, L_W, L_W_pow;
	cv::Mat L_S(rowsize, colsize, CV_32FC1);
	cv::Mat L_T(rowsize, colsize, CV_32FC1);
	cv::Mat	L_N(rowsize, colsize, CV_32FC1);
	double delta = 0.1;
	double alpha, alpha_power;
	double minValue = 100, maxValue = 0;
	cv::Point minLoc, maxLoc;
	double Lmax = 100 - 0.01;
	double Lmin = 0.01;
	//double LminValue, LmaxValue;
	cv::Point LTminLoc, LTmaxLoc;
	cv::minMaxLoc(img_L, &Lmin, &Lmax, &LTminLoc, &LTmaxLoc);
	std::cout << "img_L:" << Lmin << "\t" << Lmax << std::endl;

	img_L.convertTo(img_L, CV_32FC3, 1.0, (-Lmin ));
	cv::minMaxLoc(img_L, &Lmin, &Lmax, &LTminLoc, &LTmaxLoc);

	std::cout << "img_L:" << Lmin << "\t" << Lmax<<std::endl;

	for (int col = 0; col < colsize; col++)
	{
		for (int row = 0; row < rowsize; row++)
		{

			if (img_L.at<float>(row, col) <= -1)
			{
				L_avg += 0;
			}
			else
			{
				L_avg += log(img_L.at<float>(row, col) + delta);
			}
		}
	}
	L_avg = exp(L_avg / (rowsize * colsize));

	for (int col = 0; col < img_L.cols; col++)
	{
		for (int row = 0; row < img_L.rows; row++)
		{
			if ((img_L.at<float>(row, col) > Lmax *0.01 && img_L.at<float>(row, col) < Lmax*0.8))
			{
				if (img_L.at<float>(row, col) > maxValue)
				{
					maxValue = img_L.at<float>(row, col);
				}

				if (img_L.at<float>(row, col) < minValue)
				{
					if (img_L.at<float>(row, col) == 0)
					{
						minValue = 1;
					}
					else
					{
						minValue = img_L.at<float>(row, col);
					}
				}
			}
		}
	}

	std::cout <<"ttt:"<< minValue << "\t" << maxValue << std::endl;
	alpha_power = 2 * (2 * log2(L_avg) - log2(minValue) - log2(maxValue)) / (log2(maxValue) - log2(minValue));
	alpha = 0.1* pow(2, alpha_power);
	//alpha = 0.9;
	std::cout << L_avg << std::endl;

	L_W_pow = log2(maxValue) - log2(minValue) - log2(100) / 2;
	L_W = 1.5 * pow(2, L_W_pow);
	for (int col = 0; col < colsize; col++)
	{
		for (int row = 0; row < rowsize; row++)
		{

			//L_S.at<float>(row, col) = alpha * log2(img_L.at<float>(row, col)+1) / L_avg;
			L_S.at<float>(row, col) = alpha * img_L.at<float>(row, col) / L_avg;
			L_T.at<float>(row, col) = L_S.at<float>(row, col) * (1 + (L_S.at<float>(row, col) / pow(L_W, 2))) / (1 + L_S.at<float>(row, col));
		}
	}



	//////normal
	double LTminValue, LTmaxValue;
	//cv::Point LTminLoc, LTmaxLoc;
	cv::minMaxLoc(L_T, &LTminValue, &LTmaxValue, &LTminLoc, &LTmaxLoc);
	for (int col = 0; col < colsize; col++)
	{
		for (int row = 0; row < rowsize; row++)
		{
			L_N.at<float>(row, col) = (L_T.at<float>(row, col) - LTminValue) * 100 / (LTmaxValue - LTminValue);
		}
	}



	cv::minMaxLoc(L_N, &LTminValue, &LTmaxValue, &LTminLoc, &LTmaxLoc);

	return L_N;
		
}

////====================================================================
////========================GammaCorrection=============================
////====================================================================
void GammaCorrection(cv::Mat& srcImg, float gamma)
{
	for (int i = 0; i < srcImg.rows; i++) {
		for (int j = 0; j < srcImg.cols; j++) {
			for (int c = 0; c < 3; c++) {
				float fpixel = (float)(srcImg.at<cv::Vec3b>(i, j)[c]) / 255.0;
				srcImg.at<cv::Vec3b>(i, j)[c] = cv::saturate_cast<uchar>(pow(fpixel, gamma) * 255.0);
			}
		}
	}
}
////====================================================================
////========================Error_diffusion_v3==========================
////====================================================================
void Error_diffusion_v1(cv::Mat& image, int way, std::vector<std::vector<int>> Downgrade)
{
	int rowsize = image.rows;
	int colsize = image.cols;
	int level_size = Downgrade[0].size();
	cv::Mat edge_img;
	edgedetcion(image, edge_img);
	if (way == 5)
	{

	}
	else
	{
		for (int row = 0; row < rowsize; row++)
		{
			for (int col = 0; col < colsize; col++)
			{
				if (edge_img.at<uchar>(row, col) == 0)
				{
					for (int c = 0; c < 3; c++)
					{
						//int level = round(((int)(temp_str.at<cv::Vec3b>(row_i, col_j)[b])+1) / level_size) * level_size - 1;
						for (int level = 0; level < level_size - 1; level++)
						{
							if ((Downgrade[c][level] < image.at<cv::Vec3f>(row, col)[c]) && (image.at<cv::Vec3f>(row, col)[c] < Downgrade[c][level + 1]))
							{
								if ((image.at<cv::Vec3f>(row, col)[c] - Downgrade[c][level]) > (Downgrade[c][level + 1] - image.at<cv::Vec3f>(row, col)[c]))
								{
									double error = image.at<cv::Vec3f>(row, col)[c] - Downgrade[c][level + 1];
									switch (way)
									{
									case 1:
										Floyd_Steinberg_algorithm(image, row, col, error, c);
										break;
									case 2:
										Atkinson_algorithm(image, row, col, error, c);
										break;
									case 3:
										Jarvis_Judice_Ninke(image, row, col, error, c);
										break;
									case 4:
										stucki(image, row, col, error, c);
										break;
									default:
										break;
									}
									image.at<cv::Vec3f>(row, col)[c] = Downgrade[c][level + 1];
								}
								else
								{
									double error = image.at<cv::Vec3f>(row, col)[c] - Downgrade[c][level];
									switch (way)
									{
									case 1:
										Floyd_Steinberg_algorithm(image, row, col, error, c);
										break;
									case 2:
										Atkinson_algorithm(image, row, col, error, c);
										break;
									case 3:
										Jarvis_Judice_Ninke(image, row, col, error, c);
										break;
									case 4:
										stucki(image, row, col, error, c);
										break;
									default:
										break;
									}
									image.at<cv::Vec3f>(row, col)[c] = Downgrade[c][level];
								}
								break;
							}
						}
					}
				}
				else
				{
					for (int c = 0; c < 3; c++)
					{
						for (int level = 0; level < level_size - 1; level++)
						{
							if (Downgrade[c][level] < image.at<cv::Vec3f>(row, col)[c] && image.at<cv::Vec3f>(row, col)[c] < Downgrade[c][level + 1])
							{
								if ((image.at<cv::Vec3f>(row, col)[c] - Downgrade[c][level]) > (Downgrade[c][level + 1] - image.at<cv::Vec3f>(row, col)[c]))
								{
									image.at<cv::Vec3f>(row, col)[c] = Downgrade[c][level + 1];
								}
								else
								{
									image.at<cv::Vec3f>(row, col)[c] = Downgrade[c][level];
								}
								break;
							}
						}
					}
				}

			}
		}
	}

}

////====================================================================
////========================Error_diffusion_v2==========================
////====================================================================
void Error_diffusion_v2(cv::Mat& image, int way, std::vector<std::vector<int>> Downgrade)
{
	int rowsize = image.rows;
	int colsize = image.cols;
	int level_size = Downgrade[0].size();
	if (way == 5)
	{

	}
	else
	{
		for (int row = 0; row < rowsize; row++)
		{
			for (int col = 0; col < colsize; col++)
			{
				for (int c = 0; c < 3; c++)
				{
					for (int level = 0; level < level_size - 1; level++)
					{
						if (Downgrade[c][level] < (int)(image.at<cv::Vec3b>(row, col)[c]) && (int)(image.at<cv::Vec3b>(row, col)[c]) < Downgrade[c][level + 1])
						{
							if (((int)(image.at<cv::Vec3b>(row, col)[c]) - Downgrade[c][level]) > (Downgrade[c][level + 1] - ((int)(image.at<cv::Vec3b>(row, col)[c]))))
							{
								int error = ((int)(image.at<cv::Vec3b>(row, col)[c]) - Downgrade[c][level + 1]);
								switch (way)
								{
								case 1:
									Floyd_Steinberg_algorithm_v2(image, row, col, error, c);
									break;
								case 2:
									Atkinson_algorithm(image, row, col, error, c);
									break;
								case 3:
									Jarvis_Judice_Ninke(image, row, col, error, c);
									break;
								case 4:
									stucki(image, row, col, error, c);
									break;
								default:
									break;
								}
								image.at<cv::Vec3b>(row, col)[c] = Downgrade[c][level + 1];
							}
							else
							{
								int error = ((int)(image.at<cv::Vec3b>(row, col)[c]) - Downgrade[c][level]);
								switch (way)
								{
								case 1:
									Floyd_Steinberg_algorithm_v2(image, row, col, error, c);
									break;
								case 2:
									Atkinson_algorithm(image, row, col, error, c);
									break;
								case 3:
									Jarvis_Judice_Ninke(image, row, col, error, c);
									break;
								case 4:
									stucki(image, row, col, error, c);
									break;
								default:
									break;
								}
								image.at<cv::Vec3b>(row, col)[c] = Downgrade[c][level];
							}
							break;
						}
					}
				}
			}
		}
	}

}

////====================================================================
////========================Error_diffusion_v3==========================
////====================================================================
void Error_diffusion_v3(cv::Mat& image, int way, std::vector<std::vector<int>> Downgrade)
{
	int rowsize = image.rows;
	int colsize = image.cols;
	int level_size = Downgrade[0].size();
	cv::Mat edge_img;
	edgedetcion(image, edge_img);
	if (way == 5)
	{

	}
	else
	{
		for (int row = 0; row < rowsize; row++)
		{
			for (int col = 0; col < colsize; col++)
			{
				if (edge_img.at<uchar>(row, col) == 0)
				{
					for (int c = 0; c < 3; c++)
					{
						//int level = round(((int)(temp_str.at<cv::Vec3b>(row_i, col_j)[b])+1) / level_size) * level_size - 1;
						for (int level = 0; level < level_size - 1; level++)
						{
							if (Downgrade[c][level] < (int)(image.at<cv::Vec3b>(row, col)[c]) && (int)(image.at<cv::Vec3b>(row, col)[c]) < Downgrade[c][level + 1])
							{
								if (((int)(image.at<cv::Vec3b>(row, col)[c]) - Downgrade[c][level]) > (Downgrade[c][level + 1] - ((int)(image.at<cv::Vec3b>(row, col)[c]))))
								{
									int error = ((int)(image.at<cv::Vec3b>(row, col)[c]) - Downgrade[c][level + 1]);
									switch (way)
									{
									case 1:
										Floyd_Steinberg_algorithm_v2(image, row, col, error, c);
										break;
									case 2:
										Atkinson_algorithm(image, row, col, error, c);
										break;
									case 3:
										Jarvis_Judice_Ninke(image, row, col, error, c);
										break;
									case 4:
										stucki(image, row, col, error, c);
										break;
									default:
										break;
									}
									image.at<cv::Vec3b>(row, col)[c] = Downgrade[c][level + 1];
								}
								else
								{
									int error = ((int)(image.at<cv::Vec3b>(row, col)[c]) - Downgrade[c][level]);
									switch (way)
									{
									case 1:
										Floyd_Steinberg_algorithm_v2(image, row, col, error, c);
										break;
									case 2:
										Atkinson_algorithm(image, row, col, error, c);
										break;
									case 3:
										Jarvis_Judice_Ninke(image, row, col, error, c);
										break;
									case 4:
										stucki(image, row, col, error, c);
										break;
									default:
										break;
									}
									image.at<cv::Vec3b>(row, col)[c] = Downgrade[c][level];
								}
								break;
							}
						}
					}
				}
				else
				{
					for (int c = 0; c < 3; c++)
					{
						for (int level = 0; level < level_size - 1; level++)
						{
							if (Downgrade[c][level] < (int)(image.at<cv::Vec3b>(row, col)[c]) && (int)(image.at<cv::Vec3b>(row, col)[c]) < Downgrade[c][level + 1])
							{
								if (((int)(image.at<cv::Vec3b>(row, col)[c]) - Downgrade[c][level]) > (Downgrade[c][level + 1] - ((int)(image.at<cv::Vec3b>(row, col)[c]))))
								{
									image.at<cv::Vec3b>(row, col)[c] = Downgrade[c][level + 1];
								}
								else
								{
									image.at<cv::Vec3b>(row, col)[c] = Downgrade[c][level];
								}
								break;
							}
						}
					}
				}

			}
		}
	}

}

////====================================================================
////========================Error_diffusion_v4==========================
////====================================================================
void Error_diffusion_v4(cv::Mat& image, int way, std::vector<std::vector<int>> Downgrade)
{
	int rowsize = image.rows;
	int colsize = image.cols;
	int level_size = Downgrade[0].size();
	cv::Mat edge_img;
	edgedetcion(image, edge_img);
	if (way == 5)
	{

	}
	else
	{
		for (int row = 0; row < rowsize; row++)
		{
			for (int col = 0; col < colsize; col++)
			{
				if (edge_img.at<uchar>(row, col) == 0)
				{
					for (int c = 0; c < 3; c++)
					{
						int level = round(((int)(image.at<cv::Vec3b>(row, col)[c])) / (256.0 / level_size));
						//int error = ((int)(image.at<cv::Vec3b>(row, col)[c]) - Downgrade[c][level]);
						int error = (Downgrade[c][level] - (int)(image.at<cv::Vec3b>(row, col)[c]));
						switch (way)
						{
						case 1:
							Floyd_Steinberg_algorithm_v2(image, row, col, error, c);
							break;
						case 2:
							Atkinson_algorithm(image, row, col, error, c);
							break;
						case 3:
							Jarvis_Judice_Ninke(image, row, col, error, c);
							break;
						case 4:
							stucki(image, row, col, error, c);
							break;
						default:
							break;
						}
						image.at<cv::Vec3b>(row, col)[c] = Downgrade[c][level];
					}
				}
				else
				{
					for (int c = 0; c < 3; c++)
					{
						int level = round(((int)(image.at<cv::Vec3b>(row, col)[c])) / (256.0 / level_size));
						image.at<cv::Vec3b>(row, col)[c] = Downgrade[c][level];
					}
				}
			}
		}
	}

}

////====================================================================
////========================Floyd_Steinberg_algorithm_2=================
////====================================================================
void Floyd_Steinberg_algorithm(cv::Mat& image, int row, int col, double err, int c)
{
	double w1 = 7.0 / 16;
	double w2 = 3.0 / 16;
	double w3 = 5.0 / 16;
	double w4 = 1.0 / 16;
	int rowsize = image.rows;
	int colsize = image.cols;
	if (col + 1 < colsize)
		if (image.at<cv::Vec3f>(row, col + 1)[c] + w1 * err > 255)
			image.at<cv::Vec3f>(row, col + 1)[c] = 255;
		else if (image.at<cv::Vec3f>(row, col + 1)[c] + w1 * err < 0)
			image.at<cv::Vec3f>(row, col + 1)[c] = 0;
		else
			image.at<cv::Vec3f>(row, col + 1)[c] = image.at<cv::Vec3f>(row, col + 1)[c] + w1 * err;

	if (row + 1 < rowsize && col - 1 > 0)
		if (image.at<cv::Vec3f>(row + 1, col - 1)[c] + w2 * err > 255)
			image.at<cv::Vec3f>(row + 1, col - 1)[c] = 255;
		else if (image.at<cv::Vec3f>(row + 1, col - 1)[c] + w2 * err < 0)
			image.at<cv::Vec3f>(row + 1, col - 1)[c] = 0;
		else
			image.at<cv::Vec3f>(row + 1, col - 1)[c] = image.at<cv::Vec3f>(row + 1, col - 1)[c] + w2 * err;

	if (row + 1 < rowsize)
		if (image.at<cv::Vec3f>(row + 1, col)[c] + w3 * err > 255)
			image.at<cv::Vec3f>(row + 1, col)[c] = 255;
		else if (image.at<cv::Vec3f>(row + 1, col)[c] + w3 * err < 0)
			image.at<cv::Vec3f>(row + 1, col)[c] = 0;
		else
			image.at<cv::Vec3f>(row + 1, col)[c] = image.at<cv::Vec3f>(row + 1, col)[c] + w3 * err;

	if (row + 1 < rowsize && col + 1 < colsize)
		if (image.at<cv::Vec3f>(row + 1, col + 1)[c] + w4 * err > 255)
			image.at<cv::Vec3f>(row + 1, col + 1)[c] = 255;
		else if (image.at<cv::Vec3f>(row + 1, col + 1)[c] + w4 * err < 0)
			image.at<cv::Vec3f>(row + 1, col + 1)[c] = 0;
		else
			image.at<cv::Vec3f>(row + 1, col + 1)[c] = image.at<cv::Vec3f>(row + 1, col + 1)[c] + w4 * err;

}


void Floyd_Steinberg_algorithm_v2(cv::Mat& image, int row, int col, double err, int c)
{
	double w1 = 7.0 / 16;
	double w2 = 3.0 / 16;
	double w3 = 5.0 / 16;
	double w4 = 1.0 / 16;
	int rowsize = image.rows;
	int colsize = image.cols;
	if (col + 1 < colsize)
		if (image.at<cv::Vec3b>(row, col + 1)[c] + w1 * err > 255)
			image.at<cv::Vec3b>(row, col + 1)[c] = 255;
		else if (image.at<cv::Vec3b>(row, col + 1)[c] + w1 * err < 0)
			image.at<cv::Vec3b>(row, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row, col + 1)[c] = image.at<cv::Vec3b>(row, col + 1)[c] + w1 * err;

	if (row + 1 < rowsize && col - 1 > 0)
		if (image.at<cv::Vec3b>(row + 1, col - 1)[c] + w2 * err > 255)
			image.at<cv::Vec3b>(row + 1, col - 1)[c] = 255;
		else if (image.at<cv::Vec3b>(row + 1, col - 1)[c] + w2 * err < 0)
			image.at<cv::Vec3b>(row + 1, col - 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col - 1)[c] = image.at<cv::Vec3b>(row + 1, col - 1)[c] + w2 * err;

	if (row + 1 < rowsize)
		if (image.at<cv::Vec3b>(row + 1, col)[c] + w3 * err > 255)
			image.at<cv::Vec3b>(row + 1, col)[c] = 255;
		else if (image.at<cv::Vec3b>(row + 1, col)[c] + w3 * err < 0)
			image.at<cv::Vec3b>(row + 1, col)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col)[c] = image.at<cv::Vec3b>(row + 1, col)[c] + w3 * err;

	if (row + 1 < rowsize && col + 1 < colsize)
		if (image.at<cv::Vec3b>(row + 1, col + 1)[c] + w4 * err > 255)
			image.at<cv::Vec3b>(row + 1, col + 1)[c] = 255;
		else if (image.at<cv::Vec3b>(row + 1, col + 1)[c] + w4 * err < 0)
			image.at<cv::Vec3b>(row + 1, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col + 1)[c] = image.at<cv::Vec3b>(row + 1, col + 1)[c] + w4 * err;

}

void Atkinson_algorithm(cv::Mat& image, int row, int col, double err, int c)
{
	double w1 = 1.0 / 8;
	int rowsize = image.rows;
	int colsize = image.cols;
	double temp;
	if (col + 1 < colsize)
	{
		temp = image.at<cv::Vec3b>(row, col + 1)[c];
		if (temp + w1 * err > 255)
			image.at<cv::Vec3b>(row, col + 1)[c] = 255;
		else if (temp + w1 * err < 0)
			image.at<cv::Vec3b>(row, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row, col + 1)[c] += w1 * err;
	}

	if (col + 2 < colsize)
	{
		temp = image.at<cv::Vec3b>(row, col + 2)[c];
		if (temp + w1 * err > 255)
			image.at<cv::Vec3b>(row, col + 2)[c] = 255;
		else if (temp + w1 * err < 0)
			image.at<cv::Vec3b>(row, col + 2)[c] = 0;
		else
			image.at<cv::Vec3b>(row, col + 2)[c] += w1 * err;
	}

	if (0 < col - 1 && row + 1 < rowsize)
	{
		temp = image.at<cv::Vec3b>(row + 1, col - 1)[c];
		if (temp + w1 * err > 255)
			image.at<cv::Vec3b>(row + 1, col - 1)[c] = 255;
		else if (temp + w1 * err < 0)
			image.at<cv::Vec3b>(row + 1, col - 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col - 1)[c] += w1 * err;

	}

	if (row + 1 < rowsize)
	{
		temp = image.at<cv::Vec3b>(row + 1, col)[c];
		if (temp + w1 * err > 255)
			image.at<cv::Vec3b>(row + 1, col)[c] = 255;
		else if (temp + w1 * err < 0)
			image.at<cv::Vec3b>(row + 1, col)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col)[c] += w1 * err;

	}

	if (col + 1 < colsize && row + 1 < rowsize)
	{
		if (image.at<cv::Vec3b>(row + 1, col + 1)[c] + w1 * err > 255)
			image.at<cv::Vec3b>(row + 1, col + 1)[c] = 255;
		else if (image.at<cv::Vec3b>(row + 1, col + 1)[c] + w1 * err < 0)
			image.at<cv::Vec3b>(row + 1, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col + 1)[c] = image.at<cv::Vec3b>(row + 1, col + 1)[c] + w1 * err;
	}

	if (row + 2 < rowsize)
	{
		temp = image.at<cv::Vec3b>(row + 2, col)[c];
		if (temp + w1 * err > 255)
			image.at<cv::Vec3b>(row + 2, col)[c] = 255;
		else if (temp + w1 * err < 0)
			image.at<cv::Vec3b>(row + 2, col)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 2, col)[c] += w1 * err;
	}
}

void Jarvis_Judice_Ninke(cv::Mat& image, int row, int col, double err, int c)
{
	double w7 = 7.0 / 48;
	double w5 = 5.0 / 48;
	double w3 = 3.0 / 48;
	double w1 = 1.0 / 48;

	int rowsize = image.rows;
	int colsize = image.cols;
	double temp;

	if (row + 1 < rowsize)
	{
		temp = image.at<cv::Vec3b>(row + 1, col)[c] + w7 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 1, col)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 1, col)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col)[c] = temp;

	}
	if (row + 2 < rowsize)
	{
		temp = image.at<cv::Vec3b>(row + 2, col)[c] + w5 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 2, col)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 2, col)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 2, col)[c] = temp;
	}

	if (0 < row - 2 && col + 1 < colsize)
	{
		temp = image.at<cv::Vec3b>(row - 2, col + 1)[c] + w3 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row - 2, col + 1)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row - 2, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row - 2, col + 1)[c] = temp;

	}
	if (0 < row - 1 && col + 1 < colsize)
	{
		temp = image.at<cv::Vec3b>(row - 1, col + 1)[c] + w5 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row - 1, col + 1)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row - 1, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row - 1, col + 1)[c] = temp;
	}
	if (col + 1 < colsize)
	{
		temp = image.at<cv::Vec3b>(row, col + 1)[c] + w7 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row, col + 1)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row, col + 1)[c] = temp;
	}
	if (row + 1 < rowsize && col + 1 < colsize)
	{
		temp = image.at<cv::Vec3b>(row + 1, col + 1)[c] + w5 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 1, col + 1)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 1, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col + 1)[c] = temp;
	}
	if (row + 2 < rowsize && col + 1 < colsize)
	{
		temp = image.at<cv::Vec3b>(row + 2, col + 1)[c] + w3 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 2, col + 1)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 2, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 2, col + 1)[c] = temp;
	}
	if (row - 2 > 0 && col + 2 < colsize)
	{
		temp = image.at<cv::Vec3b>(row - 2, col + 2)[c] + w1 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row - 2, col + 2)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row - 2, col + 2)[c] = 0;
		else
			image.at<cv::Vec3b>(row - 2, col + 2)[c] = temp;
	}
	if (row - 1 > 0 && col + 2 < colsize)
	{
		temp = image.at<cv::Vec3b>(row - 1, col + 2)[c] + w3 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row - 1, col + 2)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row - 1, col + 2)[c] = 0;
		else
			image.at<cv::Vec3b>(row - 1, col + 2)[c] = temp;
	}
	if (col + 2 < colsize)
	{
		temp = image.at<cv::Vec3b>(row, col + 2)[c] + w5 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row, col + 2)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row, col + 2)[c] = 0;
		else
			image.at<cv::Vec3b>(row, col + 2)[c] = temp;
	}
	if (row + 1 < rowsize && col + 2 < colsize)
	{
		temp = image.at<cv::Vec3b>(row + 1, col + 2)[c] + w3 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 1, col + 2)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 1, col + 2)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col + 2)[c] = temp;
	}
	if (row + 2 < rowsize && col + 2 < colsize)
	{
		temp = image.at<cv::Vec3b>(row + 2, col + 2)[c] + w1 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 2, col + 2)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 2, col + 2)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 2, col + 2)[c] = temp;
	}


}

void stucki(cv::Mat& image, int row, int col, double err, int c)
{
	double w8 = 8.0 / 42;
	double w7 = 7.0 / 42;
	double w5 = 5.0 / 42;
	double w4 = 4.0 / 42;
	double w2 = 2.0 / 42;
	double w1 = 1.0 / 42;

	int rowsize = image.rows;
	int colsize = image.cols;
	double temp;
	if (row + 1 < rowsize)
	{
		temp = image.at<cv::Vec3b>(row + 1, col)[c] + w7 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 1, col)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 1, col)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col)[c] = temp;
	}
	if (row + 2 < rowsize)
	{
		temp = image.at<cv::Vec3b>(row + 2, col)[c] + w5 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 2, col)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 2, col)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 2, col)[c] = temp;
	}
	if (0 < row - 2 && col + 1 < colsize)
	{
		temp = image.at<cv::Vec3b>(row - 2, col + 1)[c] + w2 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row - 2, col + 1)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row - 2, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row - 2, col + 1)[c] = temp;

	}
	if (0 < row - 1 && col + 1 < colsize)
	{
		temp = image.at<cv::Vec3b>(row - 1, col + 1)[c] + w4 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row - 1, col + 1)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row - 1, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row - 1, col + 1)[c] = temp;
	}
	if (col + 1 < colsize)
	{
		temp = image.at<cv::Vec3b>(row, col + 1)[c] + w8 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row, col + 1)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row, col + 1)[c] = temp;
	}
	if (row + 1 < rowsize && col + 1 < colsize)
	{
		temp = image.at<cv::Vec3b>(row + 1, col + 1)[c] + w4 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 1, col + 1)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 1, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col + 1)[c] = temp;
	}
	if (row + 2 < rowsize && col + 1 < colsize)
	{
		temp = image.at<cv::Vec3b>(row + 2, col + 1)[c] + w2 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 2, col + 1)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 2, col + 1)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 2, col + 1)[c] = temp;
	}
	if (row - 2 > 0 && col + 2 < colsize)
	{
		temp = image.at<cv::Vec3b>(row - 2, col + 2)[c] + w1 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row - 2, col + 2)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row - 2, col + 2)[c] = 0;
		else
			image.at<cv::Vec3b>(row - 2, col + 2)[c] = temp;
	}
	if (row + 1 < rowsize && col + 2 < colsize)
	{
		temp = image.at<cv::Vec3b>(row + 1, col + 2)[c] + w2 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 1, col + 2)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 1, col + 2)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col + 2)[c] = temp;
	}
	if (col + 2 < colsize)
	{
		temp = image.at<cv::Vec3b>(row, col + 2)[c] + w4 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row, col + 2)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row, col + 2)[c] = 0;
		else
			image.at<cv::Vec3b>(row, col + 2)[c] = temp;
	}
	if (row + 1 < rowsize && col + 2 < colsize)
	{
		temp = image.at<cv::Vec3b>(row + 1, col + 2)[c] + w2 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 1, col + 2)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 1, col + 2)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 1, col + 2)[c] = temp;
	}
	if (row + 2 < rowsize && col + 2 < colsize)
	{
		temp = image.at<cv::Vec3b>(row + 2, col + 2)[c] + w1 * err;
		if (temp > 255)
			image.at<cv::Vec3b>(row + 2, col + 2)[c] = 255;
		else if (temp < 0)
			image.at<cv::Vec3b>(row + 2, col + 2)[c] = 0;
		else
			image.at<cv::Vec3b>(row + 2, col + 2)[c] = temp;
	}
}

////====================================================================
////========================colorcorrection=============================
////====================================================================
cv::Mat colorcorrection_one(cv::Mat img, std::vector<double>matrix)
{
	cv::Mat dst = img.clone();
	for (int r = 0; r < img.rows; r++)
	{
		for (int c = 0; c < img.cols; c++)
		{
			for (int channel = 0; channel < 3; channel++)
			{
				int m1 = 0 + 3 * channel, m2 = 1 + 3 * channel, m3 = 2 + 3 * channel;
				double temp = matrix[m1] * img.at<cv::Vec3f>(r, c)[0] +
					matrix[m2] * img.at<cv::Vec3f>(r, c)[1] +
					matrix[m3] * img.at<cv::Vec3f>(r, c)[2];

				dst.at<cv::Vec3f>(r, c)[channel] = temp;


			}
		}
	}
	//cv::Mat out = oversaturationset(dst);

	std::vector<double> img_sort;
	for (int y = 0; y < dst.rows; y++)
	{
		for (int x = 0; x < dst.cols; x++)
		{
			for (int channel = 0; channel < 3; channel++)
			{
				img_sort.push_back(dst.at<cv::Vec3f>(y, x)[channel]);
				//std::cout << img.at<cv::Vec3f>(y, x)[channel] <<std::endl;
			}
		}
	}
	std::sort(img_sort.begin(), img_sort.end(), std::greater<double>());
	double LTminValue, LTmaxValue;
	LTmaxValue = img_sort[0];
	LTminValue = img_sort[img_sort.size()-1];
	cv::Mat out = dst.clone();
	//////normal
	std::vector<cv::Mat>ch;
	cv::split(dst, ch);
	for (int i = 0; i < 3; i++)
	{
		for (int col = 0; col < dst.cols; col++)
		{
			for (int row = 0; row < dst.rows; row++)
			{

				out.at<cv::Vec3f>(row, col)[i] = (dst.at<cv::Vec3f>(row, col)[i] - LTminValue) * 255 / (LTmaxValue - LTminValue);
			}

		}
	}


	return out;

}

////====================================================================
////========================edgedetcion=================================
////====================================================================
//void edgedetcion(cv::Mat img, cv::Mat& edge_img)
//{
//	cv::Mat gray;
//	cv::Mat src;
//	int scale = 7;
//	int delta = 0;
//	int ddepth = CV_16S;
//
//	cv::GaussianBlur(img, src, cv::Size(3, 3), 0, 0, cv::BORDER_DEFAULT);
//	cv::cvtColor(src, gray, cv::COLOR_BGR2GRAY);
//	cv::Mat grad, grad_x, grad_y;
//	cv::Mat abs_grad_x, abs_grad_y;
//	Sobel(gray, grad_x, ddepth, 1, 0, 3, scale, delta, cv::BORDER_DEFAULT);
//	convertScaleAbs(grad_x, abs_grad_x);
//	Sobel(gray, grad_y, ddepth, 0, 1, 3, scale, delta, cv::BORDER_DEFAULT);
//	convertScaleAbs(grad_y, abs_grad_y);
//	addWeighted(abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad);
//	//cv::threshold(grad, grad, 127, 255, 0);
//	cv::imwrite("grad.bmp", grad);
//	edge_img = grad.clone();
//}

////====================================================================
////========================macro_edge==================================
////====================================================================
//void macro(cv::Mat img, cv::Mat& out_img)
//{
//	cv::Mat edge_img;
//	cv::Mat	macro_img(img.rows, img.cols, CV_8UC1);
//	edgedetcion(img, edge_img);
//	//cv::imshow("edge", edge_img);
//	for (int col = 0; col < macro_img.cols; col++)
//	{
//		for (int row = 0; row < macro_img.rows; row++)
//		{
//			if ((int)edge_img.at<uchar>(row, col) == 255)
//			{
//				for (int i = -4; i < 4; i++)
//				{
//					for (int j = -4; j < 4; j++)
//					{
//						if (row + i < macro_img.rows && col + j < macro_img.cols && row + i > 0 && col + j > 0)
//						{
//							macro_img.at<uchar>(row + i, col + j) = 255;
//						}
//					}
//				}
//			}
//		}
//	}
//
//	macro_img.copyTo(out_img);
//}
//
//void macroEqualizationhis(cv::Mat img, cv::Mat& out_img)
//{
//	cv::Mat macro_img;
//	macro(img, macro_img);
//	std::vector<double> val(256, 0);
//	std::vector<cv::Mat> channel_Lab(3);
//	cv::Mat img_Lab;
//
//
//	cv::cvtColor(img, img_Lab, cv::COLOR_BGR2Lab);
//	cv::split(img_Lab, channel_Lab);
//	double sum = 0;
//	for (int x = 0; x < channel_Lab[0].cols; x++)
//	{
//		for (int y = 0; y < channel_Lab[0].rows; y++)
//		{
//			if ((int)macro_img.at<uchar>(y, x) == 255)
//			{
//				int temp = (int)channel_Lab[0].at<uchar>(y, x);
//				val[temp] += 1;
//				sum += 1;
//			}
//		}
//	}
//	std::vector<double>PDF(256);
//	for (int P = 0; P < val.size(); P++)
//	{
//		PDF[P] = val[P] / sum;
//	}
//	std::vector<double>CDF(256);
//	CDF[0] = PDF[0] * 255.0;
//
//	for (int C = 1; C < 256; C++)
//	{
//		CDF[C] = CDF[C - 1] + PDF[C] * 255.0;
//	}
//	for (int C = 0; C < 256; C++)
//	{
//		CDF[C] = round(CDF[C]);
//	}
//	for (int y = 0; y < channel_Lab[0].rows; y++)
//	{
//		for (int x = 0; x < channel_Lab[0].cols; x++)
//		{
//			int temp = (int)channel_Lab[0].at<uchar>(y, x);
//			channel_Lab[0].at<uchar>(y, x) = CDF[temp];
//			//if (channel_Lab[1].at<uchar>(y,x)>127) 
//			//{
//			//	channel_Lab[1].at<uchar>(y, x) = channel_Lab[1].at<uchar>(y, x) * 1.1;
//			//}
//			//else
//			//{
//			//	channel_Lab[1].at<uchar>(y, x) = channel_Lab[1].at<uchar>(y, x) / 1.1;
//			//}
//			//if (channel_Lab[2].at<uchar>(y, x) > 127)
//			//{
//			//	channel_Lab[2].at<uchar>(y, x) = channel_Lab[2].at<uchar>(y, x) * 1.1;
//			//}
//			//else
//			//{
//			//	channel_Lab[2].at<uchar>(y, x) = channel_Lab[2].at<uchar>(y, x) / 1.1;
//			//}
//		}
//	}
//	cv::merge(channel_Lab, img_Lab);
//	cv::cvtColor(img_Lab, img_Lab, cv::COLOR_Lab2BGR);
//	img_Lab.convertTo(img_Lab, CV_32FC3, 1 / 255.0);
//	cv::cvtColor(img_Lab, img_Lab, cv::COLOR_BGR2Lab);
//	cv::split(img_Lab, channel_Lab);
//	channel_Lab[1] = channel_Lab[1] * 1.5;
//	channel_Lab[2] = channel_Lab[2] * 1.5;
//	cv::merge(channel_Lab, img_Lab);
//	cv::cvtColor(img_Lab, img_Lab, cv::COLOR_Lab2BGR);
//	img_Lab.convertTo(img_Lab, CV_32FC3, 255.0);
//	img_Lab.convertTo(img_Lab, CV_8UC3);
//	out_img = img_Lab;
//}

void colorcorrection(cv::Mat img, cv::Mat& dst, std::vector<std::vector<double>>matirx, std::vector<std::vector<double>> center)
{
	cv::Mat img_Lab(img.rows, img.cols, CV_32FC3);
	for (int r = 0; r < img.rows; r++)
	{
		for (int c = 0; c < img.cols; c++)
		{
			img_Lab.at<cv::Vec3f>(r, c) = img.at<cv::Vec3b>(r, c) / 255.0;
		}
	}
	cv::cvtColor(img_Lab, img_Lab, cv::COLOR_RGB2Lab);
	for (int r = 0; r < img.rows; r++)
	{
		for (int c = 0; c < img.cols; c++)
		{
			std::vector<double> temp(3);
			std::vector<double> reg;
			for (int channel = 0; channel < 3; channel++)
			{
				temp[channel] = img_Lab.at<cv::Vec3f>(r, c)[channel];
			}
			for (int k = 0; k < center.size(); k++)
			{
				reg.push_back(norm2(temp, center[k]));
			}
			auto min = std::min_element(reg.begin(), reg.end());
			int min_idx = std::distance(reg.begin(), min);
			for (int channel = 0; channel < 3; channel++)
			{
				double temp2 = matirx[min_idx][0 + 3 * channel] * img.at<cv::Vec3b>(r, c)[0] +
					matirx[min_idx][1 + 3 * channel] * img.at<cv::Vec3b>(r, c)[1] +
					matirx[min_idx][2 + 3 * channel] * img.at<cv::Vec3b>(r, c)[2];
				if (temp2 > 255)
				{
					dst.at<cv::Vec3b>(r, c)[channel] = 255;
				}
				else if (temp2 < 0)
				{
					dst.at<cv::Vec3b>(r, c)[channel] = 0;
				}
				else
				{
					dst.at<cv::Vec3b>(r, c)[channel] = temp2;
				}

			}
		}
	}
}


//void macroEqualizationhis_Yuv(cv::Mat img, cv::Mat& out_img)
//{
//	cv::Mat macro_img;
//	macro(img, macro_img);
//	std::vector<double> val(256, 0);
//	std::vector<cv::Mat> channel_Lab(3);
//	cv::Mat img_Lab;
//
//
//	cv::cvtColor(img, img_Lab, cv::COLOR_BGR2YUV);
//	cv::split(img_Lab, channel_Lab);
//	double sum = 0;
//	for (int x = 0; x < channel_Lab[0].cols; x++)
//	{
//		for (int y = 0; y < channel_Lab[0].rows; y++)
//		{
//			if ((int)macro_img.at<uchar>(y, x) == 255)
//			{
//				int temp = (int)channel_Lab[0].at<uchar>(y, x);
//				val[temp] += 1;
//				sum += 1;
//			}
//		}
//	}
//	std::vector<double>PDF(256);
//	for (int P = 0; P < val.size(); P++)
//	{
//		PDF[P] = val[P] / sum;
//	}
//	std::vector<double>CDF(256);
//	CDF[0] = PDF[0] * 255.0;
//
//	for (int C = 1; C < 256; C++)
//	{
//		CDF[C] = CDF[C - 1] + PDF[C] * 255.0;
//	}
//	for (int C = 0; C < 256; C++)
//	{
//		CDF[C] = round(CDF[C]);
//	}
//	for (int y = 0; y < channel_Lab[0].rows; y++)
//	{
//		for (int x = 0; x < channel_Lab[0].cols; x++)
//		{
//			int temp = (int)channel_Lab[0].at<uchar>(y, x);
//			channel_Lab[0].at<uchar>(y, x) = CDF[temp];
//
//			if (channel_Lab[1].at<uchar>(y, x) > 128)
//			{
//
//				if (channel_Lab[1].at<uchar>(y, x) * 1.05 > 255)
//				{
//					channel_Lab[1].at<uchar>(y, x) = channel_Lab[1].at<uchar>(y, x);
//				}
//				else
//				{
//					channel_Lab[1].at<uchar>(y, x) = channel_Lab[1].at<uchar>(y, x) * 1.05;
//				}
//			}
//			else
//			{
//				channel_Lab[1].at<uchar>(y, x) = channel_Lab[1].at<uchar>(y, x) / 1.05;
//			}
//
//			if (channel_Lab[2].at<uchar>(y, x) > 128)
//			{
//				if (channel_Lab[2].at<uchar>(y, x) * 1.05 > 255)
//				{
//					channel_Lab[2].at<uchar>(y, x) = channel_Lab[2].at<uchar>(y, x);
//				}
//				else
//				{
//					channel_Lab[2].at<uchar>(y, x) = channel_Lab[2].at<uchar>(y, x) * 1.05;
//				}
//			}
//			else
//			{
//				channel_Lab[2].at<uchar>(y, x) = channel_Lab[2].at<uchar>(y, x) / 1.05;
//			}
//		}
//	}
//
//
//
//	cv::merge(channel_Lab, img_Lab);
//	cv::cvtColor(img_Lab, img_Lab, cv::COLOR_YUV2BGR);
//
//	out_img = img_Lab;
//}