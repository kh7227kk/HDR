#include "statistics.h"
#include "mydefine.h"
void imgprocessing(Filedata img, cv::Mat& out_img, int level_size, int ed_mode);
void imgprocessing_v2(cv::Mat img, cv::Mat& out_img, int level_size, int ed_mode);
void imgprocessing_v3(cv::Mat img, cv::Mat& out_img, int level_size, int ed_mode);

void GammaCorrection(cv::Mat& srcImg, float gamma);

cv::Mat balance_white(cv::Mat mat);

cv::Mat reproduction_Lab(cv::Mat img, Filedata org_img);
cv::Mat tonemapping_Lab(cv::Mat img_L,cv::Mat str);
void Floyd_Steinberg_algorithm(cv::Mat& image, int row, int col, double err, int c);
void Floyd_Steinberg_algorithm_v2(cv::Mat& image, int row, int col, double err, int c);
void Atkinson_algorithm(cv::Mat& image, int row, int col, double err, int c);
void Jarvis_Judice_Ninke(cv::Mat& image, int row, int col, double err, int c);
void stucki(cv::Mat& image, int row, int col, double err, int c);
void Error_diffusion_v1(cv::Mat& image, int way, std::vector<std::vector<int>> Downgrade);
void Error_diffusion_v2(cv::Mat& image, int way, std::vector<std::vector<int>> Downgrade);

void Error_diffusion_v3(cv::Mat& image, int way, std::vector<std::vector<int>> Downgrade);

void Error_diffusion_v4(cv::Mat& image, int way, std::vector<std::vector<int>> Downgrade);

void edgedetcion(cv::Mat img, cv::Mat& edge_img);


//void colorcorrection(cv::Mat img, cv::Mat& dst, std::vector<double>matrix);
//void macro(cv::Mat img, cv::Mat& out_img);
//void macroEqualizationhis(cv::Mat img, cv::Mat& out_img);
cv::Mat colorcorrection_one(cv::Mat img, std::vector<double>matrix);
//void colorcorrection(cv::Mat img, cv::Mat& dst, std::vector<std::vector<double>>matirx, std::vector<std::vector<double>> center);
//
//void macroEqualizationhis_Yuv(cv::Mat img, cv::Mat& out_img);
//
//
