#include "GeodeticCalculator.h"
 
/**
* @brief 计算高斯平面坐标的前向转换。
* @param gauss 输入的高斯平面坐标，包含：
* 			- B：纬度（度）
* 			- L：经度（度）
* 			- width：带宽（可选3、6；输入0则需要补全L0）
* 			- sigmaB：纬度标准差（秒）
* 			- sigmaL：经度标准差（秒）
*			- NF：北向偏移量（米）
* 			以上参数根据需要进行补全
*			以下参数在函数执行后进行赋值输出
* 			- x：输出的高斯平面坐标X（米）
* 			- y：输出的高斯平面坐标Y（米）
* 			- sigmaX：X坐标标准差（米）
* 			- sigmaY：Y坐标标准差（米）
* 			- mca：子午线收敛角（度）
* 			- sigmaMCA：子午线收敛角标准差（秒）
* 			- signeded：带号
* 			- L0：中央子午线（度）
*/
void GeodeticCalculator::GaussForward(PointGauss& gauss, Ellipsoid::Ellipsoid_para ell)
{
	if (gauss.width == 0)
	{
		caclGaussCoordForward(gauss,ell);
		if (gauss.y < 500000)
		{
			gauss.y += 500000;
			gauss.NF = 500000;
		}
		PointMCA mca(gauss);
		caclMeridianConvergenceAngle(mca, ell);
		gauss.mca = mca.mca;
		gauss.sigmaMCA = mca.sigmaMCA;
	}
	else if (gauss.width == 3|| gauss.width==6)
	{
		gauss.signeded = caclSignede(gauss.L, gauss.width);
		gauss.L0 = caclCentralMeridian(gauss.signeded, gauss.width);
		caclGaussCoordForward(gauss,ell);
		if (gauss.y < 500000)
		{
			gauss.y += 500000;
			gauss.NF = 500000;
		}
		PointMCA mca(gauss);
		caclMeridianConvergenceAngle(mca, ell);
		gauss.mca = mca.mca;
		gauss.sigmaMCA = mca.sigmaMCA;
	}
}

/**
* @brief 计算高斯平面坐标的反向转换。
* @param gauss 输入的高斯平面坐标，包含：
* 			- x：高斯平面坐标X（米）
* 			- y：高斯平面坐标Y（米）
* 			- width：带宽（可选3、6；输入0则需要补全L0）
* 			- sigmaX：X坐标标准差（米）
* 			- sigmaY：Y坐标标准差（米）
* 			- NF：北向偏移量（米）
* 			以上参数根据需要进行补全
* 			以下参数在函数执行后进行赋值输出
* 			- B：输出的纬度（度）
* 			- L：输出的经度（度）
* 			- sigmaB：纬度标准差（秒）
* 			- sigmaL：经度标准差（秒）
* 			- mca：子午线收敛角（度）
* 			- sigmaMCA：子午线收敛角标准差（秒）
* 			- signeded：带号
* 			- L0：中央子午线（度）
*/
void GeodeticCalculator::GaussBack(PointGauss& gauss, Ellipsoid::Ellipsoid_para ell)
{
	if(gauss.signeded == 0)
		gauss.signeded = gauss.y / 1000000;
	gauss.y = fmod(gauss.y, 1000000);
	if (gauss.width == 0)
	{
		caclGaussCoordBack(gauss, ell);
		PointMCA mca(gauss);
		caclMeridianConvergenceAngle(mca, ell);
		gauss.mca = mca.mca;
		gauss.sigmaMCA = mca.sigmaMCA;
	}
	else if (gauss.width == 3|| gauss.width == 6)
	{
		gauss.L0 = caclCentralMeridian(gauss.signeded, gauss.width);
		caclGaussCoordBack(gauss, ell);
		PointMCA mca(gauss);
		caclMeridianConvergenceAngle(mca, ell);
		gauss.mca = mca.mca;
		gauss.sigmaMCA = mca.sigmaMCA;
	}
}

void GeodeticCalculator::caclGaussCoordForward(PointGauss& gauss, Ellipsoid::Ellipsoid_para ell)
{
	double B = gauss.B, L = gauss.L, L0 = gauss.L0;
	double a = ell.a, f = ell.f;
	double dl = L - L0;
	B = B * pi / 180;
	L = L * pi / 180;
	dl = dl * pi / 180;
	double b = a * (1 - f);
	double e1 = sqrt(a * a - b * b) / a;
	double e2 = sqrt(a * a - b * b) / b;
	double N = a / sqrt(1 - e1 * e1 * sin(B) * sin(B));
	double t = tan(B);
	double n = e2 * cos(B);
	double m0 = a * (1 - e1 * e1);
	double m2 = 3 * e1 * e1 * m0 / 2;
	double m4 = 5 * e1 * e1 * m2 / 4;
	double m6 = 7 * e1 * e1 * m4 / 6;
	double m8 = 9 * e1 * e1 * m6 / 8;
	double a0 = m0 + m2 / 2 + 3 * m4 / 8 + 5 * m6 / 16 + 35 * m8 / 128;
	double a2 = m2 / 2 + m4 / 2 + 15 * m6 / 32 + 7 * m8 / 16;
	double a4 = m4 / 8 + 3 * m6 / 16 + 7 * m8 / 32;
	double a6 = m6 / 32 + m8 / 16;
	double a8 = m8 / 128;
	double X = a0 * B - a2 * sin(2 * B) / 2 + a4 * sin(4 * B) / 4 - a6 * sin(6 * B) / 6 + a8 * sin(8 * B) / 8;
	double x = X +
		t * N * pow(cos(B), 2) * dl * dl / 2 +
		t * N * pow(cos(B), 4) * (5 - t * t + 9 * n * n + 4 * pow(n, 4)) * pow(dl, 4) / 24 +
		t * N * pow(cos(B), 6) * (61 - 58 * t * t + pow(t, 4) + 270 * n * n - 330 * t * t * n * n) * pow(dl, 6) / 720+
		t * N * pow(cos(B), 8) * (1385 - 3111 * t * t + 543 * pow(t, 4) - pow(t, 6)) * pow(dl, 8) / 40320;
	double y = N * cos(B) * dl 
		+ N * pow(cos(B), 3) / (6) * (1 - pow(t, 2) + pow(n, 2)) * pow(dl, 3)
		+ N * pow(cos(B), 5) / (120) * (5 - 18 * pow(t, 2) + pow(t, 4) + 14 * pow(n, 2) - 58 * pow(t, 2) * pow(n, 2)) * pow(dl, 5)
		+ N * pow(cos(B), 7) * (61 - 479 * t * t + 179 * pow(t, 4) - pow(t, 6)) * pow(dl, 7) / 5040;

	gauss.x = x;
	gauss.y = y;

	L0 = gauss.L0 * pi / 180;
	double dx_dB = a * (f - 1) * (f - 1) + (a * (L - L0) * (L - L0)) / (2 * sqrt((pow(f, 2) * pow(sin(B), 2) - 2 * f * sin(B) * sin(B) + 1))) - (a * pow((L - L0), 8) * (7266 * sin(B) * sin(B) - 10920 * pow(sin(B), 4) + 5040 * pow(sin(B), 6) - 1385)) / (40320 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) - (3 * a * f * (f - 1) * (f - 1) * (f - 2)) / 4 - (a * sin(B) * sin(B) * (L - L0) * (L - L0)) / sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1)) + (45 * a * f * f * (f - 1) * (f - 1) * (f - 2) * (f - 2)) / 64 - (175 * a * f * f * f * (f - 1) * (f - 1) * (f - 2) * (f - 2) * (f - 2)) / 256 + (11025 * a * f * f * f * f * (f - 1) * (f - 1) * pow((f - 2), 4)) / 16384 + (a * sin(B) * sin(B) * pow((L - L0), 8) * (7266 * sin(B) * sin(B) - 10920 * pow(sin(B), 4) + 5040 * pow(sin(B), 6) - 1385)) / (5040 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) - (3 * a * f * cos(2 * B) * (f - 1) * (f - 1) * (735 * pow(f, 7) - 5880 * pow(f, 7) + 16940 * pow(f, 5) - 19320 * pow(f, 4) + 4000 * pow(f, 3) + 3040 * f * f + 2048 * f + 1024)) / 2048 + (15 * a * f * f * cos(4 * B) * (f * f - 3 * f + 2) * (f * f - 3 * f + 2) * (147 * pow(f, 4) - 588 * pow(f, 3) + 476 * f * f + 224 * f + 64)) / 4096 - (a * sin(B) * sin(B) * pow((L - L0), 8) * (1400 * pow(sin(B), 4) - 2436 * sin(B) * sin(B) + 1037)) / (6720 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) + (315 * a * pow(f, 4) * cos(8 * B) * (f - 1) * (f - 1) * pow((f - 2), 4)) / 16384 + (35 * a * pow(f, 3) * cos(6 * B) * (f - 1) * (f - 1) * pow((f - 2), 3) * (-9 * f * f + 18 * f + 4)) / 2048 + (a * pow(cos(B), 6) * pow((L - L0), 6) * (tan(B) * tan(B) + 1) * (pow(tan(B), 4) - 58 * tan(B) * tan(B) - (270 * f * cos(B) * cos(B) * (f - 2)) / ((f - 1) * (f - 1)) + (330 * f * sin(B) * sin(B) * (f - 2)) / ((f - 1) * (f - 1)) + 61)) / (720 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) + (a * f * (L - L0) * (L - L0) * (cos(4 * B) / 8 - 1 / 8) * (f - 2)) / (2 * pow((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1), (3 / 2))) - (a * pow(cos(B), 4) * pow((L - L0), 4) * (tan(B) * tan(B) + 1) * (tan(B) * tan(B) - (4 * f * f * pow(cos(B), 4) * (f - 2) * (f - 2)) / pow((f - 1), 4) + (9 * f * cos(B) * cos(B) * (f - 2)) / pow((f - 1), 2) - 5)) / (24 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) + (a * pow(cos(B), 6) * tan(B) * pow((L - L0), 6) * (4 * pow(tan(B), 3) * (tan(B) * tan(B) + 1) - 116 * tan(B) * (tan(B) * tan(B) + 1) + (270 * f * sin(2 * B) * (f - 2)) / ((f - 1) * (f - 1)) + (660 * f * sin(B) * (f - 2)) / (cos(B) * (f - 1) * (f - 1)) - (660 * f * pow(sin(B), 3) * (f - 2)) / (cos(B) * (f - 1) * (f - 1)))) / (720 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) - (a * pow(cos(B), 4) * tan(B) * pow((L - L0), 4) * (2 * tan(B) * (tan(B) * tan(B) + 1) - (9 * f * sin(2 * B) * (f - 2)) / ((f - 1) * (f - 1)) + (16 * f * f * pow(cos(B), 3) * sin(B) * (f - 2) * (f - 2)) / pow((f - 1), 4))) / (24 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) - (a * pow(cos(B), 5) * sin(B) * tan(B) * pow((L - L0), 6) * (pow(tan(B), 4) - 58 * tan(B) * tan(B) - (270 * f * cos(B) * cos(B) * (f - 2)) / ((f - 1) * (f - 1)) + (330 * f * sin(B) * sin(B) * (f - 2)) / ((f - 1) * (f - 1)) + 61)) / (120 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) + (a * pow(cos(B), 3) * sin(B) * tan(B) * pow((L - L0), 4) * (tan(B) * tan(B) - (4 * f * f * pow(cos(B), 4) * (f - 2) * (f - 2)) / pow((f - 1), 4) + (9 * f * cos(B) * cos(B) * (f - 2)) / ((f - 1) * (f - 1)) - 5)) / (6 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) - (a * f * sin(B) * sin(B) * pow((L - L0), 8) * (f - 2) * (18186 * pow(sin(B), 4) - 8651 * sin(B) * sin(B) - 15960 * pow(sin(B), 6) + 5040 * pow(sin(B), 8) + 1385)) / (40320 * pow((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1), (3 / 2))) - (a * f * pow(cos(B), 7) * sin(B) * tan(B) * pow((L - L0), 6) * (f - 2) * (pow(tan(B), 4) - 58 * tan(B) * tan(B) - (270 * f * cos(B) * cos(B) * (f - 2)) / ((f - 1) * (f - 1)) + (330 * f * sin(B) * sin(B) * (f - 2)) / ((f - 1) * (f - 1)) + 61)) / (720 * pow((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1), (3 / 2))) + (a * f * pow(cos(B), 5) * sin(B) * tan(B) * pow((L - L0), 4) * (f - 2) * (tan(B) * tan(B) - (4 * f * f * pow(cos(B), 4) * pow((f - 2), 2)) / pow((f - 1), 4) + (9 * f * cos(B) * cos(B) * (f - 2)) / ((f - 1) * (f - 1)) - 5)) / (24 * pow((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1), (3 / 2)));
	double dx_dL = (a * cos(B) * sin(B) * (2 * L - 2 * L0)) / (2 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) - (a * pow(cos(B), 8) * tan(B) * pow((L - L0), 7) * (3111 * tan(B) * tan(B) - 543 * pow(tan(B), 4) + pow(tan(B), 6) - 1385)) / (5040 * sqrt((f * f * pow(sin(B), 2) - 2 * f * sin(B) * sin(B) + 1))) + (a * pow(cos(B), 6) * tan(B) * pow((L - L0), 5) * (pow(tan(B), 4) - 58 * pow(tan(B), 2) - (270 * f * pow(cos(B), 2) * (f - 2)) / pow((f - 1), 2) + (330 * f * sin(B) * sin(B) * (f - 2)) / ((f - 1) * (f - 1)) + 61)) / (120 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) - (a * pow(cos(B), 4) * tan(B) * pow((L - L0), 3) * (tan(B) * tan(B) - (4 * f * f * pow(cos(B), 4) * pow((f - 2), 2)) / pow((f - 1), 4) + (9 * f * pow(cos(B), 2) * (f - 2)) / pow((f - 1), 2) - 5)) / (6 * sqrt((pow(f, 2) * pow(sin(B), 2) - 2 * f * pow(sin(B), 2) + 1)));
	double dy_dB = (a * pow(cos(B), 5) * pow((L - L0), 5) * (4 * pow(tan(B), 3) * (tan(B) * tan(B) + 1) - 36 * tan(B) * (tan(B) * tan(B) + 1) + (14 * f * sin(2 * B) * (f - 2)) / pow((f - 1), 2) + (116 * f * sin(B) * (f - 2)) / (cos(B) * pow((f - 1), 2)) - (116 * f * pow(sin(B), 3) * (f - 2)) / (cos(B) * pow((f - 1), 2)))) / (120 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) - (a * sin(B) * (L - L0)) / sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1)) - (a * sin(B) * pow((L - L0), 7) * (840 * pow(sin(B), 4) - 1316 * sin(B) * sin(B) + 479)) / (2520 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) + (a * sin(B) * pow((L - L0), 7) * (662 * sin(B) * sin(B) - 1320 * pow(sin(B), 4) + 720 * pow(sin(B), 6) - 61)) / (720 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) - (a * pow(cos(B), 3) * pow((L - L0), 3) * (2 * tan(B) * (tan(B) * tan(B) + 1) - (f * sin(2 * B) * (f - 2)) / pow((f - 1), 2))) / (6 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) + (a * cos(B) * cos(B) * sin(B) * pow((L - L0), 3) * (tan(B) * tan(B) + (f * cos(B) * cos(B) * (f - 2)) / pow((f - 1), 2) - 1)) / (2 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) - (a * pow(cos(B), 4) * sin(B) * pow((L - L0), 5) * (pow(tan(B), 4) - 18 * pow(tan(B), 2) - (14 * f * cos(B) * cos(B) * (f - 2)) / pow((f - 1), 2) + (58 * f * sin(B) * sin(B) * (f - 2)) / pow((f - 1), 2) + 5)) / (24 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) - (a * f * (L - L0) * (sin(B) - pow(sin(B), 3)) * (f - 2)) / pow((f * f * sin(B) * sin(B) - 2 * f * pow(sin(B), 2) + 1), (3 / 2)) - (a * f * sin(B) * pow((L - L0), 7) * (f - 2) * (1982 * pow(sin(B), 4) - 723 * pow(sin(B), 2) - 2040 * pow(sin(B), 6) + 720 * pow(sin(B), 8) + 61)) / (5040 * pow((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1), (3 / 2))) + (a * f * pow(cos(B), 4) * sin(B) * pow((L - L0), 3) * (f - 2) * (pow(tan(B), 2) + (f * pow(cos(B), 2) * (f - 2)) / pow((f - 1), 2) - 1)) / (6 * pow((f * f * pow(sin(B), 2) - 2 * f * pow(sin(B), 2) + 1), (3 / 2))) - (a * f * pow(cos(B), 6) * sin(B) * pow((L - L0), 5) * (f - 2) * (pow(tan(B), 4) - 18 * pow(tan(B), 2) - (14 * f * pow(cos(B), 2) * (f - 2)) / pow((f - 1), 2) + (58 * f * pow(sin(B), 2) * (f - 2)) / pow((f - 1), 2) + 5)) / (120 * pow((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1), (3 / 2)));
	double dy_dL = (a * cos(B)) / sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1)) - (a * pow(cos(B), 3) * pow((L - L0), 2) * (pow(tan(B), 2) + (f * cos(B) * cos(B) * (f - 2)) / pow((f - 1), 2) - 1)) / (2 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) - (a * pow(cos(B), 7) * pow((L - L0), 6) * (479 * tan(B) * tan(B) - 179 * pow(tan(B), 4) + pow(tan(B), 6) - 61)) / (720 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1))) + (a * pow(cos(B), 5) * pow((L - L0), 4) * (pow(tan(B), 4) - 18 * tan(B) * tan(B) - (14 * f * cos(B) * cos(B) * (f - 2)) / pow((f - 1), 2) + (58 * f * sin(B) * sin(B) * (f - 2)) / pow((f - 1), 2) + 5)) / (24 * sqrt((f * f * sin(B) * sin(B) - 2 * f * sin(B) * sin(B) + 1)));

	double sigmaB = gauss.sigmaB * pi / 180 / 3600;
	double sigmaL = gauss.sigmaL * pi / 180 / 3600;

	gauss.sigmaX = sqrt(dx_dB * dx_dB * sigmaB * sigmaB + dx_dL * dx_dL * sigmaL * sigmaL);
	gauss.sigmaY = sqrt(dy_dB * dy_dB * sigmaB * sigmaB + dy_dL * dy_dL * sigmaL * sigmaL);

}

void GeodeticCalculator::caclGaussCoordBack(PointGauss& gauss, Ellipsoid::Ellipsoid_para ell)
{
	double x = gauss.x, y = gauss.y-gauss.NF, L0 = gauss.L0;
	double a = ell.a, f = ell.f;
	L0 = L0 * pi / 180;
	double p = 3600 * 180 / pi;
	double b = a * (1 - f);
	double e1 = sqrt(a * a - b * b) / a;
	double e2 = sqrt(a * a - b * b) / b;

	double Bf = caclBf(x, a, f);
	double Vf2 = 1.0 + pow(cos(Bf), 2) * e1 * e1 / (1 - e1 * e1);
	double W = sqrt(1 - e1 * e1 * sin(Bf) * sin(Bf));
	double Nf = a / W;
	double tf = tan(Bf);
	double n2 = e2 * e2 * cos(Bf) * cos(Bf);
	double Mf = a * (1 - e1 * e1) / pow(W, 3);

	double m0 = a * (1 - e1 * e1);
	double m2 = 3 * e1 * e1 * m0 / 2;
	double m4 = 5 * e1 * e1 * m2 / 4;
	double m6 = 7 * e1 * e1 * m4 / 6;
	double m8 = 9 * e1 * e1 * m6 / 8;
	double a0 = m0 + m2 / 2 + 3 * m4 / 8 + 5 * m6 / 16 + 35 * m8 / 128;
	double a2 = m2 / 2 + m4 / 2 + 15 * m6 / 32 + 7 * m8 / 16;
	double a4 = m4 / 8 + 3 * m6 / 16 + 7 * m8 / 32;
	double a6 = m6 / 32 + m8 / 16;
	double a8 = m8 / 128;

	//double B = Bf - tf * y * y / (2 * Mf * Nf)
	//	+ tf * pow(y, 4) / (24.0 * Mf * pow(Nf, 3)) * (5 + 3 * pow(tf, 2) + n2 - 9 * n2 * pow(tf, 2))
	//	- tf * pow(y, 6) / (720 * Mf * pow(Nf, 5)) * (61 + 90 * pow(tf, 2) + 45 * pow(tf, 4));
	double B = Bf + tf * (-1 - n2) * y * y / 2 / Nf / Nf
		+ tf * pow(y, 4) * (5 + 3 * tf * tf + 6 * n2 - 6 * n2 * pow(tf, 2) - 3 * n2 * n2 - 9 * n2 * n2 * tf * tf) / 24 / pow(Nf, 4)
		+ tf * pow(y, 6) * (-61 - 90 * tf * tf - 45 * pow(tf, 4) - 107 * n2 + 162 * n2 * tf * tf + 45 * n2 * pow(tf, 4)) / 720 / pow(Nf, 6)
		+ tf * pow(y, 8) * (1385 + 3633 * tf * tf + 4095 * pow(tf, 4) + 1575 * pow(tf, 6)) / 40320 / pow(Nf, 8);
	//B = Bf - tf * y * y / (2 * Mf * Nf) + pow(tf, 3) * pow(y, 4) * (5 + 3 * tf * tf + n2 - 9 * n2 * tf * tf) / (24 * Mf * pow(Nf, 3));
	//double B = Bf - 0.5 * Vf2 * tf * (pow(y / Nf, 2) - (5.0 + 3.0 * n2 + 3.0 * tf * tf - 9.0 * n2 * tf * tf) * pow(y / Nf, 4) / 12.0 + (61.0 + 90.0 * tf * tf + 45 * pow(tf, 4)) * pow(y / Nf, 6) / 360.0);
	double dl = (y / Nf - (1 + 2 * tf * tf + n2) * pow(y / Nf, 3) / (6.0) + (5.0 + 28.0 * tf * tf + 24.0 * pow(tf, 4) + 6.0 * n2 + 8.0 * n2 * tf * tf) * pow(y / Nf, 5) / 120.0) / cos(Bf)
		- (61 + 622 * tf * tf + 1320 * pow(tf, 4) + 720 * pow(tf, 6)) * pow(y / Nf, 7) / 5040 / cos(Bf);
	gauss.B = B * 180 / pi;
	gauss.L = (L0 + dl) * 180 / pi;

	double dB_dx = ((pow(y , 8) * (tan(Bf) * tan(Bf) + 1) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , 4) * (3633 * tan(Bf) * tan(Bf) + 4095 * pow(tan(Bf) , 4) + 1575 * pow(tan(Bf) , 6) + 1385)) / (40320 * pow(a , 8)) + (pow(y , 4) * (tan(Bf) * tan(Bf) + 1) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , 2) * (3 * tan(Bf) * tan(Bf) + (9 * f * f * (cos(4 * Bf) / 8 - 1 / 8) * (f - 2) * (f - 2)) / pow((f - 1) , 4) - (3 * f * f * pow(cos(Bf) , 4) * pow((f - 2) , 2)) / pow((f - 1) , 4) - (6 * f * pow(cos(Bf) , 2) * (f - 2)) / pow((f - 1) , 2) + (6 * f * sin(Bf) * sin(Bf) * (f - 2)) / pow((f - 1) , 2) + 5)) / (24 * pow(a , 4)) + (pow(y , 4) * tan(Bf) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , 2) * (6 * tan(Bf) * (tan(Bf) * tan(Bf) + 1) - (9 * f * f * sin(2 * Bf) * (f - 2) * (f - 2)) / pow((f - 1) , 4) + (6 * f * sin(2 * Bf) * (f - 2)) / pow((f - 1) , 2) + (36 * f * f * cos(Bf) * pow(sin(Bf) , 3) * pow((f - 2) , 2)) / pow((f - 1) , 4) + (12 * f * f * pow(cos(Bf) , 3) * sin(Bf) * pow((f - 2) , 2)) / pow((f - 1) , 4) + (12 * f * sin(Bf) * (f - 2)) / (cos(Bf) * pow((f - 1) , 2)) - (12 * f * pow(sin(Bf) , 3) * (f - 2)) / (cos(Bf) * pow((f - 1) , 2)))) / (24 * pow(a , 4)) - (pow(y , 6) * (tan(Bf) * tan(Bf) + 1) * pow((f *f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , 3) * (90 * tan(Bf) * tan(Bf) + 45 * pow(tan(Bf) , 4) - (107 * f * cos(Bf) *cos(Bf) * (f - 2)) / pow((f - 1) , 2) + (162 * f * sin(Bf) * sin(Bf) * (f - 2)) / pow((f - 1) , 2) - (45 * f * pow(sin(Bf) , 4) * (f - 2)) / ((sin(Bf) * sin(Bf) - 1) * (f - 1) * (f - 1)) + 61)) / (720 * pow(a , 6)) + (pow(y , 8) * tan(Bf) * tan(Bf) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , 4) * (563 * tan(Bf) * tan(Bf) + 615 * pow(tan(Bf) , 4) + 225 * pow(tan(Bf) , 6) + 173)) / (960 * pow(a , 8)) + (y * y * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , 2)) / (2 * a * a * (sin(Bf) * sin(Bf) - 1) * (f - 1) * (f - 1)) - (2 * f * y * y * (f - 2) * (sin(Bf) * sin(Bf) + f * f * pow(sin(Bf) , 4) - 2 * f * pow(sin(Bf) , 4))) / (a * a * (f - 1) * (f - 1)) - (pow(y , 6) * sin(Bf) * tan(Bf) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , 3) * (672 * f * f * pow(sin(Bf) , 4) - 717 * f * f * sin(Bf) * sin(Bf) - 718 * f - 224 * f * f * pow(sin(Bf) , 6) + 359 * f * f + 1434 * f * sin(Bf) * sin(Bf) - 1344 * f * pow(sin(Bf) , 4) + 448 * f * pow(sin(Bf) , 6) + 90)) / (360 * pow(a , 6) * pow(cos(Bf) , 5) * pow((f - 1) , 2)) - (f * pow(y , 6) * cos(Bf) * sin(Bf) * tan(Bf) * (f - 2) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , 2) * (90 * tan(Bf) * tan(Bf) + 45 * pow(tan(Bf) , 4) - (107 * f * cos(Bf) * cos(Bf) * (f - 2)) / pow((f - 1) , 2) + (162 * f * pow(sin(Bf) , 2) * (f - 2)) / pow((f - 1) , 2) - (45 * f * pow(sin(Bf) , 4) * (f - 2)) / ((sin(Bf) * sin(Bf) - 1) * (f - 1) * (f - 1)) + 61)) / (120 * pow(a , 6)) + (f * pow(y , 4) * cos(Bf) * sin(Bf) * tan(Bf) * (f - 2) * (f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) * (3 * tan(Bf) * tan(Bf) + (9 * f * f * (cos(4 * Bf) / 8 - 1 / 8) * (f - 2) * (f - 2)) / pow((f - 1) , 4) - (3 * f * f * pow(cos(Bf) , 4) * (f - 2) * (f - 2)) / pow((f - 1) , 4) - (6 * f * cos(Bf) * cos(Bf) * (f - 2)) / pow((f - 1) , 2) + (6 * f * sin(Bf) * sin(Bf) * (f - 2)) / pow((f - 1) , 2) + 5)) / (6 * pow(a , 4)) + (f * pow(y , 8) * cos(Bf) * sin(Bf) * tan(Bf) * (f - 2) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , 3) * (3633 * tan(Bf) * tan(Bf) + 4095 * pow(tan(Bf) , 4) + 1575 * pow(tan(Bf) , 6) + 1385)) / (5040 * pow(a , 8)) + 1) / (a0 - a2 * cos(2 * Bf) + a4 * cos(4 * Bf) - a6 * cos(6 * Bf) + a8 * cos(8 * Bf));
	double dB_dy = (pow(y , 7) * tan(Bf) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , 4) * (3633 * tan(Bf) * tan(Bf) + 4095 * pow(tan(Bf) , 4) + 1575 * pow(tan(Bf) , 6) + 1385)) / (5040 * pow(a , 8)) + (pow(y , 3) * tan(Bf) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , 2) * (3 * tan(Bf) * tan(Bf) + (9 * f * f * (cos(4 * Bf) / 8 - 1 / 8) * (f - 2) * (f - 2)) / pow((f - 1) , 4) - (3 * f * f * pow(cos(Bf) , 4) * (f - 2) * (f - 2)) / pow((f - 1) , 4) - (6 * f * cos(Bf) * cos(Bf) * (f - 2)) / pow((f - 1) , 2) + (6 * f * sin(Bf) * sin(Bf) * (f - 2)) / pow((f - 1) , 2) + 5)) / (6 * pow(a , 4)) - (pow(y , 5) * tan(Bf) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , 3) * (90 * tan(Bf) * tan(Bf) + 45 * pow(tan(Bf) , 4) - (107 * f * cos(Bf) * cos(Bf) * (f - 2)) / pow((f - 1) , 2) + (162 * f * sin(Bf) * sin(Bf) * (f - 2)) / pow((f - 1) , 2) - (45 * f * pow(sin(Bf) , 4) * (f - 2)) / ((sin(Bf) * sin(Bf) - 1) * (f - 1) * (f - 1)) + 61)) / (120 * pow(a , 6)) - (y * tan(Bf) * (f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) * (f * f - f * f * cos(Bf) * cos(Bf) - 2 * f + 2 * f * cos(Bf) * cos(Bf) + 1)) / (a * a * (f - 1) * (f - 1));
	double dL_dx = -((pow(y , 7) * sin(Bf) * pow((pow(f , 2) * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , (7 / 2)) * (622 * tan(Bf) * tan(Bf) + 1320 * pow(tan(Bf) , 4) + 720 * pow(tan(Bf) , 6) + 61)) / (5040 * pow(a , 7) * pow(cos(Bf) , 2)) - (sin(Bf) * ((y * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , (1 / 2))) / a - (pow(y , 3) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , (3 / 2)) * (2 * tan(Bf) * tan(Bf) - (f * cos(Bf) * cos(Bf) * (f - 2)) / pow((f - 1) , 2) + 1)) / (6 * pow(a , 3)) + (pow(y , 5) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) *sin(Bf) + 1) , (5 / 2)) * (28 * tan(Bf) * tan(Bf) + 24 * pow(tan(Bf) , 4) - (6 * f * cos(Bf) * cos(Bf) * (f - 2)) / pow((f - 1) , 2) - (8 * f * sin(Bf) * sin(Bf) * (f - 2)) / pow((f - 1) , 2) + 5)) / (120 * pow(a , 5)))) / pow(cos(Bf) , 2) - ((pow(y , 5) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , (5 / 2)) * (56 * tan(Bf) * (tan(Bf) * tan(Bf) + 1) + 96 * pow(tan(Bf) , 3) * (tan(Bf) * tan(Bf) + 1) + (6 * f * sin(2 * Bf) * (f - 2)) / pow((f - 1) , 2) - (16 * f * sin(Bf) * (f - 2)) / (cos(Bf) * (f - 1) * (f - 1)) + (16 * f * pow(sin(Bf) , 3) * (f - 2)) / (cos(Bf) * (f - 1) * (f - 1)))) / (120 * pow(a , 5)) - (pow(y , 3) * (4 * tan(Bf) * (tan(Bf) * tan(Bf) + 1) + (f * sin(2 * Bf) * (f - 2)) / pow((f - 1) , 2)) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , (3 / 2))) / (6 * pow(a , 3)) + (f * y * sin(2 * Bf) * (f - 2)) / (2 * a * sqrt((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1))) + (f * pow(y , 5) * cos(Bf) * sin(Bf) * (f - 2) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , (3 / 2)) * (28 * tan(Bf) * tan(Bf) + 24 * pow(tan(Bf) , 4) - (6 * f * cos(Bf) * cos(Bf) * (f - 2)) / pow((f - 1) , 2) - (8 * f * sin(Bf) * sin(Bf) * (f - 2)) / pow((f - 1) , 2) + 5)) / (24 * pow(a , 5)) - (f * pow(y , 3) * cos(Bf) * sin(Bf) * (f - 2) * sqrt((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1)) * (2 * tan(Bf) * tan(Bf) - (f * cos(Bf) * cos(Bf) * (f - 2)) / pow((f - 1) , 2) + 1)) / (2 * pow(a , 3))) / cos(Bf) + (pow(y , 7) * tan(Bf) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , (7 / 2)) * (1631 * tan(Bf) * tan(Bf) + 2400 * pow(tan(Bf) , 4) + 1080 * pow(tan(Bf) , 6) + 311)) / (1260 * pow(a , 7) * cos(Bf)) + (f * pow(y , 7) * sin(Bf) * (f - 2) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , (5 / 2)) * (622 * tan(Bf) * tan(Bf) + 1320 * pow(tan(Bf) , 4) + 720 * pow(tan(Bf) , 6) + 61)) / (720 * pow(a , 7))) / (a0 - a2 * cos(2 * Bf) + a4 * cos(4 * Bf) - a6 * cos(6 * Bf) + a8 * cos(8 * Bf));
	double dL_dy = (sqrt((f *f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1)) / a - (y * y * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , (3 / 2)) * (2 * tan(Bf) * tan(Bf) - (f * cos(Bf) * cos(Bf) * (f - 2)) / pow((f - 1) , 2) + 1)) / (2 * pow(a , 3)) + (pow(y , 4) * pow((f * f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , (5 / 2)) * (28 * tan(Bf) * tan(Bf) + 24 * pow(tan(Bf) , 4) - (6 * f * cos(Bf) * cos(Bf) * (f - 2)) / pow((f - 1) , 2) - (8 * f * sin(Bf) * sin(Bf) * (f - 2)) / pow((f - 1) , 2) + 5)) / (24 * pow(a , 5))) / cos(Bf) - (pow(y , 6) * pow((f *f * sin(Bf) * sin(Bf) - 2 * f * sin(Bf) * sin(Bf) + 1) , (7 / 2)) * (622 * tan(Bf) * tan(Bf) + 1320 * pow(tan(Bf) , 4) + 720 * pow(tan(Bf) , 6) + 61)) / (720 * pow(a , 7) * cos(Bf));

	double sigmaX = gauss.sigmaX;
	double sigmaY = gauss.sigmaY;
	gauss.sigmaB = sqrt(dB_dx * dB_dx * sigmaX * sigmaX + dB_dy * dB_dy * sigmaY * sigmaY) * 180 * 3600 / pi;
	gauss.sigmaL = sqrt(dL_dx * dL_dx * sigmaX * sigmaX + dL_dy * dL_dy * sigmaY * sigmaY) * 180 * 3600 / pi;
}
//迭代求Bf
double GeodeticCalculator::caclBf(double X, double a, double f)
{
	double b = a * (1 - f);
	double e1 = sqrt(a * a - b * b) / a;
	double m0 = a * (1 - e1 * e1);
	double m2 = 3 * e1 * e1 * m0 / 2;
	double m4 = 5 * e1 * e1 * m2 / 4;
	double m6 = 7 * e1 * e1 * m4 / 6;
	double m8 = 9 * e1 * e1 * m6 / 8;
	double a0 = m0 + m2 / 2 + 3 * m4 / 8 + 5 * m6 / 16 + 35 * m8 / 128;
	double a2 = m2 / 2 + m4 / 2 + 15 * m6 / 32 + 7 * m8 / 16;
	double a4 = m4 / 8 + 3 * m6 / 16 + 7 * m8 / 32;
	double a6 = m6 / 32 + m8 / 16;
	double a8 = m8 / 128;

	double Bf0, Bf1;
	double FBf;
	Bf0 = X / a0;
	while (1)
	{
		FBf = -a2 * sin(2 * Bf0) / 2.0 + a4 * sin(4 * Bf0) / 4.0 - a6 / 6.0 * sin(6 * Bf0)+ a8 / 8.0 * sin(Bf0 * 8);
		Bf1 = (X - FBf) / a0;
		if (abs(Bf1 - Bf0) < 1e-25)
		{
			return Bf1;
		}
		else
		{
			Bf0 = Bf1;
		}
	}
}
double GeodeticCalculator::caclBf2(double X, double a, double f)
{
	double b = a * (1 - f);
	double e1 = sqrt(a * a - b * b) / a;
	double a0=1+3.0/4.0*e1*e1+45.0/64.0*e1*e1*e1*e1+350.0/512.0*e1*e1*e1*e1*e1*e1+11025.0/16384.0*e1*e1*e1*e1*e1*e1*e1*e1;
	double B0 = X / (a0 * a * (1 - e1 * e1));
	double K0=0.5*(3.0/4.0*e1*e1+45.0/64.0*e1*e1*e1*e1+350.0/512.0*e1*e1*e1*e1*e1*e1 + 11025.0 / 16384.0 * e1 * e1 * e1 * e1 * e1 * e1 * e1 * e1);
	double K2 = -1 / 3.0 * (63.0 / 64.0 * pow(e1, 4) + 1108.0 / 512.0 * pow(e1, 6) + 58239.0 / 16384.0 * pow(e1, 8));
	double K4 = 1 / 3.0 * (604 / 672 * pow(e1, 6) + 68484.0 / 16384.0 * pow(e1, 8));
    double K6 = -1 / 3.0 * (26328.0 / 16384.0 * pow(e1, 8));
	double Bf=B0+sin(2*B0)*(K0 +sin(B0)* sin(B0)*(K2+sin(B0)* sin(B0)*(K4+sin(B0)* sin(B0)*K6)));
	return Bf;
}

void GeodeticCalculator::errorPropagation(MatrixXd J, MatrixXd Dcoord, MatrixXd B, MatrixXd Dxx, MatrixXd& rst)
{
	MatrixXd rst1 = J * Dcoord * J.transpose();
	MatrixXd rst2 = B * Dxx * B.transpose();
	rst = rst1 + rst2;
}


int GeodeticCalculator::caclSignede(double L, int width)
{
	int signede = 0;
	if (width == 6)
	{
		signede = (L / 6) + 1;
	}
	else if(width==3)
	{
		signede = (L + 1.5) / 3;
	}
	
	return signede;
}
int GeodeticCalculator::caclCentralMeridian(double signede, int width)
{
	int centralMeridian = 0;
	if (width == 6)
	{
		centralMeridian = signede * 6 -3;
	}
	else if (width == 3)
	{
		centralMeridian = signede * 3;
	}
	return centralMeridian;
}
/**
 * @brief 计算子午线收敛角
 *
 * @param mca 输入输出参数，包含转换后的子午线收敛角信息，包括：
 *            - B: 纬度（单位：度）
 *            - L: 经度（单位：度）
 *            - L0: 中央子午线经度（单位：度）
 *            - sigmaB: 纬度不确定度（单位：秒）
 *            - sigmaL: 经度不确定度（单位：秒）
 * 		  - mca: 子午线收敛角（单位：度）
 * 		  - sigmaMCA: 子午线收敛角不确定度（单位：秒）
 * @param ell 椭球参数结构体，包含：
 *            - a: 椭球的长半轴（单位：米）
 *            - f: 椭球的扁率（无单位）
 */
void GeodeticCalculator::caclMeridianConvergenceAngle(PointMCA& mca, Ellipsoid::Ellipsoid_para ell)
{
	double L0 = mca.L0;
	double L = mca.L;
	double B = mca.B;
	double a = ell.a, f = ell.f;
	double b = a * (1 - f);
	double e2 = sqrt(a * a - b * b) / b;
	double l = L-L0;
	l = deg2rad(l);
	B = deg2rad(B);
	double t = tan(B);
	double n = e2 * cos(B); 
	double eta_2 = n * n;
	double gamma = sin(B) * l * ((1 + (1 / 3.0) * (1 + 3 * eta_2+2 * eta_2 * eta_2) * (cos(B) * cos(B)) * (l * l)) +(1 / 15.0) * (2 - (t * t)) * pow(cos(B), 4) * pow(l, 4));
	mca.mca = rad2deg(gamma);
	double sigmaB = mca.sigmaB * pi / 180 / 3600;
	double sigmaL = mca.sigmaL * pi / 180 / 3600;
	double dmca_dB = cos(B) * (l) * (pow((l) , 4) * (pow(sin(B) , 4 )/ 5 - pow(sin(B) , 2) / 3 + 2 / 15) + pow(cos(B) , 2) * pow((l) , 2) * (e2 *e2 * cos(B) * cos(B) + (2 * pow(e2 , 4) * pow(cos(B) , 4)) / 3 + 1 / 3) + 1) - sin(B) * (l) * ((sin(2 * B) * pow((l) , 4)) / 15 + 4 * (sin(2 * B) / 60 + sin(4 * B) / 40) * pow((l) , 4) + 2 * cos(B) * sin(B) * pow((l) , 2 )* (e2 * e2 * cos(B) * cos(B) + (2 * pow(e2 , 4) * pow(cos(B) , 4)) / 3 + 1 / 3) + (2 * e2 * e2 * pow(cos(B) , 3) * sin(B) * pow((l) , 2) * (4 * e2 * e2 * cos(B) * cos(B) + 3)) / 3);
	double dmca_dL = sin(B) * (pow((l) , 4) * (pow(sin(B) , 4) / 5 - sin(B) * sin(B) / 3 + 2 / 15) + cos(B) *cos(B) * pow((l) , 2) * (e2 * e2 * cos(B) * cos(B) + (2 * pow(e2 , 4) * pow(cos(B) , 4)) / 3 + 1 / 3) + 1) + sin(B) * (l) * (4 * pow((l) , 3) * (pow(sin(B) , 4) / 5 - sin(B) * sin(B) / 3 + 2 / 15) + cos(B) * cos(B) * (2 * l) * (e2 *e2 * cos(B) * cos(B) + (2 * pow(e2 , 4) * pow(cos(B) , 4)) / 3 + 1 / 3));
	mca.sigmaMCA = rad2deg(sqrt(dmca_dB * dmca_dB * sigmaB * sigmaB + dmca_dL * dmca_dL * sigmaL * sigmaL))/3600;
}

void GeodeticCalculator::caclAstroAngleByGeo(PointAz& pointAz)
{
	//地面
	double alpha = deg2rad(pointAz.alpha), xi = deg2rad(pointAz.xi / 3600), eta = deg2rad(pointAz.eta / 3600), Z = deg2rad(pointAz.Z);
	double L = deg2rad(pointAz.L), lamda = deg2rad(pointAz.lamda), fai = deg2rad(pointAz.fai);
	double A0 = alpha, A = alpha;
	int n = 0;
	if (eta == 0 && xi == 0)
	{
		A = alpha - (lamda - L) * sin(fai);
	}
	else
	{
		do
		{
			A0 = A;
			A = alpha - (lamda - L) * sin(fai) - (xi * sin(A0) - eta * cos(A0)) / tan(Z);
			n++;
		} while (abs(A - A0) > 1e-10 && n < 100);
	}
	pointAz.A = rad2deg(A);
	//精度
	double sigmaL = deg2rad(pointAz.sigmaL / 3600), sigmaLamda = deg2rad(pointAz.sigmaLamda / 3600), sigmaAlpha = deg2rad(pointAz.sigmaAlpha / 3600), sigmaFai = deg2rad(pointAz.sigmaFai / 3600);
	double sigmaEta = deg2rad(pointAz.sigmaEta / 3600), sigmaXi = deg2rad(pointAz.sigmaXi / 3600), sigmaZ = deg2rad(pointAz.sigmaZ / 3600);
	double DA = (sigmaL * sigmaL) * pow(sin(fai), 2.0) + (sigmaLamda * sigmaLamda) * pow(sin(fai), 2.0) + sigmaAlpha * sigmaAlpha + (sigmaFai * sigmaFai) * pow(cos(fai), 2.0) * pow(L - lamda, 2.0) + (sigmaEta * sigmaEta) * pow(cos(A), 2.0) * 1.0 / pow(tan(Z), 2.0) + (sigmaXi * sigmaXi) * pow(sin(A), 2.0) * 1.0 / pow(tan(Z), 2.0) + (sigmaZ * sigmaZ) * 1.0 / pow(tan(Z), 4.0) * pow(eta * cos(A) - xi * sin(A), 2.0) * pow(pow(tan(Z), 2.0) + 1.0, 2.0);
	pointAz.sigmaA = rad2deg(sqrt(DA))*3600;

	//海
	double dHfai = deg2rad(-0.000171 * pointAz.H1 * sin(2 * fai) / 3600);
	double fai_sea = fai + dHfai;
	pointAz.fai_sea = rad2deg(fai_sea);
	double h2 = pointAz.H2 + pointAz.Zeta;
	double dHalpha = deg2rad(- 0.000108 * h2 * cos(fai) * cos(fai) * sin(2 * alpha) / 3600);
	double alpha_sea = alpha + dHalpha;
	pointAz.alpha_sea = rad2deg(alpha_sea);

	double A_sea = alpha_sea, A0_sea = alpha_sea;
	if (eta == 0 && xi == 0)
	{
		A_sea = alpha_sea - (lamda - L) * sin(fai_sea);
	}
	else
	{
		do
		{
			A0_sea = A_sea;
			A_sea = alpha_sea - (lamda - L) * sin(fai_sea) - (xi * sin(A0_sea) - eta * cos(A0_sea)) / tan(Z);
			n++;
		} while (abs(A_sea - A0_sea) > 1e-10 && n < 100);
	}
	pointAz.A_sea = rad2deg(A_sea);

	double H1 = pointAz.H1, H2 = pointAz.H2, Zeta = pointAz.Zeta;
	double sigmaH1 = (pointAz.sigmaH1), sigmaH2 = (pointAz.sigmaH2), sigmaZeta = (pointAz.sigmaZeta);
	double DA_SEA=(sigmaAlpha * sigmaAlpha)* pow(cos(alpha * 2.0) * pow(cos(fai), 2.0) * (H2 * 1.08E-4 + Zeta * 1.08E-4) * 2.0 - 1.0, 2.0) + (sigmaFai * sigmaFai) * pow(cos(fai - H1 * sin(fai * 2.0) * 1.71E-4) * (H1 * cos(fai * 2.0) * 3.42E-4 - 1.0) * (L - lamda) - sin(alpha * 2.0) * cos(fai) * sin(fai) * (H2 * 1.08E-4 + Zeta * 1.08E-4) * 2.0, 2.0) + (sigmaL * sigmaL) * pow(sin(fai - H1 * sin(fai * 2.0) * 1.71E-4), 2.0) + (sigmaLamda * sigmaLamda) * pow(sin(fai - H1 * sin(fai * 2.0) * 1.71E-4), 2.0) + (sigmaEta * sigmaEta) * pow(cos(A_sea), 2.0) * 1.0 / pow(tan(Z), 2.0) + (sigmaXi * sigmaXi) * pow(sin(A_sea), 2.0) * 1.0 / pow(tan(Z), 2.0) + (sigmaH2 * sigmaH2) * pow(sin(alpha * 2.0), 2.0) * pow(cos(fai), 4.0) * 1.1664E-8 + (sigmaZeta * sigmaZeta) * pow(sin(alpha * 2.0), 2.0) * pow(cos(fai), 4.0) * 1.1664E-8 + (sigmaZ * sigmaZ) * 1.0 / pow(tan(Z), 4.0) * pow(eta * cos(A_sea) - xi * sin(A_sea), 2.0) * pow(pow(tan(Z), 2.0) + 1.0, 2.0) + (sigmaH1 * sigmaH1) * pow(sin(fai * 2.0), 2.0) * pow(cos(fai - H1 * sin(fai * 2.0) * 1.71E-4), 2.0) * pow(L - lamda, 2.0) * 2.9241E-8;
	pointAz.sigmaA_sea = rad2deg(sqrt(DA_SEA))*3600;
}

void GeodeticCalculator::caclGeoAngleByAstro(PointAz& pointAz)
{
	double A = deg2rad(pointAz.A), xi = deg2rad(pointAz.xi / 3600), eta = deg2rad(pointAz.eta / 3600), Z = deg2rad(pointAz.Z);
	double L = deg2rad(pointAz.L), lamda = deg2rad(pointAz.lamda), fai = deg2rad(pointAz.fai);
	double alpha = A + (lamda - L) * sin(fai) + (xi * sin(A) - eta * cos(A)) / tan(Z);
	pointAz.alpha = rad2deg(alpha);

	double sigmaL = deg2rad(pointAz.sigmaL / 3600), sigmaLamda = deg2rad(pointAz.sigmaLamda / 3600), sigmaA = deg2rad(pointAz.sigmaA / 3600), sigmaFai = deg2rad(pointAz.sigmaFai / 3600);
	double sigmaEta = deg2rad(pointAz.sigmaEta / 3600), sigmaXi = deg2rad(pointAz.sigmaXi / 3600), sigmaZ = deg2rad(pointAz.sigmaZ / 3600);
	double Dalpha = (sigmaL * sigmaL) * pow(sin(fai), 2.0) + (sigmaLamda * sigmaLamda) * pow(sin(fai), 2.0) + (sigmaA * sigmaA) * pow((eta * sin(A) + xi * cos(A)) / tan(Z) + 1.0, 2.0) + (sigmaFai * sigmaFai) * pow(cos(fai), 2.0) * pow(L - lamda, 2.0) + (sigmaEta * sigmaEta) * pow(cos(A), 2.0) * 1.0 / pow(tan(Z), 2.0) + (sigmaXi * sigmaXi) * pow(sin(A), 2.0) * 1.0 / pow(tan(Z), 2.0) + (sigmaZ * sigmaZ) * 1.0 / pow(tan(Z), 4.0) * pow(eta * cos(A) - xi * sin(A), 2.0) * pow(pow(tan(Z), 2.0) + 1.0, 2.0);
	pointAz.sigmaAlpha = rad2deg(sqrt(Dalpha)) * 3600;
}


void GeodeticCalculator::caclCoordAngleByGeo(PointCoordAz& pointCoordAz)
{
	double a = pointCoordAz.ell.a;
	double f = pointCoordAz.ell.f;
	double b = a * (1 - f);
	double c = a * (a / b);
	double e1 = sqrt((a / b) * (a / b) - 1.0);
	double B1 = deg2rad(pointCoordAz.p1.B);
	double n = e1 * cos(B1);
	double V = sqrt(1 + n * n);
	double M = c / (V * V * V) ;
	double N = c / V;
	double R = sqrt(M * N);
	double t = tan(B1);

	GaussForward(pointCoordAz.p1, pointCoordAz.ell);
	GaussForward(pointCoordAz.p2, pointCoordAz.ell);

	double y1 = pointCoordAz.p1.y - pointCoordAz.p1.NF;
	double y2 = pointCoordAz.p2.y - pointCoordAz.p2.NF;
	double x1 = pointCoordAz.p1.x;
	double x2 = pointCoordAz.p2.x;
	double ym = (y1+y2)/2.0;
	double derta = (x1 - x2) * (2.0 * y1 + y2 - pow(ym, 3) / (R * R)) / 6.0 / (R * R) + n * n * t * (y1 - y2) * ym * ym / pow(R, 3);
	
	double sigmaB1 = pointCoordAz.p1.sigmaB * pi / 180 / 3600;
	double sigmaX1 = pointCoordAz.p1.sigmaX;
	double sigmaY1 = pointCoordAz.p1.sigmaY;
	double sigmaX2 = pointCoordAz.p2.sigmaX;
	double sigmaY2 = pointCoordAz.p2.sigmaY;
	double Dderta = (sigmaB1 * sigmaB1)* pow(1.0 / (R * R * R) * pow(cos(B1), 2.0) * pow(y1 / 2.0 + y2 / 2.0, 2.0) * ((a * a) * 1.0 / (b * b) - 1.0) * (pow(tan(B1), 2.0) + 1.0) * (y1 - y2) - 1.0 / (R * R * R) * cos(B1) * sin(B1) * tan(B1) * pow(y1 / 2.0 + y2 / 2.0, 2.0) * ((a * a) * 1.0 / (b * b) - 1.0) * (y1 - y2) * 2.0, 2.0) + (sigmaY1 * sigmaY1) * pow(1.0 / (R * R) * (x1 - x2) * (1.0 / (R * R) * pow(y1 / 2.0 + y2 / 2.0, 2.0) * (3.0 / 2.0) - 2.0) * (-1.0 / 6.0) + 1.0 / (R * R * R) * pow(cos(B1), 2.0) * tan(B1) * pow(y1 / 2.0 + y2 / 2.0, 2.0) * ((a * a) * 1.0 / (b * b) - 1.0) + 1.0 / (R * R * R) * pow(cos(B1), 2.0) * tan(B1) * (y1 / 2.0 + y2 / 2.0) * ((a * a) * 1.0 / (b * b) - 1.0) * (y1 - y2), 2.0) + (sigmaY2 * sigmaY2) * pow((1.0 / (R * R) * (x1 - x2) * (1.0 / (R * R) * pow(y1 / 2.0 + y2 / 2.0, 2.0) * (3.0 / 2.0) - 1.0)) / 6.0 + 1.0 / (R * R * R) * pow(cos(B1), 2.0) * tan(B1) * pow(y1 / 2.0 + y2 / 2.0, 2.0) * ((a * a) * 1.0 / (b * b) - 1.0) - 1.0 / (R * R * R) * pow(cos(B1), 2.0) * tan(B1) * (y1 / 2.0 + y2 / 2.0) * ((a * a) * 1.0 / (b * b) - 1.0) * (y1 - y2), 2.0) + (1.0 / (R * R * R * R) * (sigmaX1 * sigmaX1) * pow(y1 * 2.0 + y2 - 1.0 / (R * R) * pow(y1 / 2.0 + y2 / 2.0, 3.0), 2.0)) / 3.6E1 + (1.0 / (R * R * R * R) * (sigmaX2 * sigmaX2) * pow(y1 * 2.0 + y2 - 1.0 / (R * R) * pow(y1 / 2.0 + y2 / 2.0, 3.0), 2.0)) / 3.6E1;
	
	pointCoordAz.alphaCoord = pointCoordAz.A - pointCoordAz.p1.mca + rad2deg(derta);
	double sigmaA = pointCoordAz.sigmaA * pi / 180 / 3600;
	double sigmaMCA = pointCoordAz.p1.sigmaMCA * pi / 180 / 3600;
	pointCoordAz.sigmaAlphaCoord = rad2deg(sqrt(pow(sigmaA, 2.0) + pow(sigmaMCA, 2.0) + Dderta))*3600;
}


void GeodeticCalculator::caclGeoAngleByCoord(PointCoordAz& pointCoordAz)
{
	double a = pointCoordAz.ell.a;
	double f = pointCoordAz.ell.f;
	double b = a * (1 - f);
	double c = a * (a / b);
	double e1 = sqrt((a / b) * (a / b) - 1.0);
	double B1 = deg2rad(pointCoordAz.p1.B);
	double n = e1 * cos(B1);
	double V = sqrt(1 + n * n);
	double M = c / (V * V * V);
	double N = c / V;
	double R = sqrt(M * N);
	double t = tan(B1);

	GaussForward(pointCoordAz.p1, pointCoordAz.ell);
	GaussForward(pointCoordAz.p2, pointCoordAz.ell);

	double y1 = pointCoordAz.p1.y - pointCoordAz.p1.NF;
	double y2 = pointCoordAz.p2.y - pointCoordAz.p2.NF;
	double x1 = pointCoordAz.p1.x;
	double x2 = pointCoordAz.p2.x;
	double ym = (y1 + y2) / 2.0;
	double derta = (x1 - x2) * (2.0 * y1 + y2 - pow(ym, 3) / (R * R)) / 6.0 / (R * R) + n * n * t * (y1 - y2) * ym * ym / pow(R, 3);

	double sigmaB1 = pointCoordAz.p1.sigmaB * pi / 180 / 3600;
	double sigmaX1 = pointCoordAz.p1.sigmaX;
	double sigmaY1 = pointCoordAz.p1.sigmaY;
	double sigmaX2 = pointCoordAz.p2.sigmaX;
	double sigmaY2 = pointCoordAz.p2.sigmaY;
	double Dderta = (sigmaB1 * sigmaB1) * pow(1.0 / (R * R * R) * pow(cos(B1), 2.0) * pow(y1 / 2.0 + y2 / 2.0, 2.0) * ((a * a) * 1.0 / (b * b) - 1.0) * (pow(tan(B1), 2.0) + 1.0) * (y1 - y2) - 1.0 / (R * R * R) * cos(B1) * sin(B1) * tan(B1) * pow(y1 / 2.0 + y2 / 2.0, 2.0) * ((a * a) * 1.0 / (b * b) - 1.0) * (y1 - y2) * 2.0, 2.0) + (sigmaY1 * sigmaY1) * pow(1.0 / (R * R) * (x1 - x2) * (1.0 / (R * R) * pow(y1 / 2.0 + y2 / 2.0, 2.0) * (3.0 / 2.0) - 2.0) * (-1.0 / 6.0) + 1.0 / (R * R * R) * pow(cos(B1), 2.0) * tan(B1) * pow(y1 / 2.0 + y2 / 2.0, 2.0) * ((a * a) * 1.0 / (b * b) - 1.0) + 1.0 / (R * R * R) * pow(cos(B1), 2.0) * tan(B1) * (y1 / 2.0 + y2 / 2.0) * ((a * a) * 1.0 / (b * b) - 1.0) * (y1 - y2), 2.0) + (sigmaY2 * sigmaY2) * pow((1.0 / (R * R) * (x1 - x2) * (1.0 / (R * R) * pow(y1 / 2.0 + y2 / 2.0, 2.0) * (3.0 / 2.0) - 1.0)) / 6.0 + 1.0 / (R * R * R) * pow(cos(B1), 2.0) * tan(B1) * pow(y1 / 2.0 + y2 / 2.0, 2.0) * ((a * a) * 1.0 / (b * b) - 1.0) - 1.0 / (R * R * R) * pow(cos(B1), 2.0) * tan(B1) * (y1 / 2.0 + y2 / 2.0) * ((a * a) * 1.0 / (b * b) - 1.0) * (y1 - y2), 2.0) + (1.0 / (R * R * R * R) * (sigmaX1 * sigmaX1) * pow(y1 * 2.0 + y2 - 1.0 / (R * R) * pow(y1 / 2.0 + y2 / 2.0, 3.0), 2.0)) / 3.6E1 + (1.0 / (R * R * R * R) * (sigmaX2 * sigmaX2) * pow(y1 * 2.0 + y2 - 1.0 / (R * R) * pow(y1 / 2.0 + y2 / 2.0, 3.0), 2.0)) / 3.6E1;

	pointCoordAz.A = pointCoordAz.alphaCoord + pointCoordAz.p1.mca - rad2deg(derta);
	double alphaCoord = pointCoordAz.alphaCoord * pi / 180 / 3600;
	double sigmaMCA = pointCoordAz.p1.sigmaMCA * pi / 180 / 3600;
	pointCoordAz.sigmaAlphaCoord = rad2deg(sqrt(pow(alphaCoord, 2.0) + pow(sigmaMCA, 2.0) + Dderta)) * 3600;

}

