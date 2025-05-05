#ifndef GEODETICCALCULATOR_H
#define GEODETICCALCULATOR_H
#include <vector>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include "AngleConversion.h"
#include "CoordSystem.h"
#include "FileOperate.h"
#define pi acos(-1)
using namespace std;
using namespace Eigen;
using namespace AngleConversion;
using namespace CoordSystem;
using namespace FileOperate;
class GeodeticCalculator
{
private:
	//其它辅助函数
	static double caclBf(double X, double a, double f);
	static double caclBf2(double X, double a, double f);
	//gauss辅助函数
	static void caclGaussCoordForward(PointGauss& gauss, Ellipsoid::Ellipsoid_para ell);
	static void caclGaussCoordBack(PointGauss& gauss, Ellipsoid::Ellipsoid_para ell);
	static int caclSignede(double L, int width);
	static int caclCentralMeridian(double signede, int width);

	//误差传递
	static void errorPropagation(MatrixXd J, MatrixXd Dcoord, MatrixXd B, MatrixXd Dxx, MatrixXd& rst);

public:
	static void caclMeridianConvergenceAngle(PointMCA& mca, Ellipsoid::Ellipsoid_para ell);
	/*高斯投影*/
	static void GaussForward(PointGauss& gauss, Ellipsoid::Ellipsoid_para ell);//高斯正算
	static void GaussBack(PointGauss& gauss, Ellipsoid::Ellipsoid_para ell);//高斯反算
	//天文方位角-大地方位角
	static void caclAstroAngleToGeo(PointAz& pointAz);
	static void caclGeoAngleToAstro(PointAz& pointAz);

	//大地方位角-坐标方位角
	static void caclCoordAngleByGeo(PointCoordAz& pointCoordAz);
	static void caclGeoAngleByCoord(PointCoordAz& pointCoordAz);
};

#endif