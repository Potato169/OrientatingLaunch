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
	//������������
	static double caclBf(double X, double a, double f);
	static double caclBf2(double X, double a, double f);
	//gauss��������
	static void caclGaussCoordForward(PointGauss& gauss, Ellipsoid::Ellipsoid_para ell);
	static void caclGaussCoordBack(PointGauss& gauss, Ellipsoid::Ellipsoid_para ell);
	static int caclSignede(double L, int width);
	static int caclCentralMeridian(double signede, int width);

	//����
	static void errorPropagation(MatrixXd J, MatrixXd Dcoord, MatrixXd B, MatrixXd Dxx, MatrixXd& rst);

public:
	static void caclMeridianConvergenceAngle(PointMCA& mca, Ellipsoid::Ellipsoid_para ell);
	/*��˹ͶӰ*/
	static void GaussForward(PointGauss& gauss, Ellipsoid::Ellipsoid_para ell);//��˹����
	static void GaussBack(PointGauss& gauss, Ellipsoid::Ellipsoid_para ell);//��˹����
	//���ķ�λ��-��ط�λ��
	static void caclAstroAngleToGeo(PointAz& pointAz);
	static void caclGeoAngleToAstro(PointAz& pointAz);

	//��ط�λ��-���귽λ��
	static void caclCoordAngleByGeo(PointCoordAz& pointCoordAz);
	static void caclGeoAngleByCoord(PointCoordAz& pointCoordAz);
};

#endif