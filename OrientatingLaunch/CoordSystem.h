#ifndef COORDSYSTEM_H
#define COORDSYSTEM_H
#include <string>
#include <vector>
using namespace std;
namespace CoordSystem
{
	struct Ellipsoid {
		struct Ellipsoid_para {
			double a = 0;  // 半长轴
			double f = 0;  // 扁率
			string name = "";
			Ellipsoid_para(double a, double f, string name) : a(a), f(f), name(name) {}
		};
		static const Ellipsoid_para bj54;
		static const Ellipsoid_para xa80;
		static const Ellipsoid_para wgs84;
		static const Ellipsoid_para cgcs2000;
		static const Ellipsoid_para grs80;
		static const Ellipsoid_para pz90;
	};
	//sigmaXYZ的单位为米
	struct PointXYZ
	{
		string PointName = "";
		double X = 0;
		double Y = 0;
		double Z = 0;
		double sigmaX = 0;
		double sigmaY = 0;
		double sigmaZ = 0;
		PointXYZ() {}
		PointXYZ(double X, double Y) :X(X), Y(Y) {};
		PointXYZ(string PointName, double X, double Y) :PointName(PointName), X(X), Y(Y) {};
		PointXYZ(double X, double Y, double sigmaX, double sigmaY) :X(X), Y(Y), sigmaX(sigmaX), sigmaY(sigmaY) {};
		PointXYZ(string PointName, double X, double Y, double sigmaX, double sigmaY) :PointName(PointName), X(X), Y(Y), sigmaX(sigmaX), sigmaY(sigmaY) {};
		PointXYZ(double X, double Y, double Z) : X(X), Y(Y), Z(Z) {}
		PointXYZ(string PointName, double X, double Y, double Z) :PointName(PointName), X(X), Y(Y), Z(Z) {}
		PointXYZ(double X, double Y, double Z, double sigmaX, double sigmaY, double sigmaZ) :X(X), Y(Y), Z(Z), sigmaX(sigmaX), sigmaY(sigmaY), sigmaZ(sigmaZ) {}
		PointXYZ(string PointName, double X, double Y, double Z, double sigmaX, double sigmaY, double sigmaZ) :PointName(PointName), X(X), Y(Y), Z(Z), sigmaX(sigmaX), sigmaY(sigmaY), sigmaZ(sigmaZ) {}
	};
	//sigmaBL的单位为秒,sigmaH的单位为米
	struct PointBLH
	{
		string PointName = "";
		double B = 0;
		double L = 0;
		double H = 0;
		double sigmaB = 0;
		double sigmaL = 0;
		double sigmaH = 0;
		PointBLH() {}
		PointBLH(double B, double L, double H = 0) :B(B), L(L), H(H) {}
		PointBLH(string PointName, double B, double L, double H = 0) :PointName(PointName), B(B), L(L), H(H) {}
		PointBLH(double B, double L, double H, double sigmaB, double sigmaL, double sigmaH) :B(B), L(L), H(H), sigmaB(sigmaB), sigmaL(sigmaL), sigmaH(sigmaH) {}
		PointBLH(string PointName, double B, double L, double H, double sigmaB, double sigmaL, double sigmaH) :PointName(PointName),B(B), L(L), H(H), sigmaB(sigmaB), sigmaL(sigmaL), sigmaH(sigmaH) {}
	};
	struct PointUTM
	{
		string PointName = "";
		double x = 0, y = 0;
		double B = 0, L = 0;
		double sigmaX = 0, sigmaY = 0;
		double sigmaB = 0, sigmaL = 0;
		string signeded="";
		PointUTM() {};
		PointUTM(string PointName,double x, double y, double B, double L, string signeded) :PointName(PointName),x(x), y(y), B(B), L(L), signeded(signeded) {};
		PointUTM(double x, double y, double B, double L, string signeded) :x(x), y(y), B(B), L(L), signeded(signeded) {};
	};
	struct PointGauss
	{
		string PointName = "";
		double B = 0, L = 0;
		double x=0, y=0;
		double sigmaX = 0, sigmaY = 0;
		double sigmaB = 0, sigmaL = 0;
		double L0 = 0; // 中央子午线
		int signeded = 0; // 带号
		double mca = 0;// 子午线收敛角
		double sigmaMCA = 0;
		double NF = 0;//北向偏移量
		int width = 0;//投影带宽度
		PointGauss() {};
		PointGauss(double B, double L, double x, double y, double L0, int signeded, double mca) 
			:B(B), L(L), x(x), y(y), L0(L0), signeded(signeded), mca(mca) {};
		PointGauss(string PointName,double B, double L, double x, double y, double L0, int signeded, double mca) 
			:PointName(PointName),B(B), L(L), x(x), y(y), L0(L0), signeded(signeded), mca(mca) {}
		PointGauss(PointUTM utm) 
			:PointName(utm.PointName), B(utm.B), L(utm.L), x(utm.x), y(utm.y), sigmaB(utm.sigmaB), sigmaL(utm.sigmaL), sigmaX(utm.sigmaX), sigmaY(utm.sigmaY) {};//构造函数
	};
	struct PointMCA
	{
		string PointName = "";
		double B = 0, L = 0, L0=0;
		double sigmaB = 0, sigmaL = 0;
		double mca = 0;
		double sigmaMCA = 0;
		PointMCA() {};
		PointMCA(string PointName, double B, double L, double L0) :PointName(PointName), B(B), L(L), L0(L0){};
		PointMCA(PointGauss gauss) 
			:PointName(gauss.PointName), B(gauss.B), L(gauss.L), L0(gauss.L0), sigmaB(gauss.sigmaB), sigmaL(gauss.sigmaL) {};
	};
	
	struct PointAz
	{
		string pointName = "";
		double alpha = 0;
		double xi = 0;
		double eta = 0;
		double Z = 90;
		double L = 0;
		double lamda = 0;
		double fai = 0;
		double A = 0;
		
		double Zeta = 0;
		double H1 = 0;
		double H2 = 0;
		double fai_sea = 0;
		double alpha_sea = 0;
		double A_sea = 0;

		double sigmaZeta = 0;
		double sigmaH1 = 0;
		double sigmaH2 = 0;
		double sigmaA_sea = 0;

		double sigmaAlpha = 0;
		double sigmaFai = 0;
		double sigmaLamda = 0;
		double sigmaL = 0;
		double sigmaZ = 0;
		double sigmaXi = 0;
		double sigmaEta = 0;
		double sigmaA = 0;
		PointAz() {};
		PointAz(string pointName, double alpha, double xi, double eta, double Z, double L, double lamda, double fai)
			: pointName(pointName), alpha(alpha), xi(xi), eta(eta), Z(Z), L(L), lamda(lamda), fai(fai) {};
	};
	struct PointCoordAz
	{
		PointGauss p1;
		PointGauss p2;
		Ellipsoid::Ellipsoid_para ell= Ellipsoid::bj54;
		double alphaCoord = 0;
		double A = 0;
		double sigmaA = 0;
		double sigmaAlphaCoord = 0;
		PointCoordAz() {};
		PointCoordAz(PointGauss p1, PointGauss p2, Ellipsoid::Ellipsoid_para ell)
			:p1(p1), p2(p2), ell(ell){};
	};
	
}
#endif