#include "AngleConversion.h"
namespace AngleConversion 
{
	double dms2degrees(string dms)
	{
		int dotPos = dms.find('.');
		string dStr = dms.substr(0, dotPos);
		int d = stoi(dStr);
		if (dotPos == -1)
		{
			return d;
		}
		string decimalPart = dms.substr(dotPos + 1);
		int m = 0; double s = 0;
		if (decimalPart.length() == 1)
			m = stoi(decimalPart.substr(0, 1)) * 10;
		else if (decimalPart.length() == 2)
			m = stoi(decimalPart.substr(0, 2));
		else if (decimalPart.length() >= 3)
		{
			m = stoi(decimalPart.substr(0, 2));
			s = stod("0." + decimalPart.substr(2)) * 100;
		}
		double degrees = 0;
		if (d > 0)
		{
			degrees = d + (m / 60.0) + (s / 3600.0);
		}
		else
		{
			degrees = d - (m / 60.0) - (s / 3600.0);
		}
		return degrees;
	}
	double degrees2dms(double degrees)
	{
		int d = degrees;
		int m = (degrees - d) * 60;
		double s = ((degrees - d) * 60 - m) * 60;
		double dms = d + m / 100.0 + s / 10000.0;
		return dms;
	}
	double deg2rad(double deg)
	{
		return deg * acos(-1) / 180;
	}
	double rad2deg(double rad)
	{
		return rad * 180 / acos(-1);
	}
}