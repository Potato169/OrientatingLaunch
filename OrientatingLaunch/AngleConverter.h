#pragma once
#ifndef ANGLECONVERTER_H
#define ANGLECONVERTER_H

#include <string>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <vector>
#include <stdexcept>

class AngleConverter {
public:
    // �ȷ���ת��Ϊʮ����
    static double parseAngleString(const std::string& angleStr);
    // ʮ����ת��Ϊ�ȷ���
    static std::string formatAngleString(double degrees, int decimalDigits = 4);

private:
    static std::vector<std::string> split(const std::string& s, char delimiter);
    static bool startsWith(const std::string& str, const std::string& prefix);
};

#endif // ANGLECONVERTER_H