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
    // 度分秒转化为十进制
    static double parseAngleString(const std::string& angleStr);
    // 十进制转化为度分秒
    static std::string formatAngleString(double degrees, int decimalDigits = 4);

private:
    static std::vector<std::string> split(const std::string& s, char delimiter);
    static bool startsWith(const std::string& str, const std::string& prefix);
};

#endif // ANGLECONVERTER_H