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
    static double parseAngleString(const std::string& angleStr);
    static std::string formatAngleString(double degrees);

private:
    static std::vector<std::string> split(const std::string& s, char delimiter);
    static bool startsWith(const std::string& str, const std::string& prefix);
};

#endif // ANGLECONVERTER_H