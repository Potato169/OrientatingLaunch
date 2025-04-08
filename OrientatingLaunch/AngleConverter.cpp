#include "AngleConverter.h"
#include <iostream>
#include <algorithm>

using namespace std;

double AngleConverter::parseAngleString(const std::string& angleStr) {
    // 1. 提取符号并移除负号
    bool isNegative = startsWith(angleStr, "-");
    string cleanStr = isNegative ? angleStr.substr(1) : angleStr;

    // 2. 分割字符串
    vector<string> parts = split(cleanStr, '.');
    if (parts.size() != 2 || parts[1].length() != 6) {
        cerr << "Invalid format. Expected 6 decimal digits." << endl;
        return 0.0;
    }

    // 3. 解析度部分
    int degrees;
    try {
        degrees = stoi(parts[0]);
    }
    catch (...) {
        cerr << "Invalid degrees value: " << parts[0] << endl;
        return 0.0;
    }
    if (degrees < 0 || degrees >= 360) {
        cerr << "Degrees out of range: " << degrees << endl;
        return 0.0;
    }

    // 4. 解析分秒部分
    string decimalPart = parts[1];
    try {
        int minutes = stoi(decimalPart.substr(0, 2));
        int secondsInt = stoi(decimalPart.substr(2, 2));
        int secondsFrac = stoi(decimalPart.substr(4, 2));

        if (minutes < 0 || minutes >= 60) {
            cerr << "Invalid minutes value: " << minutes << endl;
            return 0.0;
        }
        if (secondsInt < 0 || secondsInt >= 60) {
            cerr << "Invalid seconds integer: " << secondsInt << endl;
            return 0.0;
        }
        if (secondsFrac < 0 || secondsFrac >= 100) {
            cerr << "Invalid seconds fraction: " << secondsFrac << endl;
            return 0.0;
        }

        // 5. 计算总和并应用符号
        double fraction = minutes / 60.0 + (secondsInt + secondsFrac / 100.0) / 3600.0;
        return (isNegative ? -1 : 1) * (degrees + fraction);
    }
    catch (...) {
        cerr << "Error parsing decimal part" << endl;
        return 0.0;
    }
}

string AngleConverter::formatAngleString(double degrees) {
    bool isNegative = degrees < 0.0;
    degrees = abs(degrees);

    int deg = static_cast<int>(degrees);
    double minDecimal = (degrees - deg) * 60;
    int min = static_cast<int>(minDecimal);
    double secTotal = (minDecimal - min) * 60;

    int secInt = static_cast<int>(secTotal);
    double secFrac = (secTotal - secInt) * 100;
    secFrac = round(secFrac);

    // 进位处理
    if (secFrac >= 100) {
        secFrac -= 100;
        secInt++;
    }
    if (secInt >= 60) {
        secInt -= 60;
        min++;
    }
    if (min >= 60) {
        min -= 60;
        deg++;
    }
    if (deg >= 360) {
        deg -= 360;
    }

    ostringstream oss;
    oss << (isNegative ? "-" : "")
        << setw(3) << setfill('0') << deg << "."
        << setw(2) << setfill('0') << min
        << setw(2) << setfill('0') << secInt
        << setw(2) << setfill('0') << static_cast<int>(secFrac);

    string result = oss.str();
    // 处理首位可能的负号和零的组合
    if (result.size() > 7 && result[3] == '.') {
        if (result[0] == '0') result.erase(0, 1);
        if (result[0] == '0') result.erase(0, 1);
    }
    return result;
}

vector<string> AngleConverter::split(const string& s, char delimiter) {
    vector<string> tokens;
    string token;
    istringstream tokenStream(s);
    while (getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

bool AngleConverter::startsWith(const string& str, const string& prefix) {
    return str.compare(0, prefix.length(), prefix) == 0;
}