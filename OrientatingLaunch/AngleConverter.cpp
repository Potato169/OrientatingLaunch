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
    if (parts.size() != 2 ) {
        cerr << "Invalid format. Expected degree and decimal parts\n" << endl;
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
    if (decimalPart.length() < 2) {
        cerr << "Invalid decimal part. At least two digits required for minutes." << endl;
        return 0.0;
    }
    try {
        // 提取分的部分
        int minutes = stoi(decimalPart.substr(0, 2));

        if (minutes < 0 || minutes >= 60) {
            cerr << "Invalid minutes value: " << minutes << endl;
            return 0.0;
        }

        // 提取秒的部分
        string secondsStr = decimalPart.substr(2);
        double totalSeconds = 0.0;

        if (!secondsStr.empty()) {
            if (secondsStr.length() < 2) {
                cerr << "Invalid seconds part: need at least two digits for integer seconds." << endl;
                return 0.0;
            }
            int secondsInt = stoi(secondsStr.substr(0, 2));
            if (secondsInt < 0 || secondsInt >= 60) {
                cerr << "Invalid seconds integer: " << secondsInt << endl;
                return 0.0;
            }

            string fracPart = secondsStr.substr(2);
            double frac = 0.0;
            if (!fracPart.empty()) {
                try {
                    frac = stod("0." + fracPart);
                }
                catch (...) {
                    cerr << "Invalid seconds fraction: " << fracPart << endl;
                    return 0.0;
                }
            }
            totalSeconds = secondsInt + frac;
        }


        // 5. 计算总和并应用符号
        double fraction = minutes / 60.0 + totalSeconds / 3600.0;
        return (isNegative ? -1 : 1) * (degrees + fraction);
    }
    catch (...) {
        cerr << "Error parsing decimal part" << endl;
        return 0.0;
    }
}

string AngleConverter::formatAngleString(double degrees, int decimalDigits) {
    bool isNegative = degrees < 0.0;
    degrees = abs(degrees);

    int deg = static_cast<int>(degrees);
    double remaining = degrees - deg;

    double minutes = remaining * 60;
    int min = static_cast<int>(minutes);
    double seconds = (minutes - min) * 60;

    double factor = pow(10, decimalDigits);
    double roundedSeconds = round(seconds * factor) / factor;

    int additionalMinutes = static_cast<int>(roundedSeconds) / 60;
    roundedSeconds = fmod(roundedSeconds, 60);
    min += additionalMinutes;

    int additionalDegrees = min / 60;
    min %= 60;
    deg += additionalDegrees;
    deg %= 360;

    int secInt = static_cast<int>(roundedSeconds);
    double secFrac = roundedSeconds - secInt;

    int secFracRounded = round(secFrac * factor);
    if (secFracRounded >= factor) {
        secFracRounded -= factor;
        secInt += 1;
        if (secInt >= 60) {
            secInt -= 60;
            min += 1;
            if (min >= 60) {
                min -= 60;
                deg += 1;
                deg %= 360;
            }
        }
    }

    ostringstream oss;
    oss << (isNegative ? "-" : "")
        << setw(3) << setfill('0') << deg << "."
        << setw(2) << setfill('0') << min
        << setw(2) << setfill('0') << secInt;

    if (decimalDigits > 0) {
        oss << setw(decimalDigits) << setfill('0') << secFracRounded;
    }

    return oss.str();
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