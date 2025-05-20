#pragma once


#ifndef IN2FILEREADER_H
#define IN2FILEREADER_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>
//#include "iconv.h"
//#pragma comment(lib,"libiconv.lib")

namespace Ui {
    // 存储精度值
    struct AccuracyValue {
        double directionPrecision = -1;
        double distancePrecision1 = -1;
        double distancePrecision2 = -1;
    };

    // 存储已知点坐标
    struct KnownPoint {
        std::string id; // 点ID
        double x;       // X坐标
        double y;       // Y坐标
    };

    // 存储观测值
    struct Observation {
        std::string pointId;   // 目标点ID
        std::string type;      // 观测类型（L: 方向观测, S: 距离观测，A：方位角观测）
        std::string value;     // 观测值
        std::string accuracyNumber = "1";
    };

    // 存储观测站数据
    struct ObservationData {
        std::string stationId;              // 测站ID
        std::vector<Observation> observations; // 观测值列表
    };
}

class In2FileReader {
public:
    In2FileReader();
    ~In2FileReader();

    bool readIn2File(const std::string& filePath);

    std::vector<Ui::AccuracyValue> getAccuracyValues() const { return accuracyValues; }
    std::vector<Ui::KnownPoint> getKnownPoints() const { return knownPoints; }
    std::vector<Ui::ObservationData> getObservations() const { return observations; }

private:
    std::vector<Ui::AccuracyValue> accuracyValues;
    std::vector<Ui::KnownPoint> knownPoints;
    std::vector<Ui::ObservationData> observations;

    //std::string detectAndConvertEncoding(const std::vector<char>& data);
    //std::string convertEncoding(const char* from, const char* to, const std::string& input);
    //bool isValidUTF8(const std::string& str) const;

    // 字符串处理辅助函数
    static std::string trim(const std::string& s);
    static void replaceChineseComma(std::string& input);
    std::string removeUTF8BOM(const std::vector<char>& data);
    static std::vector<std::string> splitString(const std::string& str, char delimiter);
};

#endif // IN2FILEREADER_H