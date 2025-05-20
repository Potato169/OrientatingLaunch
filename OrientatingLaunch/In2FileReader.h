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
    // �洢����ֵ
    struct AccuracyValue {
        double directionPrecision = -1;
        double distancePrecision1 = -1;
        double distancePrecision2 = -1;
    };

    // �洢��֪������
    struct KnownPoint {
        std::string id; // ��ID
        double x;       // X����
        double y;       // Y����
    };

    // �洢�۲�ֵ
    struct Observation {
        std::string pointId;   // Ŀ���ID
        std::string type;      // �۲����ͣ�L: ����۲�, S: ����۲⣬A����λ�ǹ۲⣩
        std::string value;     // �۲�ֵ
        std::string accuracyNumber = "1";
    };

    // �洢�۲�վ����
    struct ObservationData {
        std::string stationId;              // ��վID
        std::vector<Observation> observations; // �۲�ֵ�б�
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

    // �ַ�������������
    static std::string trim(const std::string& s);
    static void replaceChineseComma(std::string& input);
    std::string removeUTF8BOM(const std::vector<char>& data);
    static std::vector<std::string> splitString(const std::string& str, char delimiter);
};

#endif // IN2FILEREADER_H