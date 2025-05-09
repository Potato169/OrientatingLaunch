#pragma once
#include "tapeFileReader.h"


tapeFileReader::tapeFileReader() = default;
tapeFileReader::tapeFileReader(const std::string& filePath) {
    readTapeFile(filePath);
}
tapeFileReader::~tapeFileReader() = default;

// ��ȡ�ļ�
bool tapeFileReader::readTapeFile(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        return false;
    }

    std::string line;
    int groupStep = 0; // 0: ��ȡ��׼��, 1: ��ȡ��׼��, 2: ��ȡ�̶�ֵ
    partTape currentTape;
    std::vector<std::string> mzdIds;
    std::vector<double> tapeValues;

    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty()) continue;

        switch (groupStep) {
        case 0: { // ��ȡ��׼��ID
            currentTape.JZId = line;
            groupStep = 1;
            break;
        }
        case 1: { // ��ȡ��׼��ID�б�
            mzdIds = splitString(line, ',');
            groupStep = 2;
            break;
        }
        case 2: { // ��ȡ�̶�ֵ
            std::vector<std::string> strValues = splitString(line, ',');
            for (const auto& str : strValues) {
                try {
                    tapeValues.push_back(std::stod(str));
                }
                catch (...) {
                    // ����ת������
                    file.close();
                    return false;
                }
            }

            // �������һ����
            if (mzdIds.size() != tapeValues.size()) {
                file.close();
                return false;
            }

            // ������ݵ��ṹ��
            for (size_t i = 0; i < mzdIds.size(); ++i) {
                mzdData data;
                data.MzdId = mzdIds[i];
                data.tapeValue = tapeValues[i];
                currentTape.MzdData.push_back(data);
            }

            // ��ӵ����ݼ���
            TapeData.push_back(currentTape);

            // ������ʱ����
            currentTape = partTape();
            mzdIds.clear();
            tapeValues.clear();
            groupStep = 0;
            break;
        }
        }
    }

    file.close();
    return true;

}

// ����������ȥ���ַ�����β�հ�
std::string tapeFileReader::trim(const std::string& s) {
    auto start = s.begin();
    while (start != s.end() && std::isspace(*start)) {
        start++;
    }
    auto end = s.end();
    while (end != start && std::isspace(*(end - 1))) {
        end--;
    }
    return std::string(start, end);
}

// �����������ָ��ַ���
std::vector<std::string> tapeFileReader::splitString(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(trim(token));
    }
    return tokens;
}

