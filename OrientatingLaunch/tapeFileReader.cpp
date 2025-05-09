#pragma once
#include "tapeFileReader.h"


tapeFileReader::tapeFileReader() = default;
tapeFileReader::tapeFileReader(const std::string& filePath) {
    readTapeFile(filePath);
}
tapeFileReader::~tapeFileReader() = default;

// 读取文件
bool tapeFileReader::readTapeFile(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        return false;
    }

    std::string line;
    int groupStep = 0; // 0: 读取基准点, 1: 读取瞄准点, 2: 读取刻度值
    partTape currentTape;
    std::vector<std::string> mzdIds;
    std::vector<double> tapeValues;

    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty()) continue;

        switch (groupStep) {
        case 0: { // 读取基准点ID
            currentTape.JZId = line;
            groupStep = 1;
            break;
        }
        case 1: { // 读取瞄准点ID列表
            mzdIds = splitString(line, ',');
            groupStep = 2;
            break;
        }
        case 2: { // 读取刻度值
            std::vector<std::string> strValues = splitString(line, ',');
            for (const auto& str : strValues) {
                try {
                    tapeValues.push_back(std::stod(str));
                }
                catch (...) {
                    // 处理转换错误
                    file.close();
                    return false;
                }
            }

            // 检查数据一致性
            if (mzdIds.size() != tapeValues.size()) {
                file.close();
                return false;
            }

            // 填充数据到结构体
            for (size_t i = 0; i < mzdIds.size(); ++i) {
                mzdData data;
                data.MzdId = mzdIds[i];
                data.tapeValue = tapeValues[i];
                currentTape.MzdData.push_back(data);
            }

            // 添加到数据集合
            TapeData.push_back(currentTape);

            // 重置临时变量
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

// 辅助函数：去除字符串首尾空白
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

// 辅助函数：分割字符串
std::vector<std::string> tapeFileReader::splitString(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(trim(token));
    }
    return tokens;
}

