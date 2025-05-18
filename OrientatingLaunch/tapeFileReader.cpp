#pragma once
#include "tapeFileReader.h"
#include <codecvt>
#include <locale>


tapeFileReader::tapeFileReader() = default;
tapeFileReader::tapeFileReader(const std::string& filePath) {

	if (!readTapeFileNew(filePath)) {
		std::cerr << "读取标尺文件失败！" << std::endl;
	}
	else {
		std::cout << "读取标尺文件成功！" << std::endl;
	}



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

// 读取文件（新版本）
bool tapeFileReader::readTapeFileNew(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    // 读取第一行（发射点和基准点ID）
    std::string line;
    if (!std::getline(inFile, line)) {
        std::cerr << "Error: Empty file." << std::endl;
        return false;
    }
    auto firstLine = splitString(line, ',');
    if (firstLine.size() != 2) {
        std::cerr << "Error: First line must contain exactly two IDs." << std::endl;
        return false;
    }
    tapeFileData.FSId = firstLine[0];
    tapeFileData.JZId = firstLine[1];

    singleTape currentTape;
    bool processingTape = false;

    while (std::getline(inFile, line)) {
        // 跳过空行
        if (line.empty()) continue;

        // 检查标尺ID行（例如"标尺1"）
        if (line.find("tape") == 0) {
            if (processingTape) {
                // 保存当前标尺
                tapeFileData.allTapeData.push_back(currentTape);
                currentTape = singleTape();
            }
            currentTape.tapeId = line;
            processingTape = true;

            // 读取刻度ID行
            std::string markLine;
            while (std::getline(inFile, markLine)) {
                if (!markLine.empty()) break;
            }
            if (markLine.empty()) {
               std::cerr << "Error: Missing mark IDs for tape " << currentTape.tapeId << std::endl;
                return false;
            }
            currentTape.markId = splitString(markLine, ',');
            if (currentTape.markId.empty()) {
                std::cerr << "Error: No valid mark IDs for tape " << currentTape.tapeId << std::endl;
                return false;
            }
        }
        else {
            if (!processingTape) {
                std::cerr << "Error: Data outside of tape section: " << line << std::endl;
                return false;
            }
            // 处理距离观测数据
            auto parts = splitString(line, ',');
            if (parts.size() != 3) {
                std::cerr << "Error: Invalid observation format: " << line << std::endl;
                return false;
            }
            distObs obs;
            auto it1 = std::find(currentTape.markId.begin(), currentTape.markId.end(), parts[0]);
			auto it2 = std::find(currentTape.markId.begin(), currentTape.markId.end(), parts[1]);
			obs.from = it1 < it2 ? parts[0] : parts[1];
			obs.to = it1 < it2 ? parts[1] : parts[0];
            try {
                obs.distObsValue = std::stod(parts[2]);
            }
            catch (const std::exception& e) {
                std::cerr << "Error: Invalid distance value in line: " << line << std::endl;
                return false;
            }
            currentTape.distObsData.push_back(obs);
        }
    }

    // 保存最后一个标尺
    if (processingTape) {
        tapeFileData.allTapeData.push_back(currentTape);
    }

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

// 辅助函数：读取整个文件到UTF-8字符串
std::string tapeFileReader::readFileAsUTF8(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file");
    }

    // 读取二进制数据
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string contents = buffer.str();

    // 若文件含BOM头（Windows常见），跳过前3字节
    if (contents.size() >= 3 &&
        (uint8_t)contents[0] == 0xEF &&
        (uint8_t)contents[1] == 0xBB &&
        (uint8_t)contents[2] == 0xBF) {
        contents = contents.substr(3);
    }
    return contents;
}

