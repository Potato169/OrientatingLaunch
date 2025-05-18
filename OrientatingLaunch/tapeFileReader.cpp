#pragma once
#include "tapeFileReader.h"
#include <codecvt>
#include <locale>


tapeFileReader::tapeFileReader() = default;
tapeFileReader::tapeFileReader(const std::string& filePath) {

	if (!readTapeFileNew(filePath)) {
		std::cerr << "��ȡ����ļ�ʧ�ܣ�" << std::endl;
	}
	else {
		std::cout << "��ȡ����ļ��ɹ���" << std::endl;
	}



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

// ��ȡ�ļ����°汾��
bool tapeFileReader::readTapeFileNew(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    // ��ȡ��һ�У������ͻ�׼��ID��
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
        // ��������
        if (line.empty()) continue;

        // �����ID�У�����"���1"��
        if (line.find("tape") == 0) {
            if (processingTape) {
                // ���浱ǰ���
                tapeFileData.allTapeData.push_back(currentTape);
                currentTape = singleTape();
            }
            currentTape.tapeId = line;
            processingTape = true;

            // ��ȡ�̶�ID��
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
            // �������۲�����
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

    // �������һ�����
    if (processingTape) {
        tapeFileData.allTapeData.push_back(currentTape);
    }

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

// ������������ȡ�����ļ���UTF-8�ַ���
std::string tapeFileReader::readFileAsUTF8(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file");
    }

    // ��ȡ����������
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string contents = buffer.str();

    // ���ļ���BOMͷ��Windows������������ǰ3�ֽ�
    if (contents.size() >= 3 &&
        (uint8_t)contents[0] == 0xEF &&
        (uint8_t)contents[1] == 0xBB &&
        (uint8_t)contents[2] == 0xBF) {
        contents = contents.substr(3);
    }
    return contents;
}

