#pragma once
#include <string>
#include <vector>
#include <iostream>

// 单个瞄准点的数据（暂定）
struct mzdData {
	std::string MzdId; // 瞄准点ID（含两端实测与中间未测瞄准点名称）
	double tapeValue; // 瞄准点刻度值集合
	double dirValue = 0.0; // 方向值集合，需要后续处理计算得到，单位默认为弧度
	double distValue = 0.0; // 距离值集合，需要后续处理计算得到，单位默认为m
	double fwValue = 0.0; // 方位角值集合，需要后续处理计算得到
};

// 两个实测瞄准点间的观测数据（暂定）
struct partTape {
    std::string JZId; // 基准点ID
	std::vector<mzdData> MzdData; // 中间点观测数据集合
};



class tapeFileReader {
public:
	std::vector<partTape> TapeData;

    tapeFileReader();
    ~tapeFileReader();

    bool readTapeFile(const std::string& filePath);

    static tapeFileReader& getDefaultTapeReader() {
        static tapeFileReader defaultReader;

        partTape testTape;
		testTape.JZId = "JZD0";
		mzdData testMzd1, testMzd2, testMzd3, testMzd4, testMzd5, testMzd6;
		testMzd1.MzdId = "MZD1";
		testMzd1.tapeValue = 0.0;
		testMzd2.MzdId = "MZD1_MZD2_1";
		testMzd2.tapeValue = 1.0;
		testMzd3.MzdId = "MZD1_MZD2_2";
		testMzd3.tapeValue = 2.0;
		testMzd4.MzdId = "MZD1_MZD2_3";
		testMzd4.tapeValue = 3.0;
		testMzd5.MzdId = "MZD1_MZD2_4";
		testMzd5.tapeValue = 4.0;
		testMzd6.MzdId = "MZD2";
		testMzd6.tapeValue = 5.0;



        testTape.MzdData.push_back(testMzd1);
		testTape.MzdData.push_back(testMzd2);
		testTape.MzdData.push_back(testMzd3);
		testTape.MzdData.push_back(testMzd4);
		testTape.MzdData.push_back(testMzd5);
		testTape.MzdData.push_back(testMzd6);

		std::cout << "程序未读入标尺数据，使用测试标尺数据进行测试" << std::endl;

		defaultReader.TapeData.push_back(testTape);

        return defaultReader;
    };

private:
    //std::vector<singleTape> tapeData;
    


    std::string detectAndConvertEncoding(const std::vector<char>& data);
    std::string convertEncoding(const char* from, const char* to, const std::string& input);
    bool isValidUTF8(const std::string& str) const;
    static std::vector<std::string> splitString(const std::string& str, char delimiter);
    // 字符串处理辅助函数
    static std::string trim(const std::string& s);
};