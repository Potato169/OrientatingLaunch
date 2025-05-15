#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
// 单个瞄准点的数据（暂定）
struct mzdData {
	std::string MzdId; // 瞄准点ID（含两端实测与中间未测瞄准点名称）
	double tapeValue; // 瞄准点刻度值集合
	double dirValue = -1.0; // 方向值集合，需要后续处理计算得到，单位默认为角度
	// 注意！这里的方向值不同于常规的方向值，虽然是从起始方向到终止方向的夹角，
	// 但是范围在0-180度之间，逆时针或顺时针取决于两条边的空间关系
	double distValue = -1.0; // 距离值集合，需要后续处理计算得到，单位默认为m
	double fwValue = 361.0; // 方位角值集合，需要后续处理计算得到（坐标方位角, 单位度）
	std::string fwValueDms; // 方位角值集合（度分秒格式）
	double reverseFwValue = 361.0; // 反方位角值集合，需要后续处理计算得到（坐标方位角）
	std::string reverseFwValueDms; // 反方位角值集合（度分秒格式）
	double fwValueGeo = 361.0; // 大地方位角值集合，需要后续处理计算得到
	std::string fwValueGeoDms; // 大地方位角值集合（度分秒格式）
	double reverseFwValueGeo = 361.0; // 反向大地方位角值集合，需要后续处理计算得到
	std::string reverseFwValueGeoDms; // 反向大地方位角值集合（度分秒格式）
	double fwValueAstro = 361.0; // 天文方位角值集合，需要后续处理计算得到
	std::string fwValueAstroDms; // 天文方位角值集合（度分秒格式）
	double reverseFwValueAstro = 361.0; // 反向天文方位角值集合，需要后续处理计算得到
	std::string reverseFwValueAstroDms; // 反向天文方位角值集合（度分秒格式）
};

// 两个实测瞄准点间的观测数据（暂定）
struct partTape {
    std::string JZId; // 基准点ID
	std::vector<mzdData> MzdData; // 中间点观测数据集合
	double jzAngleByAzi = 0.0; // 由方位角计算的起始边到终止边的夹角 单位为度 范围0-180
};

struct allTape {
	std::string JZId; // 基准点ID
	std::string FSId; // 发射点ID
	std::vector<singleTape> allTapeData; // 所有标尺数据集合
};

struct singleTape {
	std::string tapeId; // 标尺ID
	std::vector<std::string> markId; // 标尺刻度ID集合
	std::vector<double> markValue ; // 标尺刻度值集合(平差后确定)
	std::vector<distObs> distObsData; // 标尺观测值集合
};

struct distObs {
	std::string from; // 起始点ID
	std::string to; // 终止点ID
	double distObaValue; // 距离值
};

class tapeFileReader {
public:
	std::vector<partTape> TapeData;

	tapeFileReader();
    tapeFileReader(const std::string& filePath);
    ~tapeFileReader();

	// 读取标尺文件
	bool readTapeFile(const std::string& filePath);

	// 读取标尺文件（新版本）
	bool readTapeFileNew(const std::string& filePath);

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
    static std::vector<std::string> splitString(const std::string& str, char delimiter);
    // 字符串处理辅助函数
    static std::string trim(const std::string& s);

	// 检查是否为标尺行
	bool isTapeLine(const std::string& s) {
		return s.find("标尺") == 0;
	}
};