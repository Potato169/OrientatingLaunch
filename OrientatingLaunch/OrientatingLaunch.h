#pragma once
#ifndef ORIENTATINGLAUNCH_H
#define ORIENTATINGLAUNCH_H

#include "EdgeManager.h"
#include "In2FileReader.h"
#include "tapeFileReader.h"
#include <Eigen/Dense>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <utility>
#include <map>

struct AdjustmentResult {
    std::string FROM;
    std::string TO;
    std::string TYPE = "L";
    std::string VALUE;
    std::string M;
    std::string V;
    std::string RESULT;
    std::string Ri;
};

struct GeoAzimuth {
    std::string from;
    std::string to;
    double azimuth;
    std::string azimuth_dms;
};

struct StroAzimuth {
    std::string from;
    std::string to;
    double azimuth;
    std::string azimuth_dms;
};


class DebugLog {
public:
    ~DebugLog() { std::clog << buffer.str() << std::endl; } // 自动换行

    template<typename T>
    DebugLog& operator<<(T&& value) {
        buffer << std::forward<T>(value);
        return *this;
    }

private:
    std::ostringstream buffer; // 缓存输出内容
};

class LogWarning {
public:
    ~LogWarning() { std::cerr << buffer.str() << std::endl; } // 自动换行

    template<typename T>
    LogWarning& operator<<(T&& value) {
        buffer << std::forward<T>(value);
        return *this;
    }

private:
    std::ostringstream buffer;  // 缓存输出内容
};

class LogCritical {  // 专用错误日志类
public:
    ~LogCritical() { std::cerr << buffer.str() << std::endl; } // 自动换行到标准错误

    template<typename T>
    LogCritical& operator<<(T&& value) {
        buffer << std::forward<T>(value);  // 支持所有流式类型
        return *this;
    }

private:
    std::ostringstream buffer;  // 内容缓存
};

class OrientatingLaunch {
public:
    //OrientatingLaunch(EdgeManager& manager, const In2FileReader& reader);
	OrientatingLaunch(EdgeManager& manager, const In2FileReader& reader,
        tapeFileReader& tapeReader = tapeFileReader::getDefaultTapeReader());
    OrientatingLaunch(const OrientatingLaunch&) = delete;
    OrientatingLaunch& operator=(const OrientatingLaunch&) = delete;

    //   这里是射向标定的基本方法主要包括：
    //   1.初始化估计网中的每一条边的近似方位角
    //   2.构建系数矩阵、法方程、观测值向量、权阵等平差计算必要的数据
    //   3.进行平差计算与精度评定
	//   4.将结果按一定格式输出
	//   5.提供这些结果的访问接口
	//   6.当提供了标尺观测数据时，额外执行标尺观测数据的处理（可选）
    void initializeEdges();
    void buildMatrices();
    void performAdjustment();
    bool processTapeData();


    const Eigen::MatrixXd& getBMatrix() const { return B; }
    const Eigen::VectorXd& getLMatrix() const { return L; }
    const Eigen::MatrixXd& getNBBMatrix() const { return NBB; }
    const Eigen::VectorXd& getWMatrix() const { return W; }
    const std::vector<AdjustmentResult>& getResults() const { return results; }
    const std::vector<AdjustmentResult>& getVResults() const { return V_results; }
	const Eigen::VectorXd& getVFin() const { return VFin; }// 这里得到的V不包含和方程改正，且包含了没有参与平差的观测值改正，为最终结果
    const Eigen::VectorXd& getVNoSum() const { return VNoSum; }
    const Eigen::VectorXd& getV() const { return V; }
    const Eigen::VectorXd& getX() const { return x; }
    Eigen::MatrixXd getQx() const { return Qx; }
    double getSigma() const { return sigma; }
    Eigen::MatrixXd getDx() const { return Dx; }
    const std::unordered_map<std::string, int>& getEdgeColumnMap() const { return edgeColumnMap; }
    double getSigma0() const { return sigma0; }
    void printEdgeColumnMap();


	//  由于上面的平差结果为平面坐标方位角
    // 下面的模块用于将平面坐标方位角转换为天文方位角或大地方位角
	// 算法选择枚举类型
    // 方便在计算大地方位角的两种算法间切换，对比结果
    enum class GeodeticAlgorithm { Algorithm1, Algorithm2 };

	// 转换函数主体
    void convertAzimuths();

    // 读取平差坐标，存储在knownPointsCoord中
    bool readPointsFromFile(const std::string& filePath);


private:
    EdgeManager& manager;
    const In2FileReader& reader;

    tapeFileReader& tapeReader;

    std::unordered_set<std::string> processedStations;
    // 存储点的平差坐标
    std::unordered_map<std::string, std::pair<double, double>> pointsCoord;

    Eigen::MatrixXd B;
    Eigen::VectorXd L;
    std::vector<AdjustmentResult> results;
    std::vector<AdjustmentResult> V_results;
    Eigen::VectorXd x;
    Eigen::MatrixXd Qx;
    double sigma = 0.0;
    Eigen::MatrixXd Dx;
	Eigen::VectorXd V;// 这里的改正值包含和方程改正
    Eigen::VectorXd VNoSum;
	Eigen::VectorXd VFin;// 这里的改正值不包含和方程改正，且包含了没有参与平差的观测值改正，为最终结果
    Eigen::MatrixXd NBB;
    Eigen::VectorXd W; // W=BTPL
    std::vector<double> P; // 对角权阵P的向量形式
    // sigma0为先验单位权中误差
	double sigma0 = 0.0;
    std::unordered_map<std::string, int> edgeColumnMap;

    std::vector<int> sumEquationRows; // 用于记录和方程行的索引

	// 该模块用于记录不参与平差的观测值
    std::vector<int> filteredObsIndices; // 存储被过滤观测的原始索引
    std::vector<int> originalObsIndices; // 记录处理过的观测的原始索引
    std::vector<bool> isSumEquationRow; // 标记是否为和方程行
    int originalObsCounter = 0;     // 全局原始观测计数器

    // 大地以及天文结果向量
	std::vector<GeoAzimuth> geoAzimuthList;
	std::vector<StroAzimuth> stroAzimuthList;

    void processKnownPoints();
    void processObservations();
    void processEdgeQueue(std::queue<std::pair<std::string, std::string>>& edgeQueue,
        std::unordered_map<std::string, Ui::ObservationData>& unprocessedStations,
        std::queue<std::string>& stationQueue);
    void processStationQueue(std::queue<std::string>& stationQueue,
        std::queue<std::pair<std::string, std::string>>& edgeQueue,
        std::unordered_map<std::string, Ui::ObservationData>& unprocessedStations);

    bool processKnownAngel();
    static double calculateAzimuth(const Ui::KnownPoint& from, const Ui::KnownPoint& to);

    void buildEdgeColumnMapping();
    void processInitialObservation(int row, const std::string& stationId, const std::string& targetId);
    void processRegularObservation(int row, const std::string& stationId, const std::string& targetId, double obsValue);
    DirectedEdge* findInitialEdge(const std::string& stationId);

    void calculateCorrections();
    void updateEdgeAzimuths();
    void calculateAccuracyParameters();
    void generateResultList();

    // 字符串处理辅助函数
    static std::string trim(const std::string& s);
    static bool startsWith(const std::string& s, const std::string& prefix);
    static std::vector<std::string> split(const std::string& s, const std::string& delimiter);
    std::string doubleToString(const double value, int precision = 8);

	Eigen::VectorXd calVNoSum();// 这里得到的V剔除了和方程改正
    Eigen::VectorXd calVFin();// 这里在VNoSum的基础之上添加了没有参与平差的观测值改正


	// 平面坐标方位角转换天文方位角模块的辅助函数
	// 判断该边是否有对应天文信息
    bool hasAstronomicalInfo(const DirectedEdge* edge) const;
	// 平面坐标方位角转换大地方位角
    double convertToGeodeticAzimuth(double planeAzimuth) const;
	// 大地方位角转换天文方位角
    double convertToAstronomicalAzimuth(double geodeticAzimuth) const;

	// 标尺数据处理辅助函数
    // 初始化方位角
    bool processPartTapeValue(partTape& part);
    bool calPartTapeValue(partTape& part);

};

#endif // ORIENTATINGLAUNCH_H