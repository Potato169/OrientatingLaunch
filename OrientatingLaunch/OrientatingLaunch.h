#pragma once
#ifndef ORIENTATINGLAUNCH_H
#define ORIENTATINGLAUNCH_H

#include "EdgeManager.h"
#include "In2FileReader.h"
#include "tapeFileReader.h"
#include "CoordSystem.h"
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

struct edgeAzimuth {
    double geoAzimuth = -1.0;
    std::string geoAzimuthDms;
	double astroAzimuth = -1.0;
	std::string astroAzimuthDms;
};

struct pointCoord {
	double x = 0.0;
	double y = 0.0;
	double BGeo = 0.0;
	double LGeo = 0.0;
    std::string BGeoDms;
	std::string LGeoDms;
	double BAstro = 361.0; // 这里默认值为361.0，表示未初始化
	double LAstro = 361.0;
	std::string BAstroDms;
    std::string LAstroDms;
	CoordSystem::Ellipsoid::Ellipsoid_para ell = CoordSystem::Ellipsoid::bj54;
	double NF = 0.0; 
    double M = 0.0; // 子午圈曲率半径
	double N = 0.0; // 卯酉圈曲率半径
	double L0 = 0.0; // 中央子午线

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


    // 字符串处理辅助函数
    static std::string trim(const std::string& s);
    static bool startsWith(const std::string& s, const std::string& prefix);
    static std::vector<std::string> split(const std::string& s, const std::string& delimiter);
    std::string doubleToString(const double value, int precision = 8);

	// 转换函数主体
    void convertAzimuths();

    // 读取平差坐标，存储在knownPointsCoord中，并对平面坐标处理的到大地经纬度
    bool readPlainPointsFromFile(const std::string& plainCoordFilePath);

    // 读取已有的天文经纬度，将其他的点的天文经纬度通过归心计算得到
    bool readAstroBLFromFile(const std::string& AstroBLFilePath);

    // 归心改正
    bool correct4Centering();

	// 大地与天文方位角转换的函数
    void geoWstro();

	// 平面方位角与大地方位角转换的函数
    void coordWGeoAngle(CoordSystem::PointCoordAz& coordAz);

    // 高斯反算(平面坐标到经纬度),参数分别对应平面坐标、带号以及带宽
    CoordSystem::PointGauss gaussBack(double x, double y, int sign = 19, int width = 6);

    // geoAstroEdgeInfo自定义哈希函数
    struct PairHash {
        template <typename T1, typename T2>
        std::size_t operator()(const std::pair<T1, T2>& p) const {
            // 分别计算两个元素的哈希值
            std::size_t hash1 = std::hash<T1>{}(p.first);
            std::size_t hash2 = std::hash<T2>{}(p.second);

            // 组合哈希值（避免简单异或导致的冲突）
            return hash1 ^ (hash2 << 1);
        }
    };

private:
    EdgeManager& manager;
    const In2FileReader& reader;

    tapeFileReader& tapeReader;

    std::unordered_set<std::string> processedStations;

	// 存储点的坐标（近似？）、大地经纬、天文经纬、NF、中央经线、卯酉圈、子午圈曲率半径、椭球参数
    std::unordered_map<std::string, pointCoord> pointsInfo;

	// 存储有向边的大地与天文方位角(包含标尺方位角)
	std::unordered_map<std::pair<std::string, std::string>, edgeAzimuth, PairHash> geoAstroEdgeInfo;

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



	Eigen::VectorXd calVNoSum();// 这里得到的V剔除了和方程改正
    Eigen::VectorXd calVFin();// 这里在VNoSum的基础之上添加了没有参与平差的观测值改正


	// 平面坐标方位角转换天文方位角模块的辅助函数
 //   bool hasAstronomicalInfo(const DirectedEdge* edge) const;
	//// 平面坐标方位角转换大地方位角
 //   double convertToGeodeticAzimuth(double planeAzimuth) const;
	//// 大地方位角转换天文方位角
 //   double convertToAstronomicalAzimuth(double geodeticAzimuth) const;

	// 标尺数据处理辅助函数
    // 初始化方位角
    bool processPartTapeValue(partTape& part);
    bool convertPartTapeValue(partTape& part);

    // 标尺距离测量值平差
    bool distObsAdjustment();
    // 标尺读数平差
    bool adjustSingleTape(singleTape& tape);
    // 标尺平差辅助函数
    // 辅助函数：获取刻度点索引
    size_t getMarkIndex(const std::vector<std::string>& marks, const std::string& id);
    // 将标尺的数据结构进行转化
    bool convertAllTapeToPartTape();
	// 将单条标尺数据的转化
    void convertSingleTape(const singleTape& st, const std::string& id,
        std::vector<partTape>& result);


    // 由椭球参数ell以及大地纬度B计算M、N（子午圈曲率半径以及卯酉圈曲率半径）
	void calMN(const double& B, const CoordSystem::Ellipsoid::Ellipsoid_para& ell, double& M, double& N);

};

#endif // ORIENTATINGLAUNCH_H