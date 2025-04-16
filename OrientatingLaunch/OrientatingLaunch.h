#pragma once
#ifndef ORIENTATINGLAUNCH_H
#define ORIENTATINGLAUNCH_H

#include "EdgeManager.h"
#include "In2FileReader.h"
#include <Eigen/Dense>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <utility>

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
    OrientatingLaunch(EdgeManager& manager, const In2FileReader& reader);
    OrientatingLaunch(const OrientatingLaunch&) = delete;
    OrientatingLaunch& operator=(const OrientatingLaunch&) = delete;

    void initializeEdges();
    void buildMatrices();
    void performAdjustment();

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

private:
    EdgeManager& manager;
    const In2FileReader& reader;

    std::unordered_set<std::string> processedStations;

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

    // 字符串处理辅助函数
    static std::string trim(const std::string& s);
    static bool startsWith(const std::string& s, const std::string& prefix);
    static std::vector<std::string> split(const std::string& s, const std::string& delimiter);


	Eigen::VectorXd calVNoSum();// 这里得到的V剔除了和方程改正
    Eigen::VectorXd calVFin();// 这里在VNoSum的基础之上添加了没有参与平差的观测值改正

};

#endif // ORIENTATINGLAUNCH_H