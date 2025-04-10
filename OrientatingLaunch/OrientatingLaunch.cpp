#include "OrientatingLaunch.h"
#include "AngleConverter.h"
#include <cmath>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <regex>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif 

using namespace std;

// 辅助函数实现
string OrientatingLaunch::trim(const string& s) {
    auto wsfront = find_if_not(s.begin(), s.end(), [](int c) {return isspace(c);});
    auto wsback = find_if_not(s.rbegin(), s.rend(), [](int c) {return isspace(c);}).base();
    return (wsback <= wsfront ? string() : string(wsfront, wsback));
}

bool OrientatingLaunch::startsWith(const string& s, const string& prefix) {
    return s.rfind(prefix, 0) == 0;
}

vector<string> OrientatingLaunch::split(const string& s, const string& delimiter) {
    vector<string> tokens;
    size_t pos = 0, last = 0;
    while ((pos = s.find(delimiter, last)) != string::npos) {
        tokens.push_back(s.substr(last, pos - last));
        last = pos + delimiter.length();
    }
    tokens.push_back(s.substr(last));
    return tokens;
}

double OrientatingLaunch::calculateAzimuth(const Ui::KnownPoint& from,
    const Ui::KnownPoint& to) {
    double dx = to.x - from.x;
    double dy = to.y - from.y;
    double azimuth = atan2(dy, dx) * 180.0 / M_PI;
    return fmod(azimuth + 360.0, 360.0);
}

OrientatingLaunch::OrientatingLaunch(EdgeManager& manager,
    const In2FileReader& reader)
    : manager(manager), reader(reader) {
}

void OrientatingLaunch::initializeEdges() {
    try {
        if (!processKnownAngel()) {
            processKnownPoints();
        }
        processObservations();
    }
    catch (const exception& e) {
        cerr << "Initialization failed: " << e.what() << endl;
        throw;
    }
}

void OrientatingLaunch::processKnownPoints() {
    auto points = reader.getKnownPoints();

    vector<Ui::KnownPoint> twdPoints, fwdPoints;
    for (const auto& p : points) {
        string id = trim(p.id);
        if (startsWith(id, "TWD")) twdPoints.push_back(p);
        if (startsWith(id, "FWD")) fwdPoints.push_back(p);
    }

    if (twdPoints.size() != 2 || fwdPoints.size() != 2) {
        throw runtime_error("需要恰好2个TWD点和2个FWD点");
    }

    double az1 = calculateAzimuth(twdPoints[0], fwdPoints[0]);
    double az2 = calculateAzimuth(twdPoints[1], fwdPoints[1]);

    manager.addEdge(trim(twdPoints[0].id), trim(fwdPoints[0].id), az1, true);
    manager.addEdge(trim(twdPoints[1].id), trim(fwdPoints[1].id), az2, true);
}

void OrientatingLaunch::processObservations() {
    auto obsData = reader.getObservations();
    queue<string> stationQueue;
    unordered_map<string, Ui::ObservationData> unprocessedStations;

    queue<pair<string, string>> edgeQueue;
    edgeQueue.push({ "TWD1", "FWD1" });
    edgeQueue.push({ "TWD2", "FWD2" });

    for (const auto& data : obsData) {
        unprocessedStations[data.stationId] = data;
    }

    while (!edgeQueue.empty() || !stationQueue.empty()) {
        processEdgeQueue(edgeQueue, unprocessedStations, stationQueue);
        processStationQueue(stationQueue, edgeQueue, unprocessedStations);
    }
}

void OrientatingLaunch::processEdgeQueue(
    queue<pair<string, string>>& edgeQueue,
    unordered_map<string, Ui::ObservationData>& unprocessedStations,
    queue<string>& stationQueue) {

    while (!edgeQueue.empty()) {
        auto edge = edgeQueue.front();
        edgeQueue.pop();
        auto from = edge.first;
		auto to = edge.second;
        auto* edgePtr = manager.getEdge(from, to);
        if (!edgePtr) continue;

        // 查找所有以from或to为测站且未处理的测站
        for (auto it = unprocessedStations.begin(); it != unprocessedStations.end(); ) {
            if ((it->first == from || it->first == to) &&
                processedStations.find(it->first) == processedStations.end()) {
                stationQueue.push(it->first);
                // it = unprocessedStations.erase(it);
                processedStations.insert(it->first);
            }
            else {
                ++it;
            }
        }
    }
}

void OrientatingLaunch::processStationQueue(
    queue<string>& stationQueue,
    queue<pair<string, string>>& edgeQueue,
    unordered_map<string, Ui::ObservationData>& unprocessedStations) {

    while (!stationQueue.empty()) {
        string stationId = stationQueue.front();
		stationQueue.pop();
        // Ui::ObservationData data = unprocessedStations.take(stationId);
        // 从原始数据获取（数据依然存在）
        if (unprocessedStations.find(stationId) == unprocessedStations.end()) continue;

        Ui::ObservationData data = unprocessedStations[stationId];

        // 查找起始方向（观测值为0的L类型）
        vector<Ui::Observation> validObs;
        for (const auto& obs : data.observations) {
            if (obs.type == "方向观测" && stod(obs.value) == 0.0) {
                validObs.insert(validObs.begin(), obs); // 优先处理起始方向
            }
            else if (obs.type == "方向观测") {
                validObs.push_back(obs);
            }
        }

        if (validObs.empty()) continue;

        double baseAzimuth = 0.0;
        double ObsAzimuth = 0.0;
        // 处理所有观测值
        for (const auto& obs : validObs) {
            if (stod(obs.value) == 0.0) { // 起始方向
                auto* edge = manager.getEdge(stationId, trim(obs.pointId));
                if (edge) {
                    baseAzimuth = edge->azimuth;
                }
                else {
                    edge = manager.getEdge(trim(obs.pointId), stationId);
                    if (edge) baseAzimuth = fmod(edge->azimuth + 180.0, 360.0);
                }
                // 创建新边
                auto* newEdge = manager.addEdge(stationId, trim(obs.pointId), baseAzimuth);
                if (newEdge) edgeQueue.push({ stationId, trim(obs.pointId) });

            }
            else { // 非起始方向
                //这里显然使用的是未转化前的角度，这里用自己给的函数转换一下，注意不要todouble
                // double delta = obs.value.toDouble();

                double delta = AngleConverter::parseAngleString(obs.value);
                ObsAzimuth = fmod(baseAzimuth + delta, 360.0);

                // 创建新边
                auto* newEdge = manager.addEdge(stationId, trim(obs.pointId), ObsAzimuth);
                if (newEdge) edgeQueue.push({ stationId, trim(obs.pointId) });

            }

        }
    }
}




void OrientatingLaunch::buildMatrices() {
    buildEdgeColumnMapping();

    vector<pair<string, Ui::Observation>> directionObservations;
    for (const auto& stationData : reader.getObservations()) {
        string stationId = trim(stationData.stationId);
        for (const auto& obs : stationData.observations) {
            if (obs.type == "方向观测") {
                directionObservations.emplace_back(stationId, obs);
            }
        }
    }

    int m = directionObservations.size();
    int n = edgeColumnMap.size();
    B = Eigen::MatrixXd::Zero(m, n);
    L = Eigen::VectorXd::Zero(m);
    Eigen::VectorXd P(m);



    for (int i = 0; i < m; ++i) {
        string stationId = directionObservations[i].first;
        const auto& obs = directionObservations[i].second;
        string targetId = trim(obs.pointId);
        double obsValue = AngleConverter::parseAngleString(obs.value);

        if (obsValue == 0.0) {
            processInitialObservation(i, stationId, targetId);
        }
        else {
            processRegularObservation(i, stationId, targetId, obsValue);
        }
    }
}

void OrientatingLaunch::buildEdgeColumnMapping() {
    set<pair<string, string>> processedEdges;
    edgeColumnMap.clear();

    for (auto it = manager.begin(); it != manager.end(); ++it) {
        string from = it->first.first;
        string to = it->first.second;
        pair<string, string> reversePair(to, from);

        if (processedEdges.find(reversePair) == processedEdges.end()) {
            ostringstream edgeIdStream;
            if (from < to) {
                edgeIdStream << from << "→" << to;
            }
            else {
                edgeIdStream << to << "→" << from;
            }
            string edgeId = edgeIdStream.str();
            edgeColumnMap[edgeId] = edgeColumnMap.size();
            processedEdges.insert(it->first);
        }
    }
}

void OrientatingLaunch::processInitialObservation(int row,
    const string& stationId, const string& targetId) {
    // 查找边或反向边
    DirectedEdge* edge = manager.getEdge(stationId, targetId);
    if (!edge) {
        edge = manager.getEdge(targetId, stationId);
    }

    // 获取标准化边ID
    std::ostringstream edgeIdStream;
    if (edge->from < edge->to) {
        edgeIdStream << edge->from << "→" << edge->to;
    }
    else {
        edgeIdStream << edge->to << "→" << edge->from;
    }
    string edgeId = edgeIdStream.str();

    // 填充B矩阵
    if (edgeColumnMap.find(edgeId) != edgeColumnMap.end()) {
        B(row, edgeColumnMap[edgeId]) = 1.0;
    }
    else {
        cerr << "未找到边列映射：" << edgeId << endl;
    }
}

void OrientatingLaunch::processRegularObservation(int row,
    const string& stationId, const string& targetId, double obsValue) {
    // 查找起始边（假设已处理过初始观测）
    DirectedEdge* initialEdge = findInitialEdge(stationId);
    if (!initialEdge) {
        cerr << "未找到起始边：" << stationId << endl;
        return;
    }

    // 查找当前观测边
    DirectedEdge* currentEdge = manager.getEdge(stationId, targetId);
    if (!currentEdge) {
        currentEdge = manager.getEdge(targetId, stationId);
    }

    // 计算方位角差
    double initialAz = initialEdge->azimuth;
    double currentAz = currentEdge->azimuth;
    double delta = fmod(currentAz - initialAz + 360, 360);
    L(row) = (obsValue - delta);

    // 填充B矩阵
    string initialEdgeId = initialEdge->from < initialEdge->to ?
        (ostringstream() << initialEdge->from << "→" << initialEdge->to).str() :
        (ostringstream() << initialEdge->to << "→" << initialEdge->from).str();

    string currentEdgeId = currentEdge->from < currentEdge->to ?
        (ostringstream() << currentEdge->from << "→" << currentEdge->to).str() :
        (ostringstream() << currentEdge->to << "→" << currentEdge->from).str();

    if (edgeColumnMap.find(initialEdgeId) != edgeColumnMap.end()) {
        B(row, edgeColumnMap[initialEdgeId]) = -1.0;
    }
    if (edgeColumnMap.find(currentEdgeId) != edgeColumnMap.end()) {
        B(row, edgeColumnMap[currentEdgeId]) = 1.0;
    }
}

DirectedEdge* OrientatingLaunch::findInitialEdge(const string& stationId) {
    string cleanedStationId = trim(stationId);
    const vector<Ui::ObservationData>& obsList = reader.getObservations();

    // 查找对应的测站观测数据
    for (const auto& obsData : obsList) {
        if (trim(obsData.stationId) == cleanedStationId) {
            // 遍历该测站的所有观测值
            for (const auto& obs : obsData.observations) {
                if (obs.type == "方向观测") {
                    double value = AngleConverter::parseAngleString(obs.value);
                    if (fabs(value) < std::numeric_limits<double>::epsilon()) {
                        string targetId = trim(obs.pointId);

                        // 优先查找正向边
                        DirectedEdge* edge = manager.getEdge(cleanedStationId, targetId);
                        if (!edge) {
                            edge = manager.getEdge(targetId, cleanedStationId);
                        }

                        if (edge) {
                            DebugLog() << "Intialized edge found: " << edge->from
                                << "→" << edge->to
                                << " (azimuth " << std::fixed << std::setprecision(4)
                                << edge->azimuth << "°)";
                        }
                        return edge; // 返回找到的第一个有效边
                    }
                }
            }
        }
    }

    LogWarning() << "No initialized edge found for station " << cleanedStationId;
    return nullptr;
}




//平差的整体流程
void OrientatingLaunch::performAdjustment() {
    try {
        buildMatrices();        // 先构建B、P和L矩阵
        calculateCorrections(); // 计算平差参数
        updateEdgeAzimuths();   // 更新边方位角
        // calculateObservationCorrections(); // 计算观测值改正
        calculateAccuracyParameters();     // 计算精度参数
        generateResultList();   // 生成结果列表
    }
    catch (const std::exception& e) {
        LogCritical() << "平差计算失败：" << e.what();
    }
}

void OrientatingLaunch::calculateCorrections() {
    int n = B.cols();
    int m = B.rows();


    // 计算法方程矩阵 NBB = B^T * P * B (P暂时设为单位矩阵)
    NBB = B.transpose() * B;

    const Eigen::IOFormat fmt(6, 0, ", ", "\n", "", "");
    cout << "B:\n" << B.format(fmt) << endl;
    cout << "NBB:\n" << NBB.format(fmt) << endl;

    Eigen::MatrixXd NBB_inv = NBB.inverse();

    // 计算改正数 x = NBB^{-1} * B^T * P * L
    x = NBB_inv * B.transpose() * L;
    cout << "x:\n" << x.format(fmt) << endl;
    // 计算协因数矩阵 Qx = NBB^{-1}
    Qx = NBB_inv;
    cout << "Qx:\n" << Qx.format(fmt) << endl;
    // 计算观测值改正数 V = Bx - L
    V = B * x - L;
    cout << "V:\n" << V.format(fmt) << endl;

}

void OrientatingLaunch::updateEdgeAzimuths() {

    for (auto it = edgeColumnMap.begin(); it != edgeColumnMap.end(); ++it) {
        string edgeId = it->first;
        vector<string> parts = split(edgeId, "→");
        string from = parts[0];
        string to = parts[1];
        int idx = it->second;
        // 使用成员变量x
        double correction = x(idx);

        // 获取并更新正向边
        DirectedEdge* edge = manager.getEdge(from, to);
        if (edge) {
            edge->setAzimuth(fmod(edge->azimuth + correction, 360.0));
        }
        else {
            // edge = manager.getEdge(to, from);
            // if (edge) {
            //     edge->setAzimuth(edge->azimuth - correction);
            // }
            throw std::invalid_argument("Invalid edge variation when updating EdgeAzimuths");
        }
    }
}

void OrientatingLaunch::calculateAccuracyParameters() {
    int m = V.size();
    int n = edgeColumnMap.size();

    // 计算单位权中误差
    sigma = std::sqrt(V.squaredNorm() / (m - n));

    // 计算协方差矩阵 Dx = Qx * sigma^2
    Dx = Qx * std::pow(sigma, 2);
}

void OrientatingLaunch::generateResultList() {
    results.clear();

    // 填充边结果
    for (const auto& edgeIdPair : edgeColumnMap) {
        const std::string& edgeId = edgeIdPair.first;
        vector<string> parts = split(edgeId, "→");
        string from = parts[0];
        string to = parts[1];

        DirectedEdge* edge = manager.getEdge(from, to);
        if (!edge) edge = manager.getEdge(to, from);

        if (edge) {
            AdjustmentResult res;
            res.FROM = edge->from;
            res.TO = edge->to;
            // res.VALUE = QString::number(edge->azimuth, 'f', 4);
            res.VALUE = 
                AngleConverter::formatAngleString(edge->azimuth);
            // res.M = QString::number(std::sqrt(Dx.diagonal()(edgeColumnMap[edgeId])), 'f', 6);
            res.M = 
                AngleConverter::formatAngleString(std::sqrt(Dx.diagonal()(edgeColumnMap[edgeId])));
            results.push_back(res);
        }
    }

    // 填充观测值结果
    int obsIdx = 0;
    for (const auto& stationData : reader.getObservations()) {
        for (const auto& obs : stationData.observations) {
            if (obs.type == "方向观测" && obsIdx < V.size()) {
                AdjustmentResult res;
                res.FROM = stationData.stationId;
                res.TO = obs.pointId;
                res.V = AngleConverter::formatAngleString(V(obsIdx));
                // V_results[obsIdx].V = AngleConverter::formatAngleString(V(obsIdx));

                // V_results[obsIdx].RESULT =
                //     AngleConverter::formatAngleString(AngleConverter::parseAngleString(obs.value)
                //                                                           + V(obsIdx));
                res.VALUE = AngleConverter::formatAngleString(
                    AngleConverter::parseAngleString(obs.value));
                res.RESULT = AngleConverter::formatAngleString(
                    AngleConverter::parseAngleString(obs.value)+ V(obsIdx));
                V_results.push_back(res);
                ++obsIdx;
            }
        }
    }
}

void OrientatingLaunch::printEdgeColumnMap() {
    cout << "Edge Column Mapping Table:" << endl;
    if (edgeColumnMap.empty()) {
        cout << "  [Empty]" << endl;
        return;
    }

    for (const auto& pair : edgeColumnMap) {
		const auto& key = pair.first;
        const auto& value = pair.second;
        cout << "  " << key << " → Column " << value << endl;
    }
}

bool OrientatingLaunch::processKnownAngel() {
    vector<Ui::ObservationData> obsData = reader.getObservations();

    vector<pair<string, Ui::Observation>> fwObservations;

    for (const auto& stationData : obsData) {
        string stationId = trim(stationData.stationId);
        for (const auto& obs : stationData.observations) {
            if (obs.type == "方位角") {
                fwObservations.push_back({ stationId, obs });
            }
        }

    }

    //若没有已知方位角则采用已知点反算已知方位角
    if (fwObservations.empty()) {
        cout << "没有已知方位角!从已知点反算方位角";
        return false; // 已知方位角数据为空直接返回
    }

    for (int i = 0; i < fwObservations.size(); i++) {
        manager.addEdge(fwObservations[i].first,
            trim(fwObservations[i].second.pointId),
            AngleConverter::parseAngleString(fwObservations[i].second.value), true);
    }

    return true;
}