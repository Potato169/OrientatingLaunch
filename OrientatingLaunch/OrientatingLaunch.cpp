#include "OrientatingLaunch.h"
#include "AngleConverter.h"
#include <cmath>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <regex>

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

double OrientatingLaunch::calculateAzimuth(const Ui::KnownPoint& from, const Ui::KnownPoint& to) {
    double dx = to.x - from.x;
    double dy = to.y - from.y;
    double azimuth = atan2(dy, dx) * 180.0 / M_PI;
    return fmod(azimuth + 360.0, 360.0);
}

OrientatingLaunch::OrientatingLaunch(EdgeManager& manager, const In2FileReader& reader)
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

        if (validObs.isEmpty()) continue;

        double baseAzimuth = 0.0;
        double ObsAzimuth = 0.0;
        // 处理所有观测值
        for (const auto& obs : validObs) {
            if (obs.value.toDouble() == 0.0) { // 起始方向
                auto* edge = manager.getEdge(stationId, obs.pointId.trimmed());
                if (edge) {
                    baseAzimuth = edge->azimuth;
                }
                else {
                    edge = manager.getEdge(obs.pointId.trimmed(), stationId);
                    if (edge) baseAzimuth = fmod(edge->azimuth + 180.0, 360.0);
                }
                // 创建新边
                auto* newEdge = manager.addEdge(stationId, obs.pointId.trimmed(), baseAzimuth);
                if (newEdge) edgeQueue.enqueue({ stationId, obs.pointId.trimmed() });

            }
            else { // 非起始方向
                //这里显然使用的是未转化前的角度，这里用自己给的函数转换一下，注意不要todouble
                // double delta = obs.value.toDouble();

                double delta = AngleConverter::parseAngleString(obs.value);
                ObsAzimuth = fmod(baseAzimuth + delta, 360.0);

                // 创建新边
                auto* newEdge = manager.addEdge(stationId, obs.pointId.trimmed(), ObsAzimuth);
                if (newEdge) edgeQueue.enqueue({ stationId, obs.pointId.trimmed() });

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

// 其他方法的实现需要类似地转换Qt部分到标准库...

void OrientatingLaunch::printEdgeColumnMap() {
    cout << "Edge Column Mapping Table:" << endl;
    if (edgeColumnMap.empty()) {
        cout << "  [Empty]" << endl;
        return;
    }

    for (const auto& [key, value] : edgeColumnMap) {
        cout << "  " << key << " → Column " << value << endl;
    }
}

// 注意：需要实现所有其他成员函数的转换，特别是涉及字符串处理和容器操作的部分