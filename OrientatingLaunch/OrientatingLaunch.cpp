#pragma once
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
#include <limits>
#include <Eigen/Dense>
#include <string>
#include "GeodeticCalculator.h" 
#include "CoordSystem.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif 

#ifndef RAD_TO_ARCSEC
#define RAD_TO_ARCSEC 206264.806
#endif // !RAD_TO_ARCSEC


using namespace std;

// 定义极小值 epsilon
const double epsilon = std::numeric_limits<double>::epsilon() * 10; // 例如 10 倍 epsilon

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

string OrientatingLaunch::doubleToString(const double value, int precision) { // 默认保留8位小数
    ostringstream oss;
    oss << fixed << setprecision(precision) << value;
    return oss.str();
}



//OrientatingLaunch::OrientatingLaunch(EdgeManager& manager,
//    const In2FileReader& reader)
//    : manager(manager), reader(reader) {
//    sigma0 = reader.getAccuracyValues()[0].directionPrecision;
//}

OrientatingLaunch::OrientatingLaunch(EdgeManager& manager,
    const In2FileReader& reader, tapeFileReader& tapeReader)
    : manager(manager), reader(reader), tapeReader(tapeReader) {
    sigma0 = reader.getAccuracyValues()[0].directionPrecision;
    std::cout << "存在标尺数据，载入数据"
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

    //vector<Ui::KnownPoint> twdPoints, fwdPoints;
    //for (const auto& p : points) {
    //    string id = trim(p.id);
    //    if (startsWith(id, "TWD")) twdPoints.push_back(p);
    //    if (startsWith(id, "FWD")) fwdPoints.push_back(p);
    //}

    //if (twdPoints.size() != 2 || fwdPoints.size() != 2) {
    //    throw runtime_error("需要恰好2个TWD点和2个FWD点");
    //}

      // 构建点ID集合和映射
    std::unordered_set<std::string> point_ids;
    std::unordered_map<std::string, Ui::KnownPoint> point_map;
    for (const auto& point : points) {
        point_ids.insert(trim(point.id));
        point_map[trim(point.id)] = point;
    }

    // 遍历观测值
    for (const auto& stationData : reader.getObservations()) {
        string from_id = trim(stationData.stationId);
        for (const auto& obs : stationData.observations) {
            if (obs.type == "方向观测") {
                string to_id = trim(obs.pointId);
                // 检查两个端点是否存在
                if (point_ids.count(from_id) && point_ids.count(to_id)) {
                    const auto& from_point = point_map.at(from_id);
                    const auto& to_point = point_map.at(to_id);

                    // 计算坐标差
                    double dx = to_point.x - from_point.x;
                    double dy = to_point.y - from_point.y;

                    // 计算方位角（单位：弧度）
                    double azimuth_rad = std::atan2(dy, dx); // 注意参数顺序为dx, dy
                    double azimuth_deg = azimuth_rad * 180.0 / M_PI;

                    // 调整到0-360度范围
                    if (azimuth_deg < 0) {
                        azimuth_deg += 360.0;
                    }

                    // 添加边
                    auto knownEdge = manager.addEdge(from_id, to_id, azimuth_deg, true);
                }
            }
        }   
    }

    

    //double az1 = calculateAzimuth(twdPoints[0], fwdPoints[0]);
    //double az2 = calculateAzimuth(twdPoints[1], fwdPoints[1]);

    //manager.addEdge(trim(twdPoints[0].id), trim(fwdPoints[0].id), az1, true);
    //manager.addEdge(trim(twdPoints[1].id), trim(fwdPoints[1].id), az2, true);
}

void OrientatingLaunch::processObservations() {
    auto obsData = reader.getObservations();
    queue<string> stationQueue;
    unordered_map<string, Ui::ObservationData> unprocessedStations;

	// 遍历manager对象，查找已建立的边，这些边即作为起始边
    // 用来初始化edgeQueue进行后续处理
    queue<pair<string, string>> edgeQueue;
    for (auto it = manager.begin(); it != manager.end(); ++it) {
        string from = it->first.first;
        string to = it->first.second;
        edgeQueue.push({ from, to });
    }


    for (const auto& data : obsData) {
        unprocessedStations[data.stationId] = data;
    }

    while (!edgeQueue.empty() || !stationQueue.empty()) {
        processEdgeQueue(edgeQueue, unprocessedStations, stationQueue);
        processStationQueue(stationQueue, edgeQueue, unprocessedStations);
    }

	// 处理完所有测站后，检查是否有未处理的测站
    for (auto it = processedStations.begin(); it != processedStations.end(); ++it) {
		auto stationId = *it;
		if (unprocessedStations.find(stationId) != unprocessedStations.end()) {
			unprocessedStations.erase(stationId);
        }
        else {
			std::cerr << "方位角初始化异常！ " << endl;
        }
    }

    // 检查是否所有测站都已处理, 若存在未处理的测站，表明可能存在后方交会观测，需进一步结合近似坐标处理
    if (unprocessedStations.empty()) {
        std::cout << "所有测站方位角初始化均已处理完成！" << endl;
    }
    else {

        for (auto it = unprocessedStations.begin(); it != unprocessedStations.end(); ++it) {
            std::string unprocessedStationId = it->first;
            auto unprocessedStationObs = it->second;
            for (const auto& obs : unprocessedStationObs.observations) {
                if (obs.type == "方向观测") {
                    std::string unprocessedPointId = obs.pointId;
                    auto edge = manager.getEdge(unprocessedStationId, unprocessedPointId);
                    if (!edge) {
                        // 计算坐标差
                        double dx = pointsInfo[unprocessedPointId].x - pointsInfo[unprocessedStationId].x;
                        double dy = pointsInfo[unprocessedPointId].y - pointsInfo[unprocessedStationId].y;

                        // 计算方位角
                        double azimuth_rad = std::atan2(dy, dx);
                        double azimuth_deg = azimuth_rad * 180.0 / M_PI;

                        // 调整到0-360度范围
                        if (azimuth_deg < 0) {
                            azimuth_deg += 360.0;
                        }

                        // 添加边
                        auto knownEdge = manager.addEdge(unprocessedStationId, unprocessedPointId, azimuth_deg, false);

                    }
                }
            }
        }

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

        // 查找起始方向（观测值为0的L类型）(待改正)
        vector<Ui::Observation> validObs;
        int obsNum = 0;
        for (const auto& obs : data.observations) {
            if (obs.type == "方向观测" && obsNum == 0) {
                validObs.insert(validObs.begin(), obs); // 优先处理起始方向
                obsNum++;
            }
            else if (obs.type == "方向观测") {
                validObs.push_back(obs);
                obsNum++;
            }
            
        }

        if (validObs.empty()) continue;

        // 步骤1：查找存在的起始边
        DirectedEdge* startEdge = nullptr;
        double startObsValue = 0.0;
        std::string startToPoint;

        for (const auto& obs : validObs) {
            if (obs.type != "方向观测") continue; // 确保是方向观测
            std::string to = obs.pointId;
            DirectedEdge* edge = manager.getEdge(stationId, to);
            if (edge) {
                try {
                    startObsValue = AngleConverter::parseAngleString(obs.value);
                }
                catch (const std::exception& e) {
                    throw std::runtime_error("无效的观测值: " + obs.value);
                }
                startEdge = edge;
                startToPoint = to;
                break;
            }
        }

        if (!startEdge) {
            throw std::runtime_error("在validObs中未找到存在于EdgeManager的边");
        }

        // 步骤2：归化所有观测值
        std::vector<std::pair<std::string, double>> normalizedObs;
        for (const auto& obs : validObs) {
            if (obs.type != "方向观测") continue;
            string value;
            try {
                value = obs.value;
            }
            catch (const std::exception& e) {
                throw std::runtime_error("无效的观测值: " + obs.value);
            }
            // 计算相对于起始观测值的差值
			double normalized = AngleConverter::parseAngleString(value) - startObsValue;
            // 调整到0-360范围
            normalized = std::fmod(normalized, 360.0);
            if (normalized < 0) {
                normalized += 360.0;
            }
            normalizedObs.emplace_back(obs.pointId, normalized);
        }

        // 步骤3：获取起始边的实际方位角
        double startAzimuth = startEdge->azimuth;

        // 步骤4：计算并更新方位角
        for (const auto& normalizedPair : normalizedObs) {
            const auto& toPoint = normalizedPair.first;
            const auto& normalized = normalizedPair.second;
            double azimuth = startAzimuth + normalized;
            azimuth = std::fmod(azimuth, 360.0);
            if (azimuth < 0) {
                azimuth += 360.0;
            }

            // 获取或创建边
            DirectedEdge* edge = manager.getEdge(stationId, toPoint);
            if (!edge) {
                // 添加，设为未固定
                edge = manager.addEdge(stationId, toPoint, azimuth, false);
                edgeQueue.push({ stationId, toPoint });

            }
            else {
                // 仅更新未固定的边
                //if (!edge->isAzimuthFixed()) {
                //    edge->setAzimuth(azimuth);
                //}
                std::cout << "Warning: Edge " << stationId << "→" << toPoint
                    << " already exists " << std::endl;
            }

        }

    }
}

void OrientatingLaunch::buildMatrices() {

    sumEquationRows.clear(); // 清空之前的记录
	// 构造每个参数的列索引，每个方位角对应一个参数
    buildEdgeColumnMapping();

    // n为参数的个数
    const size_t n = edgeColumnMap.size();
    // 逐个观测值进行处理构造矩阵NBB与W
    int row = 0;
    // B、NBB、L与W矩阵的向量形式，最后使用eigen的矩阵存储
    vector<vector<double>> NBB_temp(n, vector<double>(n, 0.0));

    vector<double> W_temp(n, 0.0);

    vector<vector<double>> B_temp;

    vector<double> L_temp;

	vector<double> P_temp;

    //int m = 17;

 
    //L = Eigen::VectorXd::Zero(m);
    //Eigen::VectorXd P(m);

    //vector<pair<string, Ui::Observation>> directionObservations;
    //for (const auto& stationData : reader.getObservations()) {
    //    string stationId = trim(stationData.stationId);
    //    for (const auto& obs : stationData.observations) {
    //        if (obs.type == "方向观测") {
    //            directionObservations.emplace_back(stationId, obs);
    //        }
    //    }
    //}

    //for (int i = 0; i < m; ++i) {
    //    string stationId = directionObservations[i].first;
    //    const auto& obs = directionObservations[i].second;
    //    string targetId = trim(obs.pointId);
    //    double obsValue = AngleConverter::parseAngleString(obs.value);

    //    if (obsValue == 0.0) {
    //        processInitialObservation(i, stationId, targetId);
    //    }
    //    else {
    //        processRegularObservation(i, stationId, targetId, obsValue);
    //    }
    //}
    //const Eigen::IOFormat fmt(6, 0, ", ", "\n", "", "");
    //cout << "B:\n" << B.format(fmt) << endl;


    originalObsCounter = 0;
    // 遍历所有观测值
    // 下面这个循环处理每个测站
    for (const auto& stationData : reader.getObservations()) {
        string stationId = trim(stationData.stationId);

        // obsNum用来标记当前处理的观测值为本测站的第几个观测值
        // 通常第一个观测值的改正值等于对应同名边方位角的改正值
        // 因此需要特殊处理
        int obsNum = 0;
        //initAziKey用来记录当前测站的第一个方向观测值对应在edgeColumnMap中的位置
        int initEdgeKey = 0;

        // 记录当前测站的起始以及其他方向观测值
        DirectedEdge* initialEdge = nullptr;
        DirectedEdge* currentEdge = nullptr;

        // 当前测站的起始方向观测值
		double initObsValue = 0.0;

        // 定义当前测站和方程观测值的权值
        double p_sum = 0.0;
        // 定义当前观测值的和方程余项
        double L_sum = 0.0;
        // 定义当前测站的和方程系数
        vector<double> B_row_sum(n, 0.0);
        

        // 下面这个循环用来处理当前测站的所有观测值
        for (const auto& obs : stationData.observations) {
            if (obs.type == "方向观测") {
                obsNum++;
                row++;

                // 目标点ID去空格
                string targetId = trim(obs.pointId);

                // 按字典序生成边标识符
                string from, to;
                if (stationId < targetId) {
                    from = stationId;
                    to = targetId;
                }
                else {
                    from = targetId;
                    to = stationId;
                }
                string edgeKey = from + "→" + to;

                // 判断边的同名边方位角是否存在参数中
                if (edgeColumnMap.find(edgeKey) != edgeColumnMap.end()) {
                    originalObsIndices.push_back(originalObsCounter);
                    isSumEquationRow.push_back(false);

                    double p = 0.0;
                    // 当该观测值对应的边为已知边时默认观测精度为0.0001秒，几乎不存在误差
                    // 对应权重也会变大，相当于没有改正
                    if (manager.getEdge(trim(stationId), trim(obs.pointId))->isAzimuthFixed()) {   
                        p = (sigma0 * sigma0) /
                            (0.0001 * 0.0001);
                    }
                    // 当该观测值对应的边为未知边时根据读取的文件正常计算即可
                    else {
                        p = (sigma0 * sigma0) /
                            (stod(trim(obs.accuracyNumber)) * stod(trim(obs.accuracyNumber)));
                    }

					

                    if (p == 0) {
                        cerr << "矩阵构建失败，权值为0！" << endl;
                        return;
                    }
                    P_temp.push_back(p);
                    p_sum += 1.0 / p;

                    // 定义当前观测值的常数项
                    double L = 0.0;


                    // 定义B矩阵的第k行
                    vector<double> B_row(n, 0.0);


                    // 当观测值为该测站的第一个方向观测时，表示零方向
                    // 只有边的同名边方位角系数设为 1.0
                    if (obsNum == 1) {
                        // 保留起始观测值，便于后续计算
                        initialEdge = manager.getEdge(trim(stationId), trim(obs.pointId));
                        if (!initialEdge) {
                            cerr << "未找到起始边：" << stationId << endl;
                            return;
                        }
                        // 同时保留起始边的矩阵列号
                        initEdgeKey = edgeColumnMap[edgeKey];
                        // 同时保留起始方向观测值
                        initObsValue = stod(trim(obs.value));
                        // 计算当前误差方程的常数项，零方向常数项为0
                        L = 0.0;
                        L_sum += L;

                        // 检查 n 是否在有效范围内
                        if (initEdgeKey >= 0 && initEdgeKey < n) {
                            B_row[initEdgeKey] = 1.0;  // 将边的同名边方位角系数设为 1.0
                            B_row_sum[initEdgeKey] += 1.0;
                        }
                        else {
                            std::cout << "Error: build matrix element is out of bounds!" << std::endl;
                        }
                    }
                    // 当观测值为该测站的其他方向观测时，零方向的系数设为 -1.0，
                    // 其他边的同名边方位角系数设为 1.0
                    else {

                        currentEdge = manager.getEdge(trim(stationId), trim(obs.pointId));
                        // 计算方位角差
                        double initialAz = initialEdge->azimuth;
                        double currentAz = currentEdge->azimuth;
                        double delta = fmod(currentAz - initialAz + 360, 360);
						double delta2 = AngleConverter::parseAngleString(trim(obs.value));
                        L = (delta2 - delta);
                        L_sum += L;

                        // 填充B矩阵
                        string initialEdgeId = initialEdge->from < initialEdge->to ?
                            (ostringstream() << initialEdge->from << "→" << initialEdge->to).str() :
                            (ostringstream() << initialEdge->to << "→" << initialEdge->from).str();

                        string currentEdgeId = currentEdge->from < currentEdge->to ?
                            (ostringstream() << currentEdge->from << "→" << currentEdge->to).str() :
                            (ostringstream() << currentEdge->to << "→" << currentEdge->from).str();

                        if (edgeColumnMap.find(initialEdgeId) != edgeColumnMap.end()) {
                            B_row[initEdgeKey] = -1.0;  // 将起始边的同名边方位角系数设为 -1.0
                            B_row_sum[initEdgeKey] += -1.0;
                        }


                        if (edgeColumnMap.find(currentEdgeId) != edgeColumnMap.end()) {
                            B_row[edgeColumnMap[currentEdgeId]] = 1.0; // 将非起始边的同名边方位角系数设为 1.0
                            B_row_sum[edgeColumnMap[currentEdgeId]] += 1.0;
                        }
                    }

                    // 计算B的第k行对NBB的贡献
                    for (int i = 0; i < n; ++i) {
                        for (int j = 0; j < n; ++j) {
                            NBB_temp[i][j] += p * B_row[i] * B_row[j];
                        }
                    }

                    // 计算W的当前行贡献
                    double scalar = p * L;
                    for (int j = 0; j < n; ++j) {
                        W_temp[j] += scalar * B_row[j];
                    }

                    B_temp.push_back(B_row); //将该行系数添加到系数阵B中
					L_temp.push_back(L); //将该行常数项添加到常数项阵L中
                }
                else {
                    cerr << "build Matrix error！Please check the process when initialing edge" << endl;


      //              filteredObsIndices.push_back(originalObsCounter);
      //              if (obsNum == 1) {
      //                  // 保留起始观测值，便于后续计算
      //                  initialEdge = manager.getEdge(trim(stationId), trim(obs.pointId));
      //                  if (!initialEdge) {
      //                      cerr << "未找到起始边：" << stationId << endl;
      //                      return;
      //                  }

      //                  initEdgeKey = -1;
      //              }
      //              else {
						//// 这里还有一个逻辑未完善，因为可能存在某个测站的非第一个观测值
      //                  // 不在参数列表的情形，
      //              }
                }

                originalObsCounter++;

            }


        }

        if (fabs(p_sum) < epsilon) {
			cerr << "矩阵构建失败，权值为0！" << endl;
			return;
        }
        p_sum = 1.0 / p_sum;

        // 计算当前测站和方程对NBB的贡献
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                NBB_temp[i][j] += p_sum * B_row_sum[i] * B_row_sum[j];
            }
        }

        // 计算当前测站和方程对W的贡献
        double scalar = p_sum * L_sum;
        for (int j = 0; j < n; ++j) {
            W_temp[j] += scalar * B_row_sum[j];
        }

        P_temp.push_back(p_sum);
        B_temp.push_back(B_row_sum); //将该行和方程系数添加到矩阵B中
		L_temp.push_back(L_sum); //将该行和方程常数项添加到矩阵L中
        sumEquationRows.push_back(B_temp.size() - 1); // 记录和方程行的索引
        isSumEquationRow.push_back(true); // 标记为和方程行
    }

	// NBB矩阵赋值
    NBB = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            NBB(i, j) = NBB_temp[i][j];
        }
    }

    // W赋值
	W = Eigen::VectorXd::Zero(n);
	for (int i = 0; i < n; ++i) {
		W(i) = W_temp[i];
	}

    // B矩阵赋值
    const size_t m = B_temp.size();
    B = Eigen::MatrixXd::Zero(m, n);
    for (int i = 0; i < m; ++i) {
        if (B_temp[i].size() != n) {
            throw std::runtime_error("行维度不一致！");
        }
        for (int j = 0; j < n; ++j) {
            B(i, j) = B_temp[i][j];
        }
    }

	 //L赋值
	L = Eigen::VectorXd::Zero(m);
	for (int i = 0; i < m; ++i) {
		L(i) = L_temp[i];
	}


	Eigen::MatrixXd PMatrix = Eigen::Map<const Eigen::VectorXd>
        (P_temp.data(), P_temp.size()).asDiagonal();
    P = P_temp;
	Eigen::MatrixXd NBBTest = B.transpose() * PMatrix * B;
	Eigen::MatrixXd WTest = B.transpose() * PMatrix * L;



	// 打印矩阵
    const Eigen::IOFormat fnt(1, 0, ", ", "\n", "", "");
    //cout << "B:\n" << B.format(fnt) << endl;
    //cout << "P:\n" << PMatrix.format(fnt) << endl;
	const Eigen::IOFormat fmt(8, 0, ", ", "\n", "", "");	
	//cout << "L:\n" << L.format(fmt) << endl;	
	//cout << "W:\n" << W.format(fmt) << endl;
	//cout << "WTest:\n" << WTest.format(fmt) << endl;

 //   cout << "NBB:\n" << NBB.format(fmt) << endl;
	//cout << "NBBTest:\n" << NBBTest.format(fmt) << endl;
}

void OrientatingLaunch::buildEdgeColumnMapping() {
    set<pair<string, string>> processedEdges;
    edgeColumnMap.clear();

    for (auto it = manager.begin(); it != manager.end(); ++it) {
        string from = it->first.first;
        string to = it->first.second;
        pair<string, string> reversePair(to, from);

        // 当边为固定边时不参与平差
        //if (it->second->isAzimuthFixed()) continue;

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
    //int n = B.cols();
    //int m = B.rows();




    // 计算法方程矩阵 NBB = B^T * P * B (P暂时设为单位矩阵)
    //NBB = B.transpose() * B;

    const Eigen::IOFormat fmt(8, 0, ", ", "\n", "", "");


    Eigen::MatrixXd NBB_inv = NBB.inverse();

    // 计算改正数 x = NBB^{-1} * B^T * P * L
    x = NBB_inv * W;
    cout << "x:\n" << x.format(fmt) << endl;
    // 计算协因数矩阵 Qx = NBB^{-1}
    Qx = NBB_inv;
    //cout << "Qx:\n" << Qx.format(fmt) << endl;
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

    
    for (int i = 0; i < V.size(); ++i) {
        sigma += V(i) * P[i] * V(i);
    }

    // 计算单位权中误差
    sigma = std::sqrt(sigma / (m - n));

    // 计算协方差矩阵 Dx = Qx * sigma^2
    Dx = Qx * std::pow(sigma, 2);
}

void OrientatingLaunch::generateResultList() {
    results.clear();
	Eigen::VectorXd V1 = calVNoSum();
	Eigen::VectorXd V2 = calVFin();


    const Eigen::IOFormat fmt(8, 0, ", ", "\n", "", "");
    cout << "V:\n" << V.format(fmt) << endl;
    cout << "VNoSum:\n" << V1.format(fmt) << endl;

    cout << "VFin:\n" << V2.format(fmt) << endl;

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
            //res.M = 
            //    AngleConverter::formatAngleString(std::sqrt(Dx.diagonal()(edgeColumnMap[edgeId])));
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
                
                // V_results[obsIdx].V = AngleConverter::formatAngleString(V(obsIdx));

                // V_results[obsIdx].RESULT =
                //     AngleConverter::formatAngleString(AngleConverter::parseAngleString(obs.value)
                //                                                           + V(obsIdx));

                // 以度分秒形式展示
                res.V = AngleConverter::formatAngleString(VFin(obsIdx));
                res.VALUE = AngleConverter::formatAngleString(
                    AngleConverter::parseAngleString(obs.value));
                res.RESULT = AngleConverter::formatAngleString(
                    AngleConverter::parseAngleString(obs.value)+ VFin(obsIdx));

                // 以十进制展示   
                //res.V = doubleToString(VFin(obsIdx));
                //res.VALUE = doubleToString(AngleConverter::parseAngleString(obs.value));
                //res.RESULT = doubleToString(AngleConverter::parseAngleString(obs.value) + VFin(obsIdx));

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

    //若没有已知方位角则采用已知点反算已知方位角，跳出当前函数
    if (fwObservations.empty()) {
        cout << "没有已知方位角!从已知点反算方位角";
        return false; // 已知方位角数据为空直接返回
    }

	// 遍历已知方位角数据，向edgemanager中添加边
    for (int i = 0; i < fwObservations.size(); i++) {
        manager.addEdge(fwObservations[i].first,
            trim(fwObservations[i].second.pointId),
            AngleConverter::parseAngleString(fwObservations[i].second.value), true);
    }

    return true;
}

Eigen::VectorXd OrientatingLaunch::calVNoSum(){
    std::unordered_set<int> sum_rows_set(sumEquationRows.begin(), sumEquationRows.end());
    vector<int> keep_indices;
    keep_indices.reserve(V.size() - sumEquationRows.size());

    for (int i = 0; i < V.size(); ++i) {
        if (!sum_rows_set.count(i)) {
            keep_indices.push_back(i);
        }
    }

    VNoSum = Eigen::VectorXd::Zero(keep_indices.size());
    for (size_t j = 0; j < keep_indices.size(); ++j) {
        VNoSum[j] = V[keep_indices[j]];
    }

    return VNoSum;
}

Eigen::VectorXd OrientatingLaunch::calVFin(){
    // 初始化VFin为全零，长度等于所有原始方向观测数
    VFin = Eigen::VectorXd::Zero(originalObsCounter);

    // 遍历所有处理过的行，填充有效改正值
    int vIndex = 0; // V向量中的索引
    for (size_t i = 0; i < isSumEquationRow.size(); ++i) {
        if (!isSumEquationRow[i]) {
            // 普通行，从originalObsIndices获取原始索引
            int originalIdx = originalObsIndices[vIndex];
            VFin[originalIdx] = VNoSum[vIndex];
            vIndex++;
        }
        // 和方程行直接跳过
    }

    return VFin;
}


// 将平面坐标方位角转换为大地方位角以及天文方位角的模块
void OrientatingLaunch::convertAzimuths() {

    for (auto it = manager.begin(); it != manager.end(); ++it) {
        const auto& edgeKey = it->first;
        const std::string& from = edgeKey.first;
        const std::string& to = edgeKey.second;
        const auto reverseKey = std::make_pair(to, from);

        DirectedEdge* edge = it->second;
        double originalAzimuth = edge->azimuth;

        // 平面坐标方位角转大地方位角

        PointCoordAz coordAz;
		coordAz.p1.B = pointsInfo[from].BGeo;
		coordAz.p1.L = pointsInfo[from].LGeo;
		coordAz.p2.B = pointsInfo[to].BGeo;
		coordAz.p2.L = pointsInfo[to].LGeo;
        coordAz.alphaCoord = originalAzimuth;


        coordWGeoAngle(coordAz);
		geoAstroEdgeInfo[edgeKey].geoAzimuth = coordAz.A;
		geoAstroEdgeInfo[edgeKey].geoAzimuthDms = AngleConverter::formatAngleString(coordAz.A);

        //大地方位角转换为天文方位角
        PointAz pointAz;
        pointAz.pointName = from;
		pointAz.A = coordAz.A; // 大地方位角
		pointAz.xi = pointsInfo[from].BAstro - pointsInfo[from].BGeo; // 垂线偏差南北分量
        pointAz.eta = (pointsInfo[from].LAstro - pointsInfo[from].LGeo) * cos(pointsInfo[from].BAstro * M_PI / 180.0); // 垂线偏差东西分量
        pointAz.Z = 42.5268; // 天顶距 （用高程计算，高程文件获取）
        pointAz.L = pointsInfo[from].LGeo; // 大地经度 
        pointAz.lamda = pointsInfo[from].LAstro; // 地面天文经度 
        pointAz.fai = pointsInfo[from].BAstro;  // 地面天文纬度
		pointAz.sigmaAlpha = 1; // 方位角精度

        pointAz.H1 = 10; // 起点正常高（高程文件获取）
        pointAz.Zeta = 0; // 高程异常 (可以取值为0（忽略不计））
        pointAz.H2 = 21; // 终点正常高（高程文件获取）
        pointAz.sigmaH1 = 0.1;

        GeodeticCalculator::caclGeoAngleToAstro(pointAz);
		geoAstroEdgeInfo[edgeKey].astroAzimuth = pointAz.alpha_sea;
		geoAstroEdgeInfo[edgeKey].astroAzimuthDms = AngleConverter::formatAngleString(pointAz.alpha_sea);
    }


}


// 标尺数据处理模块
bool OrientatingLaunch::processTapeData() {
    

    for (partTape& part : tapeReader.TapeData) {
        if (!processPartTapeValue(part)) {
			cerr << "处理部分标尺数据失败！" << endl;
			return false;
        }
        if (!calPartTapeValue(part)) {
			cerr << "计算部分标尺数据失败！" << endl;
			return false;
        }
        
    }

    return true;
}

bool OrientatingLaunch::processPartTapeValue(partTape& part) {
    string jzdId = part.JZId;
	string firstMzdId = part.MzdData[0].MzdId;
	int finMzdIndex = part.MzdData.size() - 1;
	string finalMzdId = part.MzdData[finMzdIndex].MzdId;

	// 先将两个观测瞄准点的方位角初始化，使用平差之后的方位角，在manager对象中检索
	DirectedEdge* firstMzdEdge = manager.getEdge(jzdId, firstMzdId);
	DirectedEdge* finalMzdEdge = manager.getEdge(jzdId, finalMzdId);
    if (!firstMzdEdge || !finalMzdEdge) {
		throw std::invalid_argument("Invalid edge variation when initializing part tape");
		return false;
    }
    else {
		part.MzdData[0].fwValue = firstMzdEdge->azimuth;
		part.MzdData[finMzdIndex].fwValue = finalMzdEdge->azimuth;
    }

	// 计算起始于终止边的平面距离

    part.MzdData[0].distValue = sqrt(pow((pointsInfo[firstMzdId].x - pointsInfo[jzdId].x), 2)
        + pow((pointsInfo[firstMzdId].y - pointsInfo[jzdId].y), 2));
       
	part.MzdData[finMzdIndex].distValue = sqrt(pow((pointsInfo[finalMzdId].x - pointsInfo[jzdId].x), 2)
		+ pow((pointsInfo[finalMzdId].y - pointsInfo[jzdId].y), 2));

    // 用两个变量存
	double distJz2MzdFst = part.MzdData[0].distValue;
	double distJz2MzdFin = part.MzdData[finMzdIndex].distValue;


	// 计算从起始方向到终止方向的方向值，范围在0-360之间()
	double absDeltaFw = fmod(
		part.MzdData[finMzdIndex].fwValue -
		part.MzdData[0].fwValue + 360.0, 360.0);


    // 计算起始边到终止边的夹角，这里的角度的范围必然是0-180之间，存到变量vertexAngle
    if (absDeltaFw >= 360.0 || absDeltaFw < 0.0) {
		cerr << "方向值超出范围！" << endl;
    }
    else if (absDeltaFw > 180.0) {
        part.jzAngleByAzi = 360.0 - absDeltaFw;
	}else {
        part.jzAngleByAzi = absDeltaFw;
	}
    
    // 由余弦定理计算基准点顶点角
    double distMzdFst2MzdFin = part.MzdData[finMzdIndex].tapeValue - part.MzdData[0].tapeValue;
	// 这里的vertexAngleByCos仅仅是顶点角的余弦值，还需要进一步计算
	double jzAngleByCos = (pow(distJz2MzdFst, 2) + pow(distJz2MzdFin, 2) - pow(distMzdFst2MzdFin, 2))
        / (2 * distJz2MzdFst * distJz2MzdFin);
    // 判断输入合法性
    if (jzAngleByCos < -1.0 || jzAngleByCos > 1.0) {
        std::cerr << "Error: Input must be in [-1, 1]." << std::endl;
        return 1;
    }
	jzAngleByCos = acos(jzAngleByCos) * 180.0 / M_PI; // 弧度转化为角度
    // 计算闭合差，用于后续的平差
	double closingError = part.jzAngleByAzi - jzAngleByCos; // 闭合差
    part.MzdData[finMzdIndex].dirValue = jzAngleByCos + closingError; // 模拟平差的过程

	// 计算起始边对应顶点的夹角的余弦值
	double cosMzdFstAngle = (pow(distJz2MzdFst, 2) + pow(distMzdFst2MzdFin, 2) - pow(distJz2MzdFin, 2)) 
        / (2 * distJz2MzdFst * distMzdFst2MzdFin);
    // 判断结果合理性
    if (cosMzdFstAngle < -1.0 || cosMzdFstAngle > 1.0) {
        std::cerr << "Error: Input must be in [-1, 1]." << std::endl;
        return 1;
    }
    
    // 计算剩余瞄准点的方向值
    for (int curMzdIndex = 1; curMzdIndex < finMzdIndex; curMzdIndex++) {
        // 刻度边的长度
        double distMzdFst2MzdCur = part.MzdData[curMzdIndex].tapeValue - part.MzdData[0].tapeValue;

        // 计算基准点到当前瞄准点的距离
        double distJzd2MzdCur = pow(distJz2MzdFst, 2) + pow(distMzdFst2MzdCur, 2)
            - (2 * distJz2MzdFst * distMzdFst2MzdCur * cosMzdFstAngle);
        distJzd2MzdCur = sqrt(distJzd2MzdCur);
        part.MzdData[curMzdIndex].distValue = distJzd2MzdCur;

        // 计算当前瞄准点的方向值
		double MzdCurAngle = (pow(distJz2MzdFst, 2) + pow(distJzd2MzdCur, 2) - pow(distMzdFst2MzdCur, 2))
			/ (2 * distJz2MzdFst * distJzd2MzdCur);
		// 判断输入合法性
		if (MzdCurAngle < -1.0 || MzdCurAngle > 1.0) {
			std::cerr << "Error: Input must be in [-1, 1]." << std::endl;
			return 1;
		}
		MzdCurAngle = acos(MzdCurAngle) * 180.0 / M_PI; // 弧度转化为角度
		part.MzdData[curMzdIndex].dirValue = MzdCurAngle + closingError / 2.0;
    }

    return true;
}

bool OrientatingLaunch::calPartTapeValue(partTape& part) {



    return true;
}

// 读取已知点平面坐标文件,这里同时将平面坐标进行了转化，获得了初始的大地经纬度
bool OrientatingLaunch::readPlainPointsFromFile(const string& plainCoordFilePath) {
    ifstream file(plainCoordFilePath);

    if (!file.is_open()) {
        throw runtime_error("无法打开文件: " + plainCoordFilePath);
        return false;
    }

    string line;
    while (getline(file, line)) {
        vector<string> parts;
        string part;
        istringstream iss(line);

        // 按逗号分割并去除空格
        while (getline(iss, part, ',')) {
            parts.push_back(trim(part));
        }

        if (parts.size() != 3) continue; // 跳过格式错误行

        try {
            string  id = parts[0];
            double x = stod(parts[1]);
            double y = stod(parts[2]);
			pointsInfo[id].x = x;
			pointsInfo[id].y = y;

			// 这里把平面坐标转化为大地经纬度，
            // 默认使用bj54椭球参数，带号为19，带宽为6°
            // 这里的椭球参数、中央经线、NF采用了默认值，实际需要读取参数，注意后期更改
			PointGauss res = gaussBack(x, y);
			pointsInfo[id].BGeo = res.B;
			pointsInfo[id].LGeo = res.L;
			pointsInfo[id].BGeoDms = AngleConverter::formatAngleString(res.B);
			pointsInfo[id].LGeoDms = AngleConverter::formatAngleString(res.L);
            double& M = pointsInfo[id].M;
			double& N = pointsInfo[id].N;
			calMN(res.B, pointsInfo[id].ell, M, N);
		}
		catch (const invalid_argument&) {
			// 跳过转换失败的行
			continue;
        }
        catch (const exception&) {
            // 跳过转换失败的行
            continue;
        }
    }

    return true;
}

// 读取已知天文经纬度文件，同时将其他未知点的平面坐标进行归算得到所有点的天文经纬度
bool OrientatingLaunch::readAstroBLFromFile(const std::string& AstroBLFilePath) {
    ifstream file(AstroBLFilePath);

    if (!file.is_open()) {
        throw runtime_error("无法打开文件: " + AstroBLFilePath);
        return false;
    }

    string line;
    while (getline(file, line)) {
        vector<string> parts;
        string part;
        istringstream iss(line);

        // 按逗号分割并去除空格
        while (getline(iss, part, ',')) {
            parts.push_back(trim(part));
        }

        if (parts.size() != 3) continue; // 跳过格式错误行

        try {
            string id = parts[0];
            string BAs = parts[1];
            string LAs = parts[2];
            pointsInfo[id].BAstroDms = BAs;
            pointsInfo[id].LAstroDms = LAs;
            pointsInfo[id].BAstro = AngleConverter::parseAngleString(BAs);
            pointsInfo[id].LAstro = AngleConverter::parseAngleString(LAs);
        }
        catch (const invalid_argument&) {
            // 跳过转换失败的行
            continue;
        }
        catch (const exception&) {
            // 跳过转换失败的行
            continue;
        }
    }



    return true;

}

// 归心改正
bool OrientatingLaunch::correct4Centering() {
    bool ifFindRefStation = false;

    // 下面是用于归心转换的参考站的基本信息
    double fai0 = 0.0;
	double fai0_rad = 0.0;
	double lambda0 = 0.0;
	double lambda0_rad = 0.0;
	// 遍历所有点，找到第一个符合条件的参考站
	std::string refStationId;
    const pointCoord* refStation = nullptr;

    for (const auto& entry : pointsInfo) {
        const std::string& key = entry.first;
        const pointCoord& coord = entry.second;
        if (coord.BAstro != 361.0 || coord.LAstro != 361.0) {
            ifFindRefStation = true;
            fai0 = coord.BAstro;
			fai0_rad = fai0 * M_PI / 180.0;
			lambda0 = coord.LAstro;
			lambda0_rad = lambda0 * M_PI / 180.0;
			refStationId = key;
			refStation = &coord;
            break;
		}
    }
	if (!ifFindRefStation) {
		cerr << "没有找到归心转换的参考站，请检查已知天文经纬度文件的正确性！" << endl;
		return false;
	}

	double e = 0.0; // 两点间的距离
    double alpha = 0.0; // 两点大地方位角
    double alpha_rad = 0.0;
    for (auto& entry : pointsInfo) {
        const std::string& key = entry.first;
        pointCoord& coord = entry.second;
        if (coord.BAstro == 361.0 && coord.LAstro == 361.0) {
            PointGeoSolution pointGeoSolution;
			pointGeoSolution.B1 = refStation->BGeo;
			pointGeoSolution.L1 = refStation->LGeo;
			pointGeoSolution.B2 = coord.BGeo;
			pointGeoSolution.L2 = coord.LGeo;

            GeodeticCalculator::GeodeticSlutionBack_Bessel(pointGeoSolution, refStation->ell);
            e = pointGeoSolution.S;
            alpha = pointGeoSolution.A1;
			alpha_rad = alpha * M_PI / 180.0;

			coord.BAstro = fai0 + RAD_TO_ARCSEC * (e * cos(alpha_rad) / refStation->M)/3600;
			coord.LAstro = lambda0 + RAD_TO_ARCSEC * (e * sin(alpha_rad) / (refStation->N * cos(fai0_rad)))/3600;
			coord.BAstroDms = AngleConverter::formatAngleString(coord.BAstro);
			coord.LAstroDms = AngleConverter::formatAngleString(coord.LAstro);
        }
    }


	return true;
}

void OrientatingLaunch::geoWstro()
{
    //大地方位角与天文方位角的转换，函数接收PointAz结构体，
    // 该结构体内需要填充 alpha地面天文方位角、x（垂线偏差子午分量）i、eta（垂线偏差东西分量）、Z（天顶距）、B（大地纬度）、L（大地经度）、lamda（天文经度）、fai（地面天文纬度），
    // 如需计算海面天文方位角计算大地方位角，需要填充H1、H2（正常高）、Zeta（高程异常），如不填充，计算结果与使用地面天文方位角结果相同
    //caclAstroAngleByGeo将大地方位角转为天文方位角,PointAz结构体的地面大地方位角A和sigmaA（标准差）；海面大地方位角A_sea和sigmaA_sea（标准差）；海面天文方位角alpha_sea;海面天文纬度fai_sea
    // caclGeoAngleByAstro将天文方位角转为大地方位角将天文方位角转为大地方位角，PointAz结构体的天文方位角alpha和sigmaAlpha（标准差）将被填充
    //结构体里的经纬度、天顶距单位为度，垂线偏差及标准差的单位为秒，高度、高程异常单位为米
    PointAz pointAz;
    pointAz.pointName = "PointA";
    pointAz.alpha = 135.4685;
    pointAz.xi = 7.2;
    pointAz.eta = 3.4;
    pointAz.Z = 42.5268;
    pointAz.L = 102.7216;
    pointAz.lamda = 102.7221;
    pointAz.fai = 25.0396;
    pointAz.sigmaAlpha = 1;

    pointAz.H1 = 10;
    pointAz.Zeta = 0;
    pointAz.H2 = 21;
    pointAz.sigmaH1 = 0.1;
    GeodeticCalculator::caclAstroAngleToGeo(pointAz);
    cout << "Point Name: " << pointAz.pointName << endl;
    cout << "A: " << fixed << setprecision(10) << pointAz.A << " sigmaA: " << pointAz.sigmaA << endl;
    cout << endl;
    cout << "A_sea " << pointAz.A_sea << endl;
    cout << "fai_sea " << pointAz.fai_sea << endl;
    cout << "alpha_sea " << pointAz.alpha_sea << endl;

    pointAz.alpha = 0;
    GeodeticCalculator::caclGeoAngleToAstro(pointAz);
    cout << "alpha: " << fixed << setprecision(10) << pointAz.alpha << endl;
    cout << "sigmaAlpha: " << pointAz.sigmaAlpha << endl;

}

void OrientatingLaunch::coordWGeoAngle(CoordSystem::PointCoordAz &coordAz)
{
    //坐标方位角与大地方位角的转换，函数接收PointCoordAz结构体，
    // 该结构体内包含两个gauss结构p1(起点)\p2（终点），需要填充B\L\L0（中央子午线），结构体内存储的值单位为度，标准差单位为秒
    // 还包含大地方位角A和sigmaA（标准差），函数计算后向结构体内的坐标方位角alpha或大地方位角A进行填充
    // caclCoordAngleByGeo将大地方位角转为坐标方位角，caclGeoAngleByCoord将坐标方位角转为大地方位角
    // 默认椭球参数为bj54，若需要使用其他椭球参数，请修改结构体ell参数
    //coordAz.p1.B = dms2degrees("22.145661630");
    coordAz.p1.sigmaB = 0.1;
    //coordAz.p1.L = dms2degrees("119.251490021");
    coordAz.p1.L0 = dms2degrees("120");
    //coordAz.p2.B = dms2degrees("22.582945240");
    //coordAz.p2.L = dms2degrees("119.453851692");
    coordAz.p2.L0 = dms2degrees("120");
    //coordAz.A = dms2degrees("23.260222760");
    coordAz.sigmaA = 0.1;
    coordAz.ell = Ellipsoid::cgcs2000;


    //GeodeticCalculator::caclCoordAngleByGeo(coordAz);
    //cout << fixed << setprecision(15) << degrees2dms(coordAz.alphaCoord) << endl;
    //cout << fixed << setprecision(15) << coordAz.sigmaAlphaCoord << endl;

    GeodeticCalculator::caclGeoAngleByCoord(coordAz);
    cout << fixed << setprecision(15) << degrees2dms(coordAz.A) << endl;
    cout << fixed << setprecision(15) << coordAz.sigmaA << endl;


}

// 椭球参数：椭球长半轴、扁率

CoordSystem::PointGauss OrientatingLaunch::gaussBack(double x, double y, int sign, int width)
{
    PointGauss gauss;
    gauss.x = x;
    gauss.y = y;
	gauss.signeded = sign;
    //gauss.L0 = 111;
    gauss.width = width;// 6为6度带，3为3度带，0为自定应中央子午线，即指定L0进行计算
    gauss.sigmaB = 1;
    gauss.sigmaL = 1;
    //GeodeticCalculator::GaussForward(gauss, Ellipsoid::bj54);
    //cout << "6度带" << endl;
    //cout << setprecision(15) << gauss.x << "  " << gauss.y << "  " << gauss.L0 << "  " << degrees2dms(gauss.mca) << " " << degrees2dms(gauss.sigmaMCA) << "  " << gauss.signeded << endl;
    //cout << setprecision(15) << gauss.sigmaX << "  " << gauss.sigmaY << endl;

    cout << endl << "高斯反算结果：" << endl;
    // 实际得到椭球参数：长半轴、扁率
    GeodeticCalculator::GaussBack(gauss, Ellipsoid::bj54);
	cout << setprecision(15) << "  x：" << gauss.x << "  y：" << gauss.y  << "  " <<  " 带号： " << gauss.signeded << " 带宽： " << gauss.width << endl;
    cout << setprecision(15) << "  B：" << degrees2dms(gauss.B) << " L： " << degrees2dms(gauss.L) << "  " << endl;

	return gauss;
}


void OrientatingLaunch::calMN(const double& B, const CoordSystem::Ellipsoid::Ellipsoid_para& ell, double& M, double& N) {
    double a = ell.a;
    double f = ell.f;
    double b = a * (1 - f);
    double c = a * (a / b);
    double e1 = sqrt((a / b) * (a / b) - 1.0);
    double B1 = deg2rad(B);
    double n = e1 * cos(B1);
    double V = sqrt(1 + n * n);
    M = c / (V * V * V);
    N = c / V;
}