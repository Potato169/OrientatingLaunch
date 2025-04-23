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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif 

using namespace std;

// 定义极小值 epsilon）
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



OrientatingLaunch::OrientatingLaunch(EdgeManager& manager,
    const In2FileReader& reader)
    : manager(manager), reader(reader) {
    sigma0 = reader.getAccuracyValues()[0].directionPrecision;
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

	// 遍历manager对象，查找以建立的边，这些边即作为起始边
    // 用来初始化edgeQueue进行后续处理
    queue<pair<string, string>> edgeQueue;
    for (auto it = manager.begin(); it != manager.end(); ++it) {
        string from = it->first.first;
        string to = it->first.second;
        edgeQueue.push({ from, to });
    }

    //edgeQueue.push({ "TWD1", "FWD1" });
    //edgeQueue.push({ "TWD2", "FWD2" });

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

        //double baseAzimuth = 0.0;
        //double ObsAzimuth = 0.0;
        //// 处理所有观测值
        //for (const auto& obs : validObs) {
        //    if (stod(obs.value) == 0.0) {  
        //        //起始方向, 这里逻辑存在问题，如果
        //        //起始方向未提前在edge_map中定义，那么将初始化失败
        //        auto* edge = manager.getEdge(stationId, trim(obs.pointId));
        //        if (edge) {
        //            baseAzimuth = edge->azimuth;
        //        }
        //        else {
        //            edge = manager.getEdge(trim(obs.pointId), stationId);
        //            if (edge) {
        //                baseAzimuth = fmod(edge->azimuth + 180.0, 360.0);
        //            }
        //            else {

        //            }
        //        }
        //        // 创建新边
        //        auto* newEdge = manager.addEdge(stationId, trim(obs.pointId), baseAzimuth);
        //        if (newEdge) edgeQueue.push({ stationId, trim(obs.pointId) });

        //    }
        //    else { // 非起始方向
        //        //这里显然使用的是未转化前的角度，这里用自己给的函数转换一下，注意不要todouble
        //        // double delta = obs.value.toDouble();

        //        double delta = AngleConverter::parseAngleString(obs.value);
        //        ObsAzimuth = fmod(baseAzimuth + delta, 360.0);

        //        // 创建新边
        //        auto* newEdge = manager.addEdge(stationId, trim(obs.pointId), ObsAzimuth);
        //        if (newEdge) edgeQueue.push({ stationId, trim(obs.pointId) });

        //    }

        //}

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
    cout << "B:\n" << B.format(fnt) << endl;
    cout << "P:\n" << PMatrix.format(fnt) << endl;
	const Eigen::IOFormat fmt(8, 0, ", ", "\n", "", "");	
	cout << "L:\n" << L.format(fmt) << endl;	
	cout << "W:\n" << W.format(fmt) << endl;
	cout << "WTest:\n" << WTest.format(fmt) << endl;

    cout << "NBB:\n" << NBB.format(fmt) << endl;
	cout << "NBBTest:\n" << NBBTest.format(fmt) << endl;
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


void OrientatingLaunch::convertAzimuths(GeodeticAlgorithm algorithm) {
    std::set<std::pair<std::string, std::string>> processedEdges;

    for (auto it = manager.begin(); it != manager.end(); ++it) {
        const auto& edgeKey = it->first;
        const std::string& from = edgeKey.first;
        const std::string& to = edgeKey.second;
        const auto reverseKey = std::make_pair(to, from);

        // 跳过已处理的反向边
        if (processedEdges.find(reverseKey) != processedEdges.end()) continue;
        processedEdges.insert(edgeKey);

        DirectedEdge* edge = it->second;
        double originalAzimuth = edge->azimuth;

        // 平面坐标方位角转大地方位角
        double geodeticAzimuth = 0.0;
        switch (algorithm) {
        case GeodeticAlgorithm::Algorithm1:
            geodeticAzimuth = convertToGeodeticAzimuth(originalAzimuth, algorithm);
            break;
        case GeodeticAlgorithm::Algorithm2:
            geodeticAzimuth = convertToGeodeticAzimuth(originalAzimuth, algorithm);
            break;
        }

        // 天文方位角转换判断
        if (hasAstronomicalInfo(edge)) {
            double astronomicalAzimuth = convertToAstronomicalAzimuth(geodeticAzimuth);
            edge->setAzimuth(astronomicalAzimuth);
        }
        else {
            edge->setAzimuth(geodeticAzimuth);
        }
    }
}

// 辅助函数占位实现（待补充）
bool OrientatingLaunch::hasAstronomicalInfo(const DirectedEdge* edge) const {
    // 待补充天文方位角判断条件
    return false;
}

double OrientatingLaunch::convertToGeodeticAzimuth(double planeAzimuth, GeodeticAlgorithm algorithm) const {
    // 待补充具体转换算法
    return planeAzimuth;
}

double OrientatingLaunch::convertToAstronomicalAzimuth(double geodeticAzimuth) const {
    // 待补充天文方位角转换算法
    return geodeticAzimuth;
}