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

// ���弫Сֵ epsilon
const double epsilon = std::numeric_limits<double>::epsilon() * 10; // ���� 10 �� epsilon

// ��������ʵ��
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

string OrientatingLaunch::doubleToString(const double value, int precision) { // Ĭ�ϱ���8λС��
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
    std::cout << "���ڱ�����ݣ���������"
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
    //    throw runtime_error("��Ҫǡ��2��TWD���2��FWD��");
    //}

      // ������ID���Ϻ�ӳ��
    std::unordered_set<std::string> point_ids;
    std::unordered_map<std::string, Ui::KnownPoint> point_map;
    for (const auto& point : points) {
        point_ids.insert(trim(point.id));
        point_map[trim(point.id)] = point;
    }

    // �����۲�ֵ
    for (const auto& stationData : reader.getObservations()) {
        string from_id = trim(stationData.stationId);
        for (const auto& obs : stationData.observations) {
            if (obs.type == "����۲�") {
                string to_id = trim(obs.pointId);
                // ��������˵��Ƿ����
                if (point_ids.count(from_id) && point_ids.count(to_id)) {
                    const auto& from_point = point_map.at(from_id);
                    const auto& to_point = point_map.at(to_id);

                    // ���������
                    double dx = to_point.x - from_point.x;
                    double dy = to_point.y - from_point.y;

                    // ���㷽λ�ǣ���λ�����ȣ�
                    double azimuth_rad = std::atan2(dy, dx); // ע�����˳��Ϊdx, dy
                    double azimuth_deg = azimuth_rad * 180.0 / M_PI;

                    // ������0-360�ȷ�Χ
                    if (azimuth_deg < 0) {
                        azimuth_deg += 360.0;
                    }

                    // ��ӱ�
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

	// ����manager���󣬲����ѽ����ıߣ���Щ�߼���Ϊ��ʼ��
    // ������ʼ��edgeQueue���к�������
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

	// ���������в�վ�󣬼���Ƿ���δ����Ĳ�վ
    for (auto it = processedStations.begin(); it != processedStations.end(); ++it) {
		auto stationId = *it;
		if (unprocessedStations.find(stationId) != unprocessedStations.end()) {
			unprocessedStations.erase(stationId);
        }
        else {
			std::cerr << "��λ�ǳ�ʼ���쳣�� " << endl;
        }
    }

    // ����Ƿ����в�վ���Ѵ���, ������δ����Ĳ�վ���������ܴ��ں󷽽���۲⣬���һ����Ͻ������괦��
    if (unprocessedStations.empty()) {
        std::cout << "���в�վ��λ�ǳ�ʼ�����Ѵ�����ɣ�" << endl;
    }
    else {

        for (auto it = unprocessedStations.begin(); it != unprocessedStations.end(); ++it) {
            std::string unprocessedStationId = it->first;
            auto unprocessedStationObs = it->second;
            for (const auto& obs : unprocessedStationObs.observations) {
                if (obs.type == "����۲�") {
                    std::string unprocessedPointId = obs.pointId;
                    auto edge = manager.getEdge(unprocessedStationId, unprocessedPointId);
                    if (!edge) {
                        // ���������
                        double dx = pointsInfo[unprocessedPointId].x - pointsInfo[unprocessedStationId].x;
                        double dy = pointsInfo[unprocessedPointId].y - pointsInfo[unprocessedStationId].y;

                        // ���㷽λ��
                        double azimuth_rad = std::atan2(dy, dx);
                        double azimuth_deg = azimuth_rad * 180.0 / M_PI;

                        // ������0-360�ȷ�Χ
                        if (azimuth_deg < 0) {
                            azimuth_deg += 360.0;
                        }

                        // ��ӱ�
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

        // ����������from��toΪ��վ��δ����Ĳ�վ
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
        // ��ԭʼ���ݻ�ȡ��������Ȼ���ڣ�
        if (unprocessedStations.find(stationId) == unprocessedStations.end()) continue;

        Ui::ObservationData data = unprocessedStations[stationId];

        // ������ʼ���򣨹۲�ֵΪ0��L���ͣ�(������)
        vector<Ui::Observation> validObs;
        int obsNum = 0;
        for (const auto& obs : data.observations) {
            if (obs.type == "����۲�" && obsNum == 0) {
                validObs.insert(validObs.begin(), obs); // ���ȴ�����ʼ����
                obsNum++;
            }
            else if (obs.type == "����۲�") {
                validObs.push_back(obs);
                obsNum++;
            }
            
        }

        if (validObs.empty()) continue;

        // ����1�����Ҵ��ڵ���ʼ��
        DirectedEdge* startEdge = nullptr;
        double startObsValue = 0.0;
        std::string startToPoint;

        for (const auto& obs : validObs) {
            if (obs.type != "����۲�") continue; // ȷ���Ƿ���۲�
            std::string to = obs.pointId;
            DirectedEdge* edge = manager.getEdge(stationId, to);
            if (edge) {
                try {
                    startObsValue = AngleConverter::parseAngleString(obs.value);
                }
                catch (const std::exception& e) {
                    throw std::runtime_error("��Ч�Ĺ۲�ֵ: " + obs.value);
                }
                startEdge = edge;
                startToPoint = to;
                break;
            }
        }

        if (!startEdge) {
            throw std::runtime_error("��validObs��δ�ҵ�������EdgeManager�ı�");
        }

        // ����2���黯���й۲�ֵ
        std::vector<std::pair<std::string, double>> normalizedObs;
        for (const auto& obs : validObs) {
            if (obs.type != "����۲�") continue;
            string value;
            try {
                value = obs.value;
            }
            catch (const std::exception& e) {
                throw std::runtime_error("��Ч�Ĺ۲�ֵ: " + obs.value);
            }
            // �����������ʼ�۲�ֵ�Ĳ�ֵ
			double normalized = AngleConverter::parseAngleString(value) - startObsValue;
            // ������0-360��Χ
            normalized = std::fmod(normalized, 360.0);
            if (normalized < 0) {
                normalized += 360.0;
            }
            normalizedObs.emplace_back(obs.pointId, normalized);
        }

        // ����3����ȡ��ʼ�ߵ�ʵ�ʷ�λ��
        double startAzimuth = startEdge->azimuth;

        // ����4�����㲢���·�λ��
        for (const auto& normalizedPair : normalizedObs) {
            const auto& toPoint = normalizedPair.first;
            const auto& normalized = normalizedPair.second;
            double azimuth = startAzimuth + normalized;
            azimuth = std::fmod(azimuth, 360.0);
            if (azimuth < 0) {
                azimuth += 360.0;
            }

            // ��ȡ�򴴽���
            DirectedEdge* edge = manager.getEdge(stationId, toPoint);
            if (!edge) {
                // ��ӣ���Ϊδ�̶�
                edge = manager.addEdge(stationId, toPoint, azimuth, false);
                edgeQueue.push({ stationId, toPoint });

            }
            else {
                // ������δ�̶��ı�
                //if (!edge->isAzimuthFixed()) {
                //    edge->setAzimuth(azimuth);
                //}
                std::cout << "Warning: Edge " << stationId << "��" << toPoint
                    << " already exists " << std::endl;
            }

        }

    }
}

void OrientatingLaunch::buildMatrices() {

    sumEquationRows.clear(); // ���֮ǰ�ļ�¼
	// ����ÿ����������������ÿ����λ�Ƕ�Ӧһ������
    buildEdgeColumnMapping();

    // nΪ�����ĸ���
    const size_t n = edgeColumnMap.size();
    // ����۲�ֵ���д��������NBB��W
    int row = 0;
    // B��NBB��L��W�����������ʽ�����ʹ��eigen�ľ���洢
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
    //        if (obs.type == "����۲�") {
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
    // �������й۲�ֵ
    // �������ѭ������ÿ����վ
    for (const auto& stationData : reader.getObservations()) {
        string stationId = trim(stationData.stationId);

        // obsNum������ǵ�ǰ����Ĺ۲�ֵΪ����վ�ĵڼ����۲�ֵ
        // ͨ����һ���۲�ֵ�ĸ���ֵ���ڶ�Ӧͬ���߷�λ�ǵĸ���ֵ
        // �����Ҫ���⴦��
        int obsNum = 0;
        //initAziKey������¼��ǰ��վ�ĵ�һ������۲�ֵ��Ӧ��edgeColumnMap�е�λ��
        int initEdgeKey = 0;

        // ��¼��ǰ��վ����ʼ�Լ���������۲�ֵ
        DirectedEdge* initialEdge = nullptr;
        DirectedEdge* currentEdge = nullptr;

        // ��ǰ��վ����ʼ����۲�ֵ
		double initObsValue = 0.0;

        // ���嵱ǰ��վ�ͷ��̹۲�ֵ��Ȩֵ
        double p_sum = 0.0;
        // ���嵱ǰ�۲�ֵ�ĺͷ�������
        double L_sum = 0.0;
        // ���嵱ǰ��վ�ĺͷ���ϵ��
        vector<double> B_row_sum(n, 0.0);
        

        // �������ѭ����������ǰ��վ�����й۲�ֵ
        for (const auto& obs : stationData.observations) {
            if (obs.type == "����۲�") {
                obsNum++;
                row++;

                // Ŀ���IDȥ�ո�
                string targetId = trim(obs.pointId);

                // ���ֵ������ɱ߱�ʶ��
                string from, to;
                if (stationId < targetId) {
                    from = stationId;
                    to = targetId;
                }
                else {
                    from = targetId;
                    to = stationId;
                }
                string edgeKey = from + "��" + to;

                // �жϱߵ�ͬ���߷�λ���Ƿ���ڲ�����
                if (edgeColumnMap.find(edgeKey) != edgeColumnMap.end()) {
                    originalObsIndices.push_back(originalObsCounter);
                    isSumEquationRow.push_back(false);

                    double p = 0.0;
                    // ���ù۲�ֵ��Ӧ�ı�Ϊ��֪��ʱĬ�Ϲ۲⾫��Ϊ0.0001�룬�������������
                    // ��ӦȨ��Ҳ�����൱��û�и���
                    if (manager.getEdge(trim(stationId), trim(obs.pointId))->isAzimuthFixed()) {   
                        p = (sigma0 * sigma0) /
                            (0.0001 * 0.0001);
                    }
                    // ���ù۲�ֵ��Ӧ�ı�Ϊδ֪��ʱ���ݶ�ȡ���ļ��������㼴��
                    else {
                        p = (sigma0 * sigma0) /
                            (stod(trim(obs.accuracyNumber)) * stod(trim(obs.accuracyNumber)));
                    }

					

                    if (p == 0) {
                        cerr << "���󹹽�ʧ�ܣ�ȨֵΪ0��" << endl;
                        return;
                    }
                    P_temp.push_back(p);
                    p_sum += 1.0 / p;

                    // ���嵱ǰ�۲�ֵ�ĳ�����
                    double L = 0.0;


                    // ����B����ĵ�k��
                    vector<double> B_row(n, 0.0);


                    // ���۲�ֵΪ�ò�վ�ĵ�һ������۲�ʱ����ʾ�㷽��
                    // ֻ�бߵ�ͬ���߷�λ��ϵ����Ϊ 1.0
                    if (obsNum == 1) {
                        // ������ʼ�۲�ֵ�����ں�������
                        initialEdge = manager.getEdge(trim(stationId), trim(obs.pointId));
                        if (!initialEdge) {
                            cerr << "δ�ҵ���ʼ�ߣ�" << stationId << endl;
                            return;
                        }
                        // ͬʱ������ʼ�ߵľ����к�
                        initEdgeKey = edgeColumnMap[edgeKey];
                        // ͬʱ������ʼ����۲�ֵ
                        initObsValue = stod(trim(obs.value));
                        // ���㵱ǰ���̵ĳ�����㷽������Ϊ0
                        L = 0.0;
                        L_sum += L;

                        // ��� n �Ƿ�����Ч��Χ��
                        if (initEdgeKey >= 0 && initEdgeKey < n) {
                            B_row[initEdgeKey] = 1.0;  // ���ߵ�ͬ���߷�λ��ϵ����Ϊ 1.0
                            B_row_sum[initEdgeKey] += 1.0;
                        }
                        else {
                            std::cout << "Error: build matrix element is out of bounds!" << std::endl;
                        }
                    }
                    // ���۲�ֵΪ�ò�վ����������۲�ʱ���㷽���ϵ����Ϊ -1.0��
                    // �����ߵ�ͬ���߷�λ��ϵ����Ϊ 1.0
                    else {

                        currentEdge = manager.getEdge(trim(stationId), trim(obs.pointId));
                        // ���㷽λ�ǲ�
                        double initialAz = initialEdge->azimuth;
                        double currentAz = currentEdge->azimuth;
                        double delta = fmod(currentAz - initialAz + 360, 360);
						double delta2 = AngleConverter::parseAngleString(trim(obs.value));
                        L = (delta2 - delta);
                        L_sum += L;

                        // ���B����
                        string initialEdgeId = initialEdge->from < initialEdge->to ?
                            (ostringstream() << initialEdge->from << "��" << initialEdge->to).str() :
                            (ostringstream() << initialEdge->to << "��" << initialEdge->from).str();

                        string currentEdgeId = currentEdge->from < currentEdge->to ?
                            (ostringstream() << currentEdge->from << "��" << currentEdge->to).str() :
                            (ostringstream() << currentEdge->to << "��" << currentEdge->from).str();

                        if (edgeColumnMap.find(initialEdgeId) != edgeColumnMap.end()) {
                            B_row[initEdgeKey] = -1.0;  // ����ʼ�ߵ�ͬ���߷�λ��ϵ����Ϊ -1.0
                            B_row_sum[initEdgeKey] += -1.0;
                        }


                        if (edgeColumnMap.find(currentEdgeId) != edgeColumnMap.end()) {
                            B_row[edgeColumnMap[currentEdgeId]] = 1.0; // ������ʼ�ߵ�ͬ���߷�λ��ϵ����Ϊ 1.0
                            B_row_sum[edgeColumnMap[currentEdgeId]] += 1.0;
                        }
                    }

                    // ����B�ĵ�k�ж�NBB�Ĺ���
                    for (int i = 0; i < n; ++i) {
                        for (int j = 0; j < n; ++j) {
                            NBB_temp[i][j] += p * B_row[i] * B_row[j];
                        }
                    }

                    // ����W�ĵ�ǰ�й���
                    double scalar = p * L;
                    for (int j = 0; j < n; ++j) {
                        W_temp[j] += scalar * B_row[j];
                    }

                    B_temp.push_back(B_row); //������ϵ����ӵ�ϵ����B��
					L_temp.push_back(L); //�����г�������ӵ���������L��
                }
                else {
                    cerr << "build Matrix error��Please check the process when initialing edge" << endl;


      //              filteredObsIndices.push_back(originalObsCounter);
      //              if (obsNum == 1) {
      //                  // ������ʼ�۲�ֵ�����ں�������
      //                  initialEdge = manager.getEdge(trim(stationId), trim(obs.pointId));
      //                  if (!initialEdge) {
      //                      cerr << "δ�ҵ���ʼ�ߣ�" << stationId << endl;
      //                      return;
      //                  }

      //                  initEdgeKey = -1;
      //              }
      //              else {
						//// ���ﻹ��һ���߼�δ���ƣ���Ϊ���ܴ���ĳ����վ�ķǵ�һ���۲�ֵ
      //                  // ���ڲ����б�����Σ�
      //              }
                }

                originalObsCounter++;

            }


        }

        if (fabs(p_sum) < epsilon) {
			cerr << "���󹹽�ʧ�ܣ�ȨֵΪ0��" << endl;
			return;
        }
        p_sum = 1.0 / p_sum;

        // ���㵱ǰ��վ�ͷ��̶�NBB�Ĺ���
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                NBB_temp[i][j] += p_sum * B_row_sum[i] * B_row_sum[j];
            }
        }

        // ���㵱ǰ��վ�ͷ��̶�W�Ĺ���
        double scalar = p_sum * L_sum;
        for (int j = 0; j < n; ++j) {
            W_temp[j] += scalar * B_row_sum[j];
        }

        P_temp.push_back(p_sum);
        B_temp.push_back(B_row_sum); //�����кͷ���ϵ����ӵ�����B��
		L_temp.push_back(L_sum); //�����кͷ��̳�������ӵ�����L��
        sumEquationRows.push_back(B_temp.size() - 1); // ��¼�ͷ����е�����
        isSumEquationRow.push_back(true); // ���Ϊ�ͷ�����
    }

	// NBB����ֵ
    NBB = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            NBB(i, j) = NBB_temp[i][j];
        }
    }

    // W��ֵ
	W = Eigen::VectorXd::Zero(n);
	for (int i = 0; i < n; ++i) {
		W(i) = W_temp[i];
	}

    // B����ֵ
    const size_t m = B_temp.size();
    B = Eigen::MatrixXd::Zero(m, n);
    for (int i = 0; i < m; ++i) {
        if (B_temp[i].size() != n) {
            throw std::runtime_error("��ά�Ȳ�һ�£�");
        }
        for (int j = 0; j < n; ++j) {
            B(i, j) = B_temp[i][j];
        }
    }

	 //L��ֵ
	L = Eigen::VectorXd::Zero(m);
	for (int i = 0; i < m; ++i) {
		L(i) = L_temp[i];
	}


	Eigen::MatrixXd PMatrix = Eigen::Map<const Eigen::VectorXd>
        (P_temp.data(), P_temp.size()).asDiagonal();
    P = P_temp;
	Eigen::MatrixXd NBBTest = B.transpose() * PMatrix * B;
	Eigen::MatrixXd WTest = B.transpose() * PMatrix * L;



	// ��ӡ����
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

        // ����Ϊ�̶���ʱ������ƽ��
        //if (it->second->isAzimuthFixed()) continue;

        if (processedEdges.find(reversePair) == processedEdges.end()) {
            ostringstream edgeIdStream;
            if (from < to) {
                edgeIdStream << from << "��" << to;
            }
            else {
                edgeIdStream << to << "��" << from;
            }
            string edgeId = edgeIdStream.str();
            edgeColumnMap[edgeId] = edgeColumnMap.size();
            processedEdges.insert(it->first);
        }
    }
}

void OrientatingLaunch::processInitialObservation(int row,
    const string& stationId, const string& targetId) {
    // ���ұ߻����
    DirectedEdge* edge = manager.getEdge(stationId, targetId);
    if (!edge) {
        edge = manager.getEdge(targetId, stationId);
    }

    // ��ȡ��׼����ID
    std::ostringstream edgeIdStream;
    if (edge->from < edge->to) {
        edgeIdStream << edge->from << "��" << edge->to;
    }
    else {
        edgeIdStream << edge->to << "��" << edge->from;
    }
    string edgeId = edgeIdStream.str();

    // ���B����
    if (edgeColumnMap.find(edgeId) != edgeColumnMap.end()) {
        B(row, edgeColumnMap[edgeId]) = 1.0;
    }
    else {
        cerr << "δ�ҵ�����ӳ�䣺" << edgeId << endl;
    }
}

void OrientatingLaunch::processRegularObservation(int row,
    const string& stationId, const string& targetId, double obsValue) {
    // ������ʼ�ߣ������Ѵ������ʼ�۲⣩
    DirectedEdge* initialEdge = findInitialEdge(stationId);
    if (!initialEdge) {
        cerr << "δ�ҵ���ʼ�ߣ�" << stationId << endl;
        return;
    }

    // ���ҵ�ǰ�۲��
    DirectedEdge* currentEdge = manager.getEdge(stationId, targetId);
    if (!currentEdge) {
        currentEdge = manager.getEdge(targetId, stationId);
    }

    // ���㷽λ�ǲ�
    double initialAz = initialEdge->azimuth;
    double currentAz = currentEdge->azimuth;
    double delta = fmod(currentAz - initialAz + 360, 360);
    L(row) = (obsValue - delta);

    // ���B����
    string initialEdgeId = initialEdge->from < initialEdge->to ?
        (ostringstream() << initialEdge->from << "��" << initialEdge->to).str() :
        (ostringstream() << initialEdge->to << "��" << initialEdge->from).str();

    string currentEdgeId = currentEdge->from < currentEdge->to ?
        (ostringstream() << currentEdge->from << "��" << currentEdge->to).str() :
        (ostringstream() << currentEdge->to << "��" << currentEdge->from).str();

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

    // ���Ҷ�Ӧ�Ĳ�վ�۲�����
    for (const auto& obsData : obsList) {
        if (trim(obsData.stationId) == cleanedStationId) {
            // �����ò�վ�����й۲�ֵ
            for (const auto& obs : obsData.observations) {
                if (obs.type == "����۲�") {
                    double value = AngleConverter::parseAngleString(obs.value);
                    if (fabs(value) < std::numeric_limits<double>::epsilon()) {
                        string targetId = trim(obs.pointId);

                        // ���Ȳ��������
                        DirectedEdge* edge = manager.getEdge(cleanedStationId, targetId);
                        if (!edge) {
                            edge = manager.getEdge(targetId, cleanedStationId);
                        }

                        if (edge) {
                            DebugLog() << "Intialized edge found: " << edge->from
                                << "��" << edge->to
                                << " (azimuth " << std::fixed << std::setprecision(4)
                                << edge->azimuth << "��)";
                        }
                        return edge; // �����ҵ��ĵ�һ����Ч��
                    }
                }
            }
        }
    }

    LogWarning() << "No initialized edge found for station " << cleanedStationId;
    return nullptr;
}




//ƽ�����������
void OrientatingLaunch::performAdjustment() {
    try {
        buildMatrices();        // �ȹ���B��P��L����
        calculateCorrections(); // ����ƽ�����
        updateEdgeAzimuths();   // ���±߷�λ��
        // calculateObservationCorrections(); // ����۲�ֵ����
        calculateAccuracyParameters();     // ���㾫�Ȳ���
        generateResultList();   // ���ɽ���б�
    }
    catch (const std::exception& e) {
        LogCritical() << "ƽ�����ʧ�ܣ�" << e.what();
    }
}

void OrientatingLaunch::calculateCorrections() {
    //int n = B.cols();
    //int m = B.rows();




    // ���㷨���̾��� NBB = B^T * P * B (P��ʱ��Ϊ��λ����)
    //NBB = B.transpose() * B;

    const Eigen::IOFormat fmt(8, 0, ", ", "\n", "", "");


    Eigen::MatrixXd NBB_inv = NBB.inverse();

    // ��������� x = NBB^{-1} * B^T * P * L
    x = NBB_inv * W;
    cout << "x:\n" << x.format(fmt) << endl;
    // ����Э�������� Qx = NBB^{-1}
    Qx = NBB_inv;
    //cout << "Qx:\n" << Qx.format(fmt) << endl;
    // ����۲�ֵ������ V = Bx - L
    V = B * x - L;
    cout << "V:\n" << V.format(fmt) << endl;

}

void OrientatingLaunch::updateEdgeAzimuths() {

    for (auto it = edgeColumnMap.begin(); it != edgeColumnMap.end(); ++it) {
        string edgeId = it->first;
        vector<string> parts = split(edgeId, "��");
        string from = parts[0];
        string to = parts[1];
        int idx = it->second;
        // ʹ�ó�Ա����x
        double correction = x(idx);

        // ��ȡ�����������
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

    // ���㵥λȨ�����
    sigma = std::sqrt(sigma / (m - n));

    // ����Э������� Dx = Qx * sigma^2
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

    // ���߽��
    for (const auto& edgeIdPair : edgeColumnMap) {
        const std::string& edgeId = edgeIdPair.first;
        vector<string> parts = split(edgeId, "��");
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

    // ���۲�ֵ���
    int obsIdx = 0;
    for (const auto& stationData : reader.getObservations()) {
        for (const auto& obs : stationData.observations) {
            if (obs.type == "����۲�" && obsIdx < V.size()) {
                AdjustmentResult res;
                res.FROM = stationData.stationId;
                res.TO = obs.pointId;
                
                // V_results[obsIdx].V = AngleConverter::formatAngleString(V(obsIdx));

                // V_results[obsIdx].RESULT =
                //     AngleConverter::formatAngleString(AngleConverter::parseAngleString(obs.value)
                //                                                           + V(obsIdx));

                // �Զȷ�����ʽչʾ
                res.V = AngleConverter::formatAngleString(VFin(obsIdx));
                res.VALUE = AngleConverter::formatAngleString(
                    AngleConverter::parseAngleString(obs.value));
                res.RESULT = AngleConverter::formatAngleString(
                    AngleConverter::parseAngleString(obs.value)+ VFin(obsIdx));

                // ��ʮ����չʾ   
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
        cout << "  " << key << " �� Column " << value << endl;
    }
}

bool OrientatingLaunch::processKnownAngel() {
    vector<Ui::ObservationData> obsData = reader.getObservations();

    vector<pair<string, Ui::Observation>> fwObservations;

    for (const auto& stationData : obsData) {
        string stationId = trim(stationData.stationId);
        for (const auto& obs : stationData.observations) {
            if (obs.type == "��λ��") {
                fwObservations.push_back({ stationId, obs });
            }
        }

    }

    //��û����֪��λ���������֪�㷴����֪��λ�ǣ�������ǰ����
    if (fwObservations.empty()) {
        cout << "û����֪��λ��!����֪�㷴�㷽λ��";
        return false; // ��֪��λ������Ϊ��ֱ�ӷ���
    }

	// ������֪��λ�����ݣ���edgemanager����ӱ�
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
    // ��ʼ��VFinΪȫ�㣬���ȵ�������ԭʼ����۲���
    VFin = Eigen::VectorXd::Zero(originalObsCounter);

    // �������д�������У������Ч����ֵ
    int vIndex = 0; // V�����е�����
    for (size_t i = 0; i < isSumEquationRow.size(); ++i) {
        if (!isSumEquationRow[i]) {
            // ��ͨ�У���originalObsIndices��ȡԭʼ����
            int originalIdx = originalObsIndices[vIndex];
            VFin[originalIdx] = VNoSum[vIndex];
            vIndex++;
        }
        // �ͷ�����ֱ������
    }

    return VFin;
}


// ��ƽ�����귽λ��ת��Ϊ��ط�λ���Լ����ķ�λ�ǵ�ģ��
void OrientatingLaunch::convertAzimuths() {

    for (auto it = manager.begin(); it != manager.end(); ++it) {
        const auto& edgeKey = it->first;
        const std::string& from = edgeKey.first;
        const std::string& to = edgeKey.second;
        const auto reverseKey = std::make_pair(to, from);

        DirectedEdge* edge = it->second;
        double originalAzimuth = edge->azimuth;

        // ƽ�����귽λ��ת��ط�λ��

        PointCoordAz coordAz;
		coordAz.p1.B = pointsInfo[from].BGeo;
		coordAz.p1.L = pointsInfo[from].LGeo;
		coordAz.p2.B = pointsInfo[to].BGeo;
		coordAz.p2.L = pointsInfo[to].LGeo;
        coordAz.alphaCoord = originalAzimuth;


        coordWGeoAngle(coordAz);
		geoAstroEdgeInfo[edgeKey].geoAzimuth = coordAz.A;
		geoAstroEdgeInfo[edgeKey].geoAzimuthDms = AngleConverter::formatAngleString(coordAz.A);

        //��ط�λ��ת��Ϊ���ķ�λ��
        PointAz pointAz;
        pointAz.pointName = from;
		pointAz.A = coordAz.A; // ��ط�λ��
		pointAz.xi = pointsInfo[from].BAstro - pointsInfo[from].BGeo; // ����ƫ���ϱ�����
        pointAz.eta = (pointsInfo[from].LAstro - pointsInfo[from].LGeo) * cos(pointsInfo[from].BAstro * M_PI / 180.0); // ����ƫ�������
        pointAz.Z = 42.5268; // �춥�� ���ø̼߳��㣬�߳��ļ���ȡ��
        pointAz.L = pointsInfo[from].LGeo; // ��ؾ��� 
        pointAz.lamda = pointsInfo[from].LAstro; // �������ľ��� 
        pointAz.fai = pointsInfo[from].BAstro;  // ��������γ��
		pointAz.sigmaAlpha = 1; // ��λ�Ǿ���

        pointAz.H1 = 10; // ��������ߣ��߳��ļ���ȡ��
        pointAz.Zeta = 0; // �߳��쳣 (����ȡֵΪ0�����Բ��ƣ���
        pointAz.H2 = 21; // �յ������ߣ��߳��ļ���ȡ��
        pointAz.sigmaH1 = 0.1;

        GeodeticCalculator::caclGeoAngleToAstro(pointAz);
		geoAstroEdgeInfo[edgeKey].astroAzimuth = pointAz.alpha_sea;
		geoAstroEdgeInfo[edgeKey].astroAzimuthDms = AngleConverter::formatAngleString(pointAz.alpha_sea);
    }


}


// ������ݴ���ģ��
bool OrientatingLaunch::processTapeData() {
    

    for (partTape& part : tapeReader.TapeData) {
        if (!processPartTapeValue(part)) {
			cerr << "�����ֱ������ʧ�ܣ�" << endl;
			return false;
        }
        if (!calPartTapeValue(part)) {
			cerr << "���㲿�ֱ������ʧ�ܣ�" << endl;
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

	// �Ƚ������۲���׼��ķ�λ�ǳ�ʼ����ʹ��ƽ��֮��ķ�λ�ǣ���manager�����м���
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

	// ������ʼ����ֹ�ߵ�ƽ�����

    part.MzdData[0].distValue = sqrt(pow((pointsInfo[firstMzdId].x - pointsInfo[jzdId].x), 2)
        + pow((pointsInfo[firstMzdId].y - pointsInfo[jzdId].y), 2));
       
	part.MzdData[finMzdIndex].distValue = sqrt(pow((pointsInfo[finalMzdId].x - pointsInfo[jzdId].x), 2)
		+ pow((pointsInfo[finalMzdId].y - pointsInfo[jzdId].y), 2));

    // ������������
	double distJz2MzdFst = part.MzdData[0].distValue;
	double distJz2MzdFin = part.MzdData[finMzdIndex].distValue;


	// �������ʼ������ֹ����ķ���ֵ����Χ��0-360֮��()
	double absDeltaFw = fmod(
		part.MzdData[finMzdIndex].fwValue -
		part.MzdData[0].fwValue + 360.0, 360.0);


    // ������ʼ�ߵ���ֹ�ߵļнǣ�����ĽǶȵķ�Χ��Ȼ��0-180֮�䣬�浽����vertexAngle
    if (absDeltaFw >= 360.0 || absDeltaFw < 0.0) {
		cerr << "����ֵ������Χ��" << endl;
    }
    else if (absDeltaFw > 180.0) {
        part.jzAngleByAzi = 360.0 - absDeltaFw;
	}else {
        part.jzAngleByAzi = absDeltaFw;
	}
    
    // �����Ҷ�������׼�㶥���
    double distMzdFst2MzdFin = part.MzdData[finMzdIndex].tapeValue - part.MzdData[0].tapeValue;
	// �����vertexAngleByCos�����Ƕ���ǵ�����ֵ������Ҫ��һ������
	double jzAngleByCos = (pow(distJz2MzdFst, 2) + pow(distJz2MzdFin, 2) - pow(distMzdFst2MzdFin, 2))
        / (2 * distJz2MzdFst * distJz2MzdFin);
    // �ж�����Ϸ���
    if (jzAngleByCos < -1.0 || jzAngleByCos > 1.0) {
        std::cerr << "Error: Input must be in [-1, 1]." << std::endl;
        return 1;
    }
	jzAngleByCos = acos(jzAngleByCos) * 180.0 / M_PI; // ����ת��Ϊ�Ƕ�
    // ����պϲ���ں�����ƽ��
	double closingError = part.jzAngleByAzi - jzAngleByCos; // �պϲ�
    part.MzdData[finMzdIndex].dirValue = jzAngleByCos + closingError; // ģ��ƽ��Ĺ���

	// ������ʼ�߶�Ӧ����ļнǵ�����ֵ
	double cosMzdFstAngle = (pow(distJz2MzdFst, 2) + pow(distMzdFst2MzdFin, 2) - pow(distJz2MzdFin, 2)) 
        / (2 * distJz2MzdFst * distMzdFst2MzdFin);
    // �жϽ��������
    if (cosMzdFstAngle < -1.0 || cosMzdFstAngle > 1.0) {
        std::cerr << "Error: Input must be in [-1, 1]." << std::endl;
        return 1;
    }
    
    // ����ʣ����׼��ķ���ֵ
    for (int curMzdIndex = 1; curMzdIndex < finMzdIndex; curMzdIndex++) {
        // �̶ȱߵĳ���
        double distMzdFst2MzdCur = part.MzdData[curMzdIndex].tapeValue - part.MzdData[0].tapeValue;

        // �����׼�㵽��ǰ��׼��ľ���
        double distJzd2MzdCur = pow(distJz2MzdFst, 2) + pow(distMzdFst2MzdCur, 2)
            - (2 * distJz2MzdFst * distMzdFst2MzdCur * cosMzdFstAngle);
        distJzd2MzdCur = sqrt(distJzd2MzdCur);
        part.MzdData[curMzdIndex].distValue = distJzd2MzdCur;

        // ���㵱ǰ��׼��ķ���ֵ
		double MzdCurAngle = (pow(distJz2MzdFst, 2) + pow(distJzd2MzdCur, 2) - pow(distMzdFst2MzdCur, 2))
			/ (2 * distJz2MzdFst * distJzd2MzdCur);
		// �ж�����Ϸ���
		if (MzdCurAngle < -1.0 || MzdCurAngle > 1.0) {
			std::cerr << "Error: Input must be in [-1, 1]." << std::endl;
			return 1;
		}
		MzdCurAngle = acos(MzdCurAngle) * 180.0 / M_PI; // ����ת��Ϊ�Ƕ�
		part.MzdData[curMzdIndex].dirValue = MzdCurAngle + closingError / 2.0;
    }

    return true;
}

bool OrientatingLaunch::calPartTapeValue(partTape& part) {



    return true;
}

// ��ȡ��֪��ƽ�������ļ�,����ͬʱ��ƽ�����������ת��������˳�ʼ�Ĵ�ؾ�γ��
bool OrientatingLaunch::readPlainPointsFromFile(const string& plainCoordFilePath) {
    ifstream file(plainCoordFilePath);

    if (!file.is_open()) {
        throw runtime_error("�޷����ļ�: " + plainCoordFilePath);
        return false;
    }

    string line;
    while (getline(file, line)) {
        vector<string> parts;
        string part;
        istringstream iss(line);

        // �����ŷָȥ���ո�
        while (getline(iss, part, ',')) {
            parts.push_back(trim(part));
        }

        if (parts.size() != 3) continue; // ������ʽ������

        try {
            string  id = parts[0];
            double x = stod(parts[1]);
            double y = stod(parts[2]);
			pointsInfo[id].x = x;
			pointsInfo[id].y = y;

			// �����ƽ������ת��Ϊ��ؾ�γ�ȣ�
            // Ĭ��ʹ��bj54�������������Ϊ19������Ϊ6��
            // �����������������뾭�ߡ�NF������Ĭ��ֵ��ʵ����Ҫ��ȡ������ע����ڸ���
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
			// ����ת��ʧ�ܵ���
			continue;
        }
        catch (const exception&) {
            // ����ת��ʧ�ܵ���
            continue;
        }
    }

    return true;
}

// ��ȡ��֪���ľ�γ���ļ���ͬʱ������δ֪���ƽ��������й���õ����е�����ľ�γ��
bool OrientatingLaunch::readAstroBLFromFile(const std::string& AstroBLFilePath) {
    ifstream file(AstroBLFilePath);

    if (!file.is_open()) {
        throw runtime_error("�޷����ļ�: " + AstroBLFilePath);
        return false;
    }

    string line;
    while (getline(file, line)) {
        vector<string> parts;
        string part;
        istringstream iss(line);

        // �����ŷָȥ���ո�
        while (getline(iss, part, ',')) {
            parts.push_back(trim(part));
        }

        if (parts.size() != 3) continue; // ������ʽ������

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
            // ����ת��ʧ�ܵ���
            continue;
        }
        catch (const exception&) {
            // ����ת��ʧ�ܵ���
            continue;
        }
    }



    return true;

}

// ���ĸ���
bool OrientatingLaunch::correct4Centering() {
    bool ifFindRefStation = false;

    // ���������ڹ���ת���Ĳο�վ�Ļ�����Ϣ
    double fai0 = 0.0;
	double fai0_rad = 0.0;
	double lambda0 = 0.0;
	double lambda0_rad = 0.0;
	// �������е㣬�ҵ���һ�����������Ĳο�վ
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
		cerr << "û���ҵ�����ת���Ĳο�վ��������֪���ľ�γ���ļ�����ȷ�ԣ�" << endl;
		return false;
	}

	double e = 0.0; // �����ľ���
    double alpha = 0.0; // �����ط�λ��
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
    //��ط�λ�������ķ�λ�ǵ�ת������������PointAz�ṹ�壬
    // �ýṹ������Ҫ��� alpha�������ķ�λ�ǡ�x������ƫ�����������i��eta������ƫ�����������Z���춥�ࣩ��B�����γ�ȣ���L����ؾ��ȣ���lamda�����ľ��ȣ���fai����������γ�ȣ���
    // ������㺣�����ķ�λ�Ǽ����ط�λ�ǣ���Ҫ���H1��H2�������ߣ���Zeta���߳��쳣�����粻��䣬��������ʹ�õ������ķ�λ�ǽ����ͬ
    //caclAstroAngleByGeo����ط�λ��תΪ���ķ�λ��,PointAz�ṹ��ĵ����ط�λ��A��sigmaA����׼��������ط�λ��A_sea��sigmaA_sea����׼����������ķ�λ��alpha_sea;��������γ��fai_sea
    // caclGeoAngleByAstro�����ķ�λ��תΪ��ط�λ�ǽ����ķ�λ��תΪ��ط�λ�ǣ�PointAz�ṹ������ķ�λ��alpha��sigmaAlpha����׼��������
    //�ṹ����ľ�γ�ȡ��춥�൥λΪ�ȣ�����ƫ���׼��ĵ�λΪ�룬�߶ȡ��߳��쳣��λΪ��
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
    //���귽λ�����ط�λ�ǵ�ת������������PointCoordAz�ṹ�壬
    // �ýṹ���ڰ�������gauss�ṹp1(���)\p2���յ㣩����Ҫ���B\L\L0�����������ߣ����ṹ���ڴ洢��ֵ��λΪ�ȣ���׼�λΪ��
    // ��������ط�λ��A��sigmaA����׼��������������ṹ���ڵ����귽λ��alpha���ط�λ��A�������
    // caclCoordAngleByGeo����ط�λ��תΪ���귽λ�ǣ�caclGeoAngleByCoord�����귽λ��תΪ��ط�λ��
    // Ĭ���������Ϊbj54������Ҫʹ������������������޸Ľṹ��ell����
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

// ������������򳤰��ᡢ����

CoordSystem::PointGauss OrientatingLaunch::gaussBack(double x, double y, int sign, int width)
{
    PointGauss gauss;
    gauss.x = x;
    gauss.y = y;
	gauss.signeded = sign;
    //gauss.L0 = 111;
    gauss.width = width;// 6Ϊ6�ȴ���3Ϊ3�ȴ���0Ϊ�Զ�Ӧ���������ߣ���ָ��L0���м���
    gauss.sigmaB = 1;
    gauss.sigmaL = 1;
    //GeodeticCalculator::GaussForward(gauss, Ellipsoid::bj54);
    //cout << "6�ȴ�" << endl;
    //cout << setprecision(15) << gauss.x << "  " << gauss.y << "  " << gauss.L0 << "  " << degrees2dms(gauss.mca) << " " << degrees2dms(gauss.sigmaMCA) << "  " << gauss.signeded << endl;
    //cout << setprecision(15) << gauss.sigmaX << "  " << gauss.sigmaY << endl;

    cout << endl << "��˹��������" << endl;
    // ʵ�ʵõ���������������ᡢ����
    GeodeticCalculator::GaussBack(gauss, Ellipsoid::bj54);
	cout << setprecision(15) << "  x��" << gauss.x << "  y��" << gauss.y  << "  " <<  " ���ţ� " << gauss.signeded << " ���� " << gauss.width << endl;
    cout << setprecision(15) << "  B��" << degrees2dms(gauss.B) << " L�� " << degrees2dms(gauss.L) << "  " << endl;

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