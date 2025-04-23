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

// ���弫Сֵ epsilon��
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
        throw runtime_error("��Ҫǡ��2��TWD���2��FWD��");
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

	// ����manager���󣬲����Խ����ıߣ���Щ�߼���Ϊ��ʼ��
    // ������ʼ��edgeQueue���к�������
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

        // ������ʼ���򣨹۲�ֵΪ0��L���ͣ�
        vector<Ui::Observation> validObs;
        for (const auto& obs : data.observations) {
            if (obs.type == "����۲�" && stod(obs.value) == 0.0) {
                validObs.insert(validObs.begin(), obs); // ���ȴ�����ʼ����
            }
            else if (obs.type == "����۲�") {
                validObs.push_back(obs);
            }
        }

        if (validObs.empty()) continue;

        //double baseAzimuth = 0.0;
        //double ObsAzimuth = 0.0;
        //// �������й۲�ֵ
        //for (const auto& obs : validObs) {
        //    if (stod(obs.value) == 0.0) {  
        //        //��ʼ����, �����߼��������⣬���
        //        //��ʼ����δ��ǰ��edge_map�ж��壬��ô����ʼ��ʧ��
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
        //        // �����±�
        //        auto* newEdge = manager.addEdge(stationId, trim(obs.pointId), baseAzimuth);
        //        if (newEdge) edgeQueue.push({ stationId, trim(obs.pointId) });

        //    }
        //    else { // ����ʼ����
        //        //������Ȼʹ�õ���δת��ǰ�ĽǶȣ��������Լ����ĺ���ת��һ�£�ע�ⲻҪtodouble
        //        // double delta = obs.value.toDouble();

        //        double delta = AngleConverter::parseAngleString(obs.value);
        //        ObsAzimuth = fmod(baseAzimuth + delta, 360.0);

        //        // �����±�
        //        auto* newEdge = manager.addEdge(stationId, trim(obs.pointId), ObsAzimuth);
        //        if (newEdge) edgeQueue.push({ stationId, trim(obs.pointId) });

        //    }

        //}

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
    cout << "Qx:\n" << Qx.format(fmt) << endl;
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


void OrientatingLaunch::convertAzimuths(GeodeticAlgorithm algorithm) {
    std::set<std::pair<std::string, std::string>> processedEdges;

    for (auto it = manager.begin(); it != manager.end(); ++it) {
        const auto& edgeKey = it->first;
        const std::string& from = edgeKey.first;
        const std::string& to = edgeKey.second;
        const auto reverseKey = std::make_pair(to, from);

        // �����Ѵ���ķ����
        if (processedEdges.find(reverseKey) != processedEdges.end()) continue;
        processedEdges.insert(edgeKey);

        DirectedEdge* edge = it->second;
        double originalAzimuth = edge->azimuth;

        // ƽ�����귽λ��ת��ط�λ��
        double geodeticAzimuth = 0.0;
        switch (algorithm) {
        case GeodeticAlgorithm::Algorithm1:
            geodeticAzimuth = convertToGeodeticAzimuth(originalAzimuth, algorithm);
            break;
        case GeodeticAlgorithm::Algorithm2:
            geodeticAzimuth = convertToGeodeticAzimuth(originalAzimuth, algorithm);
            break;
        }

        // ���ķ�λ��ת���ж�
        if (hasAstronomicalInfo(edge)) {
            double astronomicalAzimuth = convertToAstronomicalAzimuth(geodeticAzimuth);
            edge->setAzimuth(astronomicalAzimuth);
        }
        else {
            edge->setAzimuth(geodeticAzimuth);
        }
    }
}

// ��������ռλʵ�֣������䣩
bool OrientatingLaunch::hasAstronomicalInfo(const DirectedEdge* edge) const {
    // ���������ķ�λ���ж�����
    return false;
}

double OrientatingLaunch::convertToGeodeticAzimuth(double planeAzimuth, GeodeticAlgorithm algorithm) const {
    // ���������ת���㷨
    return planeAzimuth;
}

double OrientatingLaunch::convertToAstronomicalAzimuth(double geodeticAzimuth) const {
    // ���������ķ�λ��ת���㷨
    return geodeticAzimuth;
}