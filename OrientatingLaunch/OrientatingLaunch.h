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
    ~DebugLog() { std::clog << buffer.str() << std::endl; } // �Զ�����

    template<typename T>
    DebugLog& operator<<(T&& value) {
        buffer << std::forward<T>(value);
        return *this;
    }

private:
    std::ostringstream buffer; // �����������
};

class LogWarning {
public:
    ~LogWarning() { std::cerr << buffer.str() << std::endl; } // �Զ�����

    template<typename T>
    LogWarning& operator<<(T&& value) {
        buffer << std::forward<T>(value);
        return *this;
    }

private:
    std::ostringstream buffer;  // �����������
};

class LogCritical {  // ר�ô�����־��
public:
    ~LogCritical() { std::cerr << buffer.str() << std::endl; } // �Զ����е���׼����

    template<typename T>
    LogCritical& operator<<(T&& value) {
        buffer << std::forward<T>(value);  // ֧��������ʽ����
        return *this;
    }

private:
    std::ostringstream buffer;  // ���ݻ���
};

class OrientatingLaunch {
public:
    //OrientatingLaunch(EdgeManager& manager, const In2FileReader& reader);
	OrientatingLaunch(EdgeManager& manager, const In2FileReader& reader,
        tapeFileReader& tapeReader = tapeFileReader::getDefaultTapeReader());
    OrientatingLaunch(const OrientatingLaunch&) = delete;
    OrientatingLaunch& operator=(const OrientatingLaunch&) = delete;

    //   ����������궨�Ļ���������Ҫ������
    //   1.��ʼ���������е�ÿһ���ߵĽ��Ʒ�λ��
    //   2.����ϵ�����󡢷����̡��۲�ֵ������Ȩ���ƽ������Ҫ������
    //   3.����ƽ������뾫������
	//   4.�������һ����ʽ���
	//   5.�ṩ��Щ����ķ��ʽӿ�
	//   6.���ṩ�˱�߹۲�����ʱ������ִ�б�߹۲����ݵĴ�����ѡ��
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
	const Eigen::VectorXd& getVFin() const { return VFin; }// ����õ���V�������ͷ��̸������Ұ�����û�в���ƽ��Ĺ۲�ֵ������Ϊ���ս��
    const Eigen::VectorXd& getVNoSum() const { return VNoSum; }
    const Eigen::VectorXd& getV() const { return V; }
    const Eigen::VectorXd& getX() const { return x; }
    Eigen::MatrixXd getQx() const { return Qx; }
    double getSigma() const { return sigma; }
    Eigen::MatrixXd getDx() const { return Dx; }
    const std::unordered_map<std::string, int>& getEdgeColumnMap() const { return edgeColumnMap; }
    double getSigma0() const { return sigma0; }
    void printEdgeColumnMap();


	//  ���������ƽ����Ϊƽ�����귽λ��
    // �����ģ�����ڽ�ƽ�����귽λ��ת��Ϊ���ķ�λ�ǻ��ط�λ��
	// �㷨ѡ��ö������
    // �����ڼ����ط�λ�ǵ������㷨���л����ԱȽ��
    enum class GeodeticAlgorithm { Algorithm1, Algorithm2 };

	// ת����������
    void convertAzimuths();

    // ��ȡƽ�����꣬�洢��knownPointsCoord��
    bool readPointsFromFile(const std::string& filePath);


private:
    EdgeManager& manager;
    const In2FileReader& reader;

    tapeFileReader& tapeReader;

    std::unordered_set<std::string> processedStations;
    // �洢���ƽ������
    std::unordered_map<std::string, std::pair<double, double>> pointsCoord;

    Eigen::MatrixXd B;
    Eigen::VectorXd L;
    std::vector<AdjustmentResult> results;
    std::vector<AdjustmentResult> V_results;
    Eigen::VectorXd x;
    Eigen::MatrixXd Qx;
    double sigma = 0.0;
    Eigen::MatrixXd Dx;
	Eigen::VectorXd V;// ����ĸ���ֵ�����ͷ��̸���
    Eigen::VectorXd VNoSum;
	Eigen::VectorXd VFin;// ����ĸ���ֵ�������ͷ��̸������Ұ�����û�в���ƽ��Ĺ۲�ֵ������Ϊ���ս��
    Eigen::MatrixXd NBB;
    Eigen::VectorXd W; // W=BTPL
    std::vector<double> P; // �Խ�Ȩ��P��������ʽ
    // sigma0Ϊ���鵥λȨ�����
	double sigma0 = 0.0;
    std::unordered_map<std::string, int> edgeColumnMap;

    std::vector<int> sumEquationRows; // ���ڼ�¼�ͷ����е�����

	// ��ģ�����ڼ�¼������ƽ��Ĺ۲�ֵ
    std::vector<int> filteredObsIndices; // �洢�����˹۲��ԭʼ����
    std::vector<int> originalObsIndices; // ��¼������Ĺ۲��ԭʼ����
    std::vector<bool> isSumEquationRow; // ����Ƿ�Ϊ�ͷ�����
    int originalObsCounter = 0;     // ȫ��ԭʼ�۲������

    // ����Լ����Ľ������
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

    // �ַ�������������
    static std::string trim(const std::string& s);
    static bool startsWith(const std::string& s, const std::string& prefix);
    static std::vector<std::string> split(const std::string& s, const std::string& delimiter);
    std::string doubleToString(const double value, int precision = 8);

	Eigen::VectorXd calVNoSum();// ����õ���V�޳��˺ͷ��̸���
    Eigen::VectorXd calVFin();// ������VNoSum�Ļ���֮�������û�в���ƽ��Ĺ۲�ֵ����


	// ƽ�����귽λ��ת�����ķ�λ��ģ��ĸ�������
	// �жϸñ��Ƿ��ж�Ӧ������Ϣ
    bool hasAstronomicalInfo(const DirectedEdge* edge) const;
	// ƽ�����귽λ��ת����ط�λ��
    double convertToGeodeticAzimuth(double planeAzimuth) const;
	// ��ط�λ��ת�����ķ�λ��
    double convertToAstronomicalAzimuth(double geodeticAzimuth) const;

	// ������ݴ���������
    // ��ʼ����λ��
    bool processPartTapeValue(partTape& part);
    bool calPartTapeValue(partTape& part);

};

#endif // ORIENTATINGLAUNCH_H