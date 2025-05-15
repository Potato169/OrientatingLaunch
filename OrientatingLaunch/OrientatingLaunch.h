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
	double BAstro = 361.0; // ����Ĭ��ֵΪ361.0����ʾδ��ʼ��
	double LAstro = 361.0;
	std::string BAstroDms;
    std::string LAstroDms;
	CoordSystem::Ellipsoid::Ellipsoid_para ell = CoordSystem::Ellipsoid::bj54;
	double NF = 0.0; 
    double M = 0.0; // ����Ȧ���ʰ뾶
	double N = 0.0; // î��Ȧ���ʰ뾶
	double L0 = 0.0; // ����������

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


    // �ַ�������������
    static std::string trim(const std::string& s);
    static bool startsWith(const std::string& s, const std::string& prefix);
    static std::vector<std::string> split(const std::string& s, const std::string& delimiter);
    std::string doubleToString(const double value, int precision = 8);

	// ת����������
    void convertAzimuths();

    // ��ȡƽ�����꣬�洢��knownPointsCoord�У�����ƽ�����괦��ĵ���ؾ�γ��
    bool readPlainPointsFromFile(const std::string& plainCoordFilePath);

    // ��ȡ���е����ľ�γ�ȣ��������ĵ�����ľ�γ��ͨ�����ļ���õ�
    bool readAstroBLFromFile(const std::string& AstroBLFilePath);

    // ���ĸ���
    bool correct4Centering();

	// ��������ķ�λ��ת���ĺ���
    void geoWstro();

	// ƽ�淽λ�����ط�λ��ת���ĺ���
    void coordWGeoAngle(CoordSystem::PointCoordAz& coordAz);

    // ��˹����(ƽ�����굽��γ��),�����ֱ��Ӧƽ�����ꡢ�����Լ�����
    CoordSystem::PointGauss gaussBack(double x, double y, int sign = 19, int width = 6);

    // ��ӡ��߷�λ��ֵ�������ӡ����ƽ�淽λ�ǣ�
    void printTapeFwValue();

    // geoAstroEdgeInfo�Զ����ϣ����
    struct PairHash {
        template <typename T1, typename T2>
        std::size_t operator()(const std::pair<T1, T2>& p) const {
            // �ֱ��������Ԫ�صĹ�ϣֵ
            std::size_t hash1 = std::hash<T1>{}(p.first);
            std::size_t hash2 = std::hash<T2>{}(p.second);

            // ��Ϲ�ϣֵ�����������µĳ�ͻ��
            return hash1 ^ (hash2 << 1);
        }
    };

private:
    EdgeManager& manager;
    const In2FileReader& reader;

    tapeFileReader& tapeReader;

    std::unordered_set<std::string> processedStations;

	// �洢������꣨���ƣ�������ؾ�γ�����ľ�γ��NF�����뾭�ߡ�î��Ȧ������Ȧ���ʰ뾶���������
    std::unordered_map<std::string, pointCoord> pointsInfo;

	// �洢����ߵĴ�������ķ�λ��(������߷�λ��)
	std::unordered_map<std::pair<std::string, std::string>, edgeAzimuth, PairHash> geoAstroEdgeInfo;

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
    bool convertPartTapeValue(partTape& part);


    // ���������ell�Լ����γ��B����M��N������Ȧ���ʰ뾶�Լ�î��Ȧ���ʰ뾶��
	void calMN(const double& B, const CoordSystem::Ellipsoid::Ellipsoid_para& ell, double& M, double& N);

};

#endif // ORIENTATINGLAUNCH_H