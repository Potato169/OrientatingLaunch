#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
// ������׼������ݣ��ݶ���
struct mzdData {
	std::string MzdId; // ��׼��ID��������ʵ�����м�δ����׼�����ƣ�
	double tapeValue; // ��׼��̶�ֵ����
	double dirValue = -1.0; // ����ֵ���ϣ���Ҫ�����������õ�����λĬ��Ϊ�Ƕ�
	// ע�⣡����ķ���ֵ��ͬ�ڳ���ķ���ֵ����Ȼ�Ǵ���ʼ������ֹ����ļнǣ�
	// ���Ƿ�Χ��0-180��֮�䣬��ʱ���˳ʱ��ȡ���������ߵĿռ��ϵ
	double distValue = -1.0; // ����ֵ���ϣ���Ҫ�����������õ�����λĬ��Ϊm
	double fwValue = 361.0; // ��λ��ֵ���ϣ���Ҫ�����������õ������귽λ��, ��λ�ȣ�
	std::string fwValueDms; // ��λ��ֵ���ϣ��ȷ����ʽ��
	double reverseFwValue = 361.0; // ����λ��ֵ���ϣ���Ҫ�����������õ������귽λ�ǣ�
	std::string reverseFwValueDms; // ����λ��ֵ���ϣ��ȷ����ʽ��
	double fwValueGeo = 361.0; // ��ط�λ��ֵ���ϣ���Ҫ�����������õ�
	std::string fwValueGeoDms; // ��ط�λ��ֵ���ϣ��ȷ����ʽ��
	double reverseFwValueGeo = 361.0; // �����ط�λ��ֵ���ϣ���Ҫ�����������õ�
	std::string reverseFwValueGeoDms; // �����ط�λ��ֵ���ϣ��ȷ����ʽ��
	double fwValueAstro = 361.0; // ���ķ�λ��ֵ���ϣ���Ҫ�����������õ�
	std::string fwValueAstroDms; // ���ķ�λ��ֵ���ϣ��ȷ����ʽ��
	double reverseFwValueAstro = 361.0; // �������ķ�λ��ֵ���ϣ���Ҫ�����������õ�
	std::string reverseFwValueAstroDms; // �������ķ�λ��ֵ���ϣ��ȷ����ʽ��
};

// ����ʵ����׼���Ĺ۲����ݣ��ݶ���
struct partTape {
    std::string JZId; // ��׼��ID
	std::vector<mzdData> MzdData; // �м��۲����ݼ���
	double jzAngleByAzi = 0.0; // �ɷ�λ�Ǽ������ʼ�ߵ���ֹ�ߵļн� ��λΪ�� ��Χ0-180
};

struct allTape {
	std::string JZId; // ��׼��ID
	std::string FSId; // �����ID
	std::vector<singleTape> allTapeData; // ���б�����ݼ���
};

struct singleTape {
	std::string tapeId; // ���ID
	std::vector<std::string> markId; // ��߿̶�ID����
	std::vector<double> markValue ; // ��߿̶�ֵ����(ƽ���ȷ��)
	std::vector<distObs> distObsData; // ��߹۲�ֵ����
};

struct distObs {
	std::string from; // ��ʼ��ID
	std::string to; // ��ֹ��ID
	double distObaValue; // ����ֵ
};

class tapeFileReader {
public:
	std::vector<partTape> TapeData;

	tapeFileReader();
    tapeFileReader(const std::string& filePath);
    ~tapeFileReader();

	// ��ȡ����ļ�
	bool readTapeFile(const std::string& filePath);

	// ��ȡ����ļ����°汾��
	bool readTapeFileNew(const std::string& filePath);

    static tapeFileReader& getDefaultTapeReader() {
        static tapeFileReader defaultReader;

        partTape testTape;
		testTape.JZId = "JZD0";
		mzdData testMzd1, testMzd2, testMzd3, testMzd4, testMzd5, testMzd6;
		testMzd1.MzdId = "MZD1";
		testMzd1.tapeValue = 0.0;
		testMzd2.MzdId = "MZD1_MZD2_1";
		testMzd2.tapeValue = 1.0;
		testMzd3.MzdId = "MZD1_MZD2_2";
		testMzd3.tapeValue = 2.0;
		testMzd4.MzdId = "MZD1_MZD2_3";
		testMzd4.tapeValue = 3.0;
		testMzd5.MzdId = "MZD1_MZD2_4";
		testMzd5.tapeValue = 4.0;
		testMzd6.MzdId = "MZD2";
		testMzd6.tapeValue = 5.0;



        testTape.MzdData.push_back(testMzd1);
		testTape.MzdData.push_back(testMzd2);
		testTape.MzdData.push_back(testMzd3);
		testTape.MzdData.push_back(testMzd4);
		testTape.MzdData.push_back(testMzd5);
		testTape.MzdData.push_back(testMzd6);

		std::cout << "����δ���������ݣ�ʹ�ò��Ա�����ݽ��в���" << std::endl;

		defaultReader.TapeData.push_back(testTape);

        return defaultReader;
    };

private:
    static std::vector<std::string> splitString(const std::string& str, char delimiter);
    // �ַ�������������
    static std::string trim(const std::string& s);

	// ����Ƿ�Ϊ�����
	bool isTapeLine(const std::string& s) {
		return s.find("���") == 0;
	}
};