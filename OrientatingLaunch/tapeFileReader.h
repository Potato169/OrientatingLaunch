#pragma once
#include <string>
#include <vector>
#include <iostream>

// ������׼������ݣ��ݶ���
struct mzdData {
	std::string MzdId; // ��׼��ID��������ʵ�����м�δ����׼�����ƣ�
	double tapeValue; // ��׼��̶�ֵ����
	double dirValue = 0.0; // ����ֵ���ϣ���Ҫ�����������õ�����λĬ��Ϊ����
	double distValue = 0.0; // ����ֵ���ϣ���Ҫ�����������õ�����λĬ��Ϊm
	double fwValue = 0.0; // ��λ��ֵ���ϣ���Ҫ�����������õ�
};

// ����ʵ����׼���Ĺ۲����ݣ��ݶ���
struct partTape {
    std::string JZId; // ��׼��ID
	std::vector<mzdData> MzdData; // �м��۲����ݼ���
};



class tapeFileReader {
public:
	std::vector<partTape> TapeData;

    tapeFileReader();
    ~tapeFileReader();

    bool readTapeFile(const std::string& filePath);

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
    //std::vector<singleTape> tapeData;
    


    std::string detectAndConvertEncoding(const std::vector<char>& data);
    std::string convertEncoding(const char* from, const char* to, const std::string& input);
    bool isValidUTF8(const std::string& str) const;
    static std::vector<std::string> splitString(const std::string& str, char delimiter);
    // �ַ�������������
    static std::string trim(const std::string& s);
};