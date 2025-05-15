#include "In2FileReader.h"
#include "tapeFileReader.h"
#include <iostream>
#include <iomanip>
#include "EdgeManager.h"
#include "OrientatingLaunch.h"
#include "AngleConverter.h"
#include <windows.h>

// ��ȡexe����Ŀ¼
std::string getExeDirectory() {
    char buffer[MAX_PATH];
    GetModuleFileNameA(NULL, buffer, MAX_PATH);
    std::string exePath(buffer);
    size_t pos = exePath.find_last_of("\\/");
    return exePath.substr(0, pos);
}

// ����Ƿ�Ϊ����·��
bool isAbsolutePath(const std::string& path) {
    if (path.empty()) return false;
#ifdef _WIN32
    return (path.size() >= 3 && isalpha(path[0]) && path[1] == ':' &&
        (path[2] == '\\' || path[2] == '/')) ||
        (path.size() >= 2 && (path[0] == '\\' && path[1] == '\\'));
#else
    return path[0] == '/';
#endif
}

// ����·��
std::string resolvePath(const std::string& exeDir, const std::string& input) {
    if (input.empty()) return "";
    if (isAbsolutePath(input)) {
        return input;
    }
    else {
#ifdef _WIN32
        return exeDir + "\\" + input;
#else
        return exeDir + "/" + input;
#endif
    }
}

// ����ļ��Ƿ����
bool fileExists(const std::string& path) {
    std::ifstream file(path);
    return file.good();
}

int main(int argc, char* argv[]) {
    // ����У��
    // �ĸ������ֱ��ӦIn2�ļ�·�� ����ļ�·�� ���������ļ�·�� ���ľ�γ���ļ�·��
    // In2�ļ�·��Ϊ��Ҫ����������ļ�·��Ϊ��ѡ��������û�в�����Ĭ�ϲ������
    // ������������ģ���ȡ���ݣ�ʵ����Ӧ�ò���������������
    if (argc < 2) {
        std::cerr << "���󣺱���ָ��complex.in2�ļ�·��\n"
            << "�÷���" << argv[0]
            << " in2_file [tape_file] [points_file] [astro_file]\n"
            << "����������Ϊ��ѡ������" << std::endl;
        return 1;
    }

    // ��ȡexeĿ¼������Ĭ��·��
    std::string exeDir = getExeDirectory();
    std::string defaultTape = resolvePath(exeDir, "tapeFile.txt");
    std::string defaultPoints = resolvePath(exeDir, "pointCoord.txt");
    std::string defaultAstro = resolvePath(exeDir, "astroPoints.txt");

    // ��������
    std::string in2Path = resolvePath(exeDir, argv[1]);
    std::string tapePath = (argc >= 3) ? resolvePath(exeDir, argv[2]) : defaultTape;
    std::string pointsPath = (argc >= 4) ? resolvePath(exeDir, argv[3]) : defaultPoints;
    std::string astroPath = (argc >= 5) ? resolvePath(exeDir, argv[4]) : defaultAstro;

    // ����ļ��Ƿ����
    if (!fileExists(in2Path)) {
        std::cerr << "�����޷��ҵ��ļ� " << in2Path << std::endl;
        return 1;
    }


    In2FileReader reader;
    reader.readIn2File(in2Path);

	tapeFileReader tapeReader(tapePath);
    //tapeReader.readTapeFile();

	std::cout << std::fixed << std::setprecision(8);

    EdgeManager manager;

    OrientatingLaunch launcher(manager, reader, tapeReader);

    launcher.readPlainPointsFromFile(pointsPath);
    // �����read������ͬʱ��ʼ������Щ������ľ�γ�����ؾ�γ�ȣ�
    // ת��Ϊ��ط�λ��ʱ�õ������뾭���������������NF������Ϊ��Ĭ��ֵ��
    // ��������ʵ���������

	launcher.readAstroBLFromFile(astroPath);
    launcher.correct4Centering();
    launcher.initializeEdges();
    manager.printEdges();

    launcher.performAdjustment();
    //launcher.printEdgeColumnMap();

    manager.printEdges();

    // // ��ȡ���
    auto results = launcher.getResults();
    // for (const auto& res : results) {
    //     qDebug() << res.FROM << res.TO << res.RESULT << res.M << res.V ;
    // }`

    auto V_results = launcher.getVResults();
	std::cout << "����۲�ֵƽ�������ȷ��룩��" << std::endl;
	std::cout << "FROM\tTO\tFW_VALUE\tV\tRESULT" << std::endl;
    for (const auto& res : V_results) {
        std::cout << res.FROM << '\t' << res.TO << '\t' << res.VALUE << '\t' << res.V << '\t' << res.RESULT << std::endl;
        //std::cout << res.FROM << res.TO << AngleConverter::parseAngleString(res.VALUE)
        //    << AngleConverter::parseAngleString(res.V)
        //    << AngleConverter::parseAngleString(res.RESULT);
    }
    
    launcher.convertAzimuths();

    launcher.processTapeData();
         
	std::cout << "��߹۲�ֵ��������" << std::endl;
	launcher.printTapeFwValue();




}