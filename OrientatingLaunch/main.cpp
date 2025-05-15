#include "In2FileReader.h"
#include "tapeFileReader.h"
#include <iostream>
#include <iomanip>
#include "EdgeManager.h"
#include "OrientatingLaunch.h"
#include "AngleConverter.h"
#include <windows.h>

// 获取exe所在目录
std::string getExeDirectory() {
    char buffer[MAX_PATH];
    GetModuleFileNameA(NULL, buffer, MAX_PATH);
    std::string exePath(buffer);
    size_t pos = exePath.find_last_of("\\/");
    return exePath.substr(0, pos);
}

// 检查是否为绝对路径
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

// 解析路径
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

// 检查文件是否存在
bool fileExists(const std::string& path) {
    std::ifstream file(path);
    return file.good();
}

int main(int argc, char* argv[]) {
    // 参数校验
    // 四个参数分别对应In2文件路径 标尺文件路径 网点坐标文件路径 天文经纬度文件路径
    // In2文件路径为必要参数，标尺文件路径为可选参数，若没有参数则默认不带标尺
    // 后面两个参数模拟获取数据，实际上应该不存在这两个参数
    if (argc < 2) {
        std::cerr << "错误：必须指定complex.in2文件路径\n"
            << "用法：" << argv[0]
            << " in2_file [tape_file] [points_file] [astro_file]\n"
            << "（方括号内为可选参数）" << std::endl;
        return 1;
    }

    // 获取exe目录并设置默认路径
    std::string exeDir = getExeDirectory();
    std::string defaultTape = resolvePath(exeDir, "tapeFile.txt");
    std::string defaultPoints = resolvePath(exeDir, "pointCoord.txt");
    std::string defaultAstro = resolvePath(exeDir, "astroPoints.txt");

    // 解析参数
    std::string in2Path = resolvePath(exeDir, argv[1]);
    std::string tapePath = (argc >= 3) ? resolvePath(exeDir, argv[2]) : defaultTape;
    std::string pointsPath = (argc >= 4) ? resolvePath(exeDir, argv[3]) : defaultPoints;
    std::string astroPath = (argc >= 5) ? resolvePath(exeDir, argv[4]) : defaultAstro;

    // 检查文件是否存在
    if (!fileExists(in2Path)) {
        std::cerr << "错误：无法找到文件 " << in2Path << std::endl;
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
    // 这里的read函数，同时初始化了这些点的天文经纬度与大地经纬度，
    // 转换为大地方位角时用到的中央经线与椭球参数还有NF都设置为了默认值，
    // 后续根据实际情况更改

	launcher.readAstroBLFromFile(astroPath);
    launcher.correct4Centering();
    launcher.initializeEdges();
    manager.printEdges();

    launcher.performAdjustment();
    //launcher.printEdgeColumnMap();

    manager.printEdges();

    // // 获取结果
    auto results = launcher.getResults();
    // for (const auto& res : results) {
    //     qDebug() << res.FROM << res.TO << res.RESULT << res.M << res.V ;
    // }`

    auto V_results = launcher.getVResults();
	std::cout << "方向观测值平差结果（度分秒）：" << std::endl;
	std::cout << "FROM\tTO\tFW_VALUE\tV\tRESULT" << std::endl;
    for (const auto& res : V_results) {
        std::cout << res.FROM << '\t' << res.TO << '\t' << res.VALUE << '\t' << res.V << '\t' << res.RESULT << std::endl;
        //std::cout << res.FROM << res.TO << AngleConverter::parseAngleString(res.VALUE)
        //    << AngleConverter::parseAngleString(res.V)
        //    << AngleConverter::parseAngleString(res.RESULT);
    }
    
    launcher.convertAzimuths();

    launcher.processTapeData();
         
	std::cout << "标尺观测值计算结果：" << std::endl;
	launcher.printTapeFwValue();




}