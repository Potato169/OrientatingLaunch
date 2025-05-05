#include "In2FileReader.h"
#include "tapeFileReader.h"
#include <iostream>
#include <iomanip>
#include "EdgeManager.h"
#include "OrientatingLaunch.h"
#include "AngleConverter.h"

int main() {
    In2FileReader reader;
    reader.readIn2File("D:/data/complex.in2");

	tapeFileReader tapeReader;
    //tapeReader.readTapeFile("D:/data/sxbd_2.in2");

	std::cout << std::fixed << std::setprecision(8);

    EdgeManager manager;

    OrientatingLaunch launcher(manager, reader);

    launcher.readPlainPointsFromFile("D:/data/pointCoord.txt");
    // 这里的read函数，同时初始化了这些点的天文经纬度与大地经纬度，
    // 转换为大地方位角时用到的中央经线与椭球参数还有NF都设置为了默认值，
    // 后续根据实际情况更改
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
    for (const auto& res : V_results) {
        std::cout << res.FROM << '\t' << res.TO << '\t' << res.VALUE << '\t' << res.V << '\t' << res.RESULT << std::endl;
        //std::cout << res.FROM << res.TO << AngleConverter::parseAngleString(res.VALUE)
        //    << AngleConverter::parseAngleString(res.V)
        //    << AngleConverter::parseAngleString(res.RESULT);
    }
    
    launcher.convertAzimuths();

    launcher.processTapeData();
         

}