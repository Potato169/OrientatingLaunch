#include "In2FileReader.h"
#include <iostream>
#include <iomanip>
#include "EdgeManager.h"
#include "OrientatingLaunch.h"
#include "AngleConverter.h"

int main() {

    In2FileReader reader;
    reader.readIn2File("D:/data/sxbd_2.in2");


	std::cout << std::fixed << std::setprecision(8);

    EdgeManager manager;

    OrientatingLaunch launcher(manager, reader);
    launcher.initializeEdges();
    manager.printEdges();
	std::cout << "已知边的数量：" << manager.getFixedEdgeCount() << std::endl;
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
}