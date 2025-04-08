#include "iconv.h"
#pragma comment(lib,"libiconv.lib")
#include "In2FileReader.h"
#include <iostream>
#include <iomanip>

int main() {

    In2FileReader reader;
    reader.readIn2File("D:/data/sxbd.in2");
    auto AccuracyValues = reader.getAccuracyValues();
	auto KnownPoints = reader.getKnownPoints();
	auto Observations = reader.getObservations();

	std::cout << std::fixed << std::setprecision(6);
	// ��ӡ����ֵ
	for (const auto& av : AccuracyValues) {
		std::cout << "Direction Precision: " << av.directionPrecision << ", Distance Precision 1: " << av.distancePrecision1 << ", Distance Precision 2: " << av.distancePrecision2 << std::endl;
	}
	// ��ӡ��֪��
	for (const auto& kp : KnownPoints) {
		std::cout << "Known Point ID: " << kp.id << ", X: " << kp.x << ", Y: " << kp.y << std::endl;
	}
	// ��ӡ�۲�����
	for (const auto& obs : Observations) {
		std::cout << "Station ID: " << obs.stationId << std::endl;
		for (const auto& observation : obs.observations) {
			std::cout << "\tObservation Point ID: " << observation.pointId
				<< ", Type: " << observation.type
				<< ", Value: " << observation.value
				<< ", Accuracy Number: " << observation.accuracyNumber
				<< std::endl;
		}
	}
	return 0;
}