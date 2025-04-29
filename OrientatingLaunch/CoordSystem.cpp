#include "CoordSystem.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "FileOperate.h"
namespace CoordSystem
{
	const Ellipsoid::Ellipsoid_para Ellipsoid::bj54(6378245, 1 / 298.3,"bj54");
	const Ellipsoid::Ellipsoid_para Ellipsoid::xa80(6378140, 1 / 298.257,"xa80");
	const Ellipsoid::Ellipsoid_para Ellipsoid::wgs84(6378137, 1 / 298.257223563, "wgs84");
	const Ellipsoid::Ellipsoid_para Ellipsoid::cgcs2000(6378137, 1 / 298.257222101, "cgcs2000");
	const Ellipsoid::Ellipsoid_para Ellipsoid::grs80(6378137, 1 / 298.25722210103, "grs80");
	const Ellipsoid::Ellipsoid_para Ellipsoid::pz90(6378136, 1 / 298.25784, "pz90");
}