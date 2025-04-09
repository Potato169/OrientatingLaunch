#include "In2FileReader.h"
#include <iostream>
#include <cstring>
#include <cerrno>

using namespace std;

In2FileReader::In2FileReader() = default;
In2FileReader::~In2FileReader() = default;

// 辅助函数实现
string In2FileReader::trim(const string& s) {
    auto wsfront = find_if_not(s.begin(), s.end(), [](int c) {return isspace(c);});
    auto wsback = find_if_not(s.rbegin(), s.rend(), [](int c) {return isspace(c);}).base();
    return (wsback <= wsfront ? string() : string(wsfront, wsback));
}

bool In2FileReader::readIn2File(const std::string& filePath) {
    ifstream file(filePath, ios::binary);
    if (!file.is_open()) {
        cerr << "无法打开文件！" << endl;
        return false;
    }

    vector<char> data((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
    file.close();

    string text = detectAndConvertEncoding(data);
    text = regex_replace(text, regex("，"), ",");

    istringstream iss(text);
    string line;

    // 读取精度值
    while (getline(iss, line) && line.empty());
    if (!line.empty()) {
        auto parts = splitString(line, ',');
        if (parts.size() >= 3) {
            Ui::AccuracyValue av;
            av.directionPrecision = stod(parts[0]);
            av.distancePrecision1 = stod(parts[1]);
            av.distancePrecision2 = stod(parts[2]);
            accuracyValues.push_back(av);
        }
    }

    // 读取已知点
    while (getline(iss, line)) {
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) break;
        auto parts = splitString(line, ',');
        if (parts.size() >= 3) {
            Ui::KnownPoint point;
            point.id = parts[0];
            point.x = stod(parts[1]);
            point.y = stod(parts[2]);
            knownPoints.push_back(point);
        }
    }

    // 读取观测数据
    while (getline(iss, line)) {
        if (line.empty()) continue;

        Ui::ObservationData obsData;
        obsData.stationId = trim(line);

        while (getline(iss, line) && line.find_first_not_of(" \t\r\n") != std::string::npos) {
            auto parts = splitString(line, ',');
            if (parts.size() >= 3) {
                Ui::Observation obs;
                obs.pointId = parts[0];
                if (parts[1] == "L") obs.type = "方向观测";
                else if (parts[1] == "S") obs.type = "距离观测";
                else if (parts[1] == "A") obs.type = "方位角";
                obs.value = parts[2];
                obsData.observations.push_back(obs);
            }
        }
        observations.push_back(obsData);
    }

    cout << "文件读取成功：" << filePath << endl
        << "已知点数量：" << knownPoints.size() << endl
        << "观测站数量：" << observations.size() << endl;

    return true;
}

std::string In2FileReader::detectAndConvertEncoding(const vector<char>& data) {
    // 检查UTF-8 BOM
    bool hasBOM = (data.size() >= 3 &&
        static_cast<unsigned char>(data[0]) == 0xEF &&
        static_cast<unsigned char>(data[1]) == 0xBB &&
        static_cast<unsigned char>(data[2]) == 0xBF);

    string content(data.begin() + (hasBOM ? 3 : 0), data.end());

    if (hasBOM || isValidUTF8(content)) {
        return content;
    }

    // 尝试GBK转码
    try {
        return convertEncoding("GBK", "UTF-8//IGNORE", content);
    }
    catch (...) {
        cerr << "编码转换失败，使用原始数据" << endl;
        return content;
    }
}

std::string In2FileReader::convertEncoding(const char* from, const char* to, const std::string& input) {
    iconv_t cd = iconv_open(to, from);
    if (cd == (iconv_t)-1) {
        char errorBuffer[256];
        strerror_s(errorBuffer, sizeof(errorBuffer), errno);
        throw runtime_error("iconv_open failed: " + string(errorBuffer));
    }

    size_t inbytes = input.size();
    size_t outbytes = inbytes * 4;
    vector<char> outbuf(outbytes);

    char* inptr = const_cast<char*>(input.data());
    char* outptr = outbuf.data();

    const char* inptr_const = inptr;

    if (iconv(cd, &inptr_const, &inbytes, &outptr, &outbytes) == (size_t)-1) {
        iconv_close(cd);
        char errorBuffer[256];
        strerror_s(errorBuffer, sizeof(errorBuffer), errno);
        throw runtime_error("iconv failed: " + string(errorBuffer));
    }

    iconv_close(cd);
    return string(outbuf.data(), outptr - outbuf.data());
}

bool In2FileReader::isValidUTF8(const std::string& str) const {
    const auto* bytes = reinterpret_cast<const unsigned char*>(str.c_str());
    int num = 0;

    for (int i = 0; bytes[i]; ++i) {
        if (num) { // 后续字节必须为10xxxxxx
            if ((bytes[i] & 0xC0) != 0x80) return false;
            --num;
        }
        else {   // 首字节规则
            if (bytes[i] & 0x80) {
                if ((bytes[i] & 0xE0) == 0xC0) num = 1;
                else if ((bytes[i] & 0xF0) == 0xE0) num = 2;
                else if ((bytes[i] & 0xF8) == 0xF0) num = 3;
                else return false;
            }
        }
    }
    return num == 0;
}

vector<string> In2FileReader::splitString(const string& str, char delimiter) {
    vector<string> tokens;
    string token;
    istringstream tokenStream(str);
    while (getline(tokenStream, token, delimiter)) {
        if (!token.empty()) tokens.push_back(token);
    }
    return tokens;
}