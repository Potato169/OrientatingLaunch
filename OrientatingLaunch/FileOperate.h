#ifndef FILEOPERATE_H
#define FILEOPERATE_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <regex>
#include <functional>
#include "AngleConversion.h"
#include "CoordSystem.h"
using namespace std;
using namespace AngleConversion;
using namespace CoordSystem;
namespace FileOperate
{
	vector<string> processLine(string line);
	

	template<typename Func, typename Formatter, typename... Args>
	void outputToFileWithFormatter(const string& filename, Formatter formatter, Func func, Args&&... args)
	{
		// 首先执行函数
		func(forward<Args>(args)...);
		// 打开输出文件
		ofstream outFile(filename);
		if (!outFile) {
			cerr << "无法打开文件：" << filename << endl;
			return;
		}
		// 使用格式化器输出结果
		formatter(outFile, forward<Args>(args)...);
		outFile.close();
	}

	template<typename Func, typename Formatter, typename Reader, typename... Args>
	bool inAndOutputToFileWithFormatterOneByOne(const string& readfilename, const string& writefilename, Reader formatterIn, Formatter formatterOut, Func func, Args&... args)
	{
		ifstream inFile(readfilename);
		ofstream outFile(writefilename);
		if (!inFile)
		{
			cerr << "无法打开文件：" << readfilename << endl;
			return false;
		}
		if (!outFile) {
			cerr << "无法打开文件：" << writefilename << endl;
			return false;
		}

		while (!inFile.eof()) {
			if (!formatterIn(inFile, args...)) 
			{
				cerr << "读取数据失败" << endl;
				inFile.close();
				outFile.close();
				return false;
			}
			func(args...);
			formatterOut(outFile, (args)...);
		}
		inFile.close();
		outFile.close();
		return true;
	}

	template<typename Func, typename Formatter, typename Reader, typename... Args>
	bool inAndOutputToFileWithFormatterOnce(const string& readfilename, const string& writefilename, Reader formatterIn, Formatter formatterOut, Func func, Args&... args)
	{
		ifstream inFile(readfilename);
		ofstream outFile(writefilename);
		if (!inFile)
		{
			cerr << "无法打开文件：" << readfilename << endl;
			return false;
		}
		if (!outFile) {
			cerr << "无法打开文件：" << writefilename << endl;
			return false;
		}
		if (!formatterIn(inFile, args...))
		{
			cerr << "读取数据失败" << endl;
			inFile.close();
			outFile.close();
			return false;
		}
		func(args...);
		formatterOut(outFile, (args)...);
		inFile.close();
		outFile.close();
		return true;
	}
}

#endif
