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
		// ����ִ�к���
		func(forward<Args>(args)...);
		// ������ļ�
		ofstream outFile(filename);
		if (!outFile) {
			cerr << "�޷����ļ���" << filename << endl;
			return;
		}
		// ʹ�ø�ʽ����������
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
			cerr << "�޷����ļ���" << readfilename << endl;
			return false;
		}
		if (!outFile) {
			cerr << "�޷����ļ���" << writefilename << endl;
			return false;
		}

		while (!inFile.eof()) {
			if (!formatterIn(inFile, args...)) 
			{
				cerr << "��ȡ����ʧ��" << endl;
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
			cerr << "�޷����ļ���" << readfilename << endl;
			return false;
		}
		if (!outFile) {
			cerr << "�޷����ļ���" << writefilename << endl;
			return false;
		}
		if (!formatterIn(inFile, args...))
		{
			cerr << "��ȡ����ʧ��" << endl;
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
