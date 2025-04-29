#include "FileOperate.h"
namespace FileOperate
{
	vector<string> processLine(string line)
	{
		regex commentRegex("//.*");  // ƥ�� `//` �����ע�Ͳ���
		regex delimiter("[,\\s]+");  // ƥ�䶺�Ż�ո�ָ�
		line = regex_replace(line, std::regex("^\\s+"), "");// ȥ����ǰ�Ŀո�
		line = regex_replace(line, commentRegex, "");// ȥ��ע�Ͳ���
		if (line.empty() || regex_match(line, commentRegex))// ����ǿ��л��� `//` ��ͷ���У������� 
		{
			return vector<string>();
		}
		// ʹ��������ʽ�����Ż�ո�ָ�����
		sregex_token_iterator it(line.begin(), line.end(), delimiter, -1);
		sregex_token_iterator end;
		vector<string> values(it, end);
		return values;
	}
}
