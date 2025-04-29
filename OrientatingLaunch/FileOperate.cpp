#include "FileOperate.h"
namespace FileOperate
{
	vector<string> processLine(string line)
	{
		regex commentRegex("//.*");  // 匹配 `//` 后面的注释部分
		regex delimiter("[,\\s]+");  // 匹配逗号或空格分隔
		line = regex_replace(line, std::regex("^\\s+"), "");// 去除行前的空格
		line = regex_replace(line, commentRegex, "");// 去除注释部分
		if (line.empty() || regex_match(line, commentRegex))// 如果是空行或以 `//` 开头的行，则跳过 
		{
			return vector<string>();
		}
		// 使用正则表达式按逗号或空格分隔数据
		sregex_token_iterator it(line.begin(), line.end(), delimiter, -1);
		sregex_token_iterator end;
		vector<string> values(it, end);
		return values;
	}
}
