
#ifndef CORE_UTIL_PARAMETERFILE_H_
#define CORE_UTIL_PARAMETERFILE_H_


#include <sstream>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <list>

#include <algorithm>

#include "CVec.h"

//for linux which does not have __forceinline...
#ifndef WIN32
#define __forceinline inline 
#endif


namespace StringUtil {


	//! converts all chars of a string to lowercase
	__forceinline void toLower(std::string &str) {
		for (size_t i = 0; i < str.length(); i++) {
			if (str[i] <= 'Z' &&  str[i] >= 'A') {
				str[i] -= ('Z' - 'z');
			}
		}
	}

	//! removes all characters from a string
	__forceinline void removeChar(std::string &str, const char c) {
		str.erase(std::remove(str.begin(), str.end(), c), str.end());
	}
	__forceinline std::string removeChar(const std::string &strInput, const char c) {
		std::string str(strInput);
		str.erase(std::remove(str.begin(), str.end(), c), str.end());
		return str;
	}


	//////////////////////
	// native functions //
	//////////////////////
	__forceinline int convertStringToINT(const std::string& s) {
		return atoi(s.c_str());
	}
	__forceinline unsigned int convertStringToUINT(const std::string& s) {
		return (unsigned int)convertStringToINT(s);
	}
	__forceinline double convertStringToDOUBLE(const std::string& s) {
		return atof(s.c_str());
	}
	__forceinline float convertStringToFLOAT(const std::string& s) {
		return (float)convertStringToDOUBLE(s);
	}
	__forceinline char convertStringToCHAR(const std::string& s) {
		return s[0];
	}
/*
	template<class U> inline CVec2<U> convertStringToPoint2D(const std::string& s) {
		point3d<U> ret;
		std::stringstream ss(removeChar(s, 'f'));
		ss >> ret.x >> ret.y;
		return ret;
	}
	template<class U> inline CVec3<U> convertStringToPoint3D(const std::string& s) {
		point3d<U> ret;
		std::stringstream ss(removeChar(s, 'f'));
		ss >> ret.x >> ret.y >> ret.z;
		return ret;
	}
	template<class U> inline CVec4<U> convertStringToPoint4D(const std::string& s) {
		point4d<U> ret;
		std::stringstream ss(removeChar(s, 'f'));
		ss >> ret.x >> ret.y >> ret.z >> ret.w;
		return ret;
	}
*/
	__forceinline bool convertStringToBOOL(const std::string& s) {
		if (s == "false" || s == "0")	return false;
		else return true;
	}

	////////////////////////
	// template overloads //
	////////////////////////

	template<class T>	__forceinline void convertStringTo(const std::string& s, T& res);

	template<>	__forceinline void convertStringTo<int>(const std::string& s, int& res) {
		res = convertStringToINT(s);
	}
	template<>	__forceinline void convertStringTo<unsigned int>(const std::string& s, unsigned int& res) {
		res = convertStringToUINT(s);
	}
	template<>	__forceinline void convertStringTo<double>(const std::string& s, double& res) {
		res = convertStringToDOUBLE(s);
	}
	template<>	__forceinline void convertStringTo<float>(const std::string& s, float& res) {
		res = convertStringToFLOAT(s);
	}
	template<>	__forceinline void convertStringTo<std::string>(const std::string& s, std::string& res) {
		res = s;
	}
	template<>	__forceinline void convertStringTo<char>(const std::string& s, char& res) {
		res = convertStringToCHAR(s);
	}
	template<class U> __forceinline void convertStringTo(const std::string& s, CVec2<U>& res) {
		std::stringstream ss(removeChar(s, 'f'));
		ss >> res.x >> res.y;
	}
	template<class U> __forceinline void convertStringTo(const std::string& s, CVec3<U>& res) {
		std::stringstream ss(removeChar(s, 'f'));
		ss >> res.x >> res.y >> res.z;
	}
	template<class U> __forceinline void convertStringTo(const std::string& s, CVec4<U>& res) {
		std::stringstream ss(removeChar(s, 'f'));
		ss >> res.x >> res.y >> res.z >> res.w;
	}
	template<> __forceinline void convertStringTo<bool>(const std::string& s, bool& res) {
		res = convertStringToBOOL(s);
	}


	//! converts all chars of a string to lowercase (returns the result)
	inline std::string toLowerStr(const std::string& str) {
		std::string res(str);
		for (size_t i = 0; i < res.length(); i++) {
			if (res[i] <= 'Z' &&  res[i] >= 'A') {
				res[i] -= ('Z' - 'z');
			}
		}
		return res;
	}
}

class ParameterFile {
public:
	ParameterFile(char separator = '=', bool caseSensitive = true) {
		m_Separator = separator;
		m_CaseSensitive = caseSensitive;
	}

	ParameterFile(const std::string& filename, char separator = '=', bool caseSensitive = true) {
		m_Separator = separator;
		m_CaseSensitive = caseSensitive;
		addParameterFile(filename);
	}

	void addParameterFile(const std::string& filename) {
		std::ifstream file(filename);
		if (!file.is_open())
		{
			printf("Error to open file %s\n", filename.c_str());
			return;
		}

		while(!file.eof()) {
			std::string line;
			getline(file, line);
			removeComments(line);
			removeSpecialCharacters(line);
			if (line.length() == 0) continue;

			size_t separator = line.find(m_Separator);	//split the string at separator
			if (separator == std::string::npos)	{
				printf("No seperator found in line\n");				
				continue;
			}
			std::string attr_name = line.substr(0, separator);
			std::string attr_value = line.substr(separator + 1, line.length() - 1);
			removeSpecialCharacters(attr_name);
			removeSpecialCharacters(attr_value);
			if (attr_name.length() == 0) {
				printf("Invalid attribute or value");
				continue;
			}
			if (!m_CaseSensitive)	StringUtil::toLower(attr_name);
			m_Values[attr_name] = attr_value;
		}
		file.close();

		bool found = false;
		for (const auto &file : m_Filenames) {
			if (file == filename)	found = true;
		}
		if (!found)	m_Filenames.push_back(filename);
	}

	const void reload() {
		m_Values.clear();
		for(const auto &file : m_Filenames)
		{
			addParameterFile(file);
		}
	}

	template<class T>
	bool readParameter(const std::string& name, T& value) const {
		if (m_CaseSensitive) {
			const auto s = m_Values.find(name);
			if (s == m_Values.end()) {
                std::cout<<("parameter not found: " + name);
				return false; 
			} else {
				StringUtil::convertStringTo(s->second, value);
				return true;
			}
		} else {
			std::string lname(name);	lname = StringUtil::toLowerStr(lname);
			const auto s = m_Values.find(name);
			if (s == m_Values.end()) {
                std::cout<<("parameter not found: " + name);
				return false; 
			} else {
				StringUtil::convertStringTo(s->second, value);
				return true;
			}
		}
	} 
	template<class U>
	bool readParameter(const std::string& name, std::vector<U>& value) const {
		value.clear();
		for (size_t i = 0;; i++) {
			std::stringstream ss;	ss << i;
			std::string currName = name + "[" + ss.str() + "]";
			U currValue;
			if (readParameter(currName, currValue)) {
				value.resize(i+1);
				value[i] = currValue;
			} else {
				break;
			}
		}		
		if (value.size() == 0)	return false;
		else return true;
	}

	template<class U>
	bool readParameter(const std::string& name, std::list<U>& value) const {
		value.clear();
		for (size_t i = 0;; i++) {
			std::stringstream ss;	ss << i;
			std::string currName = name + "[" + ss.str() + "]";
			U currValue;
			if (readParameter(currName, currValue)) {
				value.push_back(currValue);
			} else {
				break;
			}
		}		
		if (value.size() == 0)	return false;
		else return true;
	}
	
	void print() const {
		for (auto iter = m_Values.begin(); iter != m_Values.end(); iter++) {
			std::cout << iter->first << " " << m_Separator << " " << iter->second << std::endl;
		}
	}

private:
	//! removes spaces and tabs at the begin and end of the string
	void removeSpecialCharacters(std::string &str) const {
		char characters[] = {' ', '\t', '\"', ';'};
		const unsigned int length = 4;
		bool found = true;
		while(str.length() && found) {
			found = false;
			for (unsigned int i = 0; i < length; i++) {
				if (*str.begin() == characters[i]) {
					str.erase(str.begin());	found = true;	break;
				}
				if (*(--str.end()) == characters[i]) {
					str.erase(--str.end()); found = true;	break;
				};
			}
		}
	}

	//! searches for comments and removes everything after the comment if there is one
	void removeComments(std::string& s) const {
		std::string comments[] = {"//", "#", ";"};
		const unsigned int length = 3;
		for (unsigned int i = 0; i < length; i++) {
			size_t curr = s.find(comments[i]);
			if (curr != std::string::npos) {
				s = s.substr(0, curr);
			}
		}
	}
	std::map<std::string, std::string> m_Values;
	char m_Separator;
	bool m_CaseSensitive;
	std::list<std::string> m_Filenames;
};



#endif  // CORE_UTIL_PARAMETERFILE_H_
