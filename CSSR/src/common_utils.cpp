#include "common_utils.h"
#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <sys/stat.h>
#include <cerrno>
#include <climits>

using namespace std;

//Determine if the path <directory>/<file> exists and (is not itself a directory)
bool fileExists(const char* const directory, const char* const file) {
	if (is_blank(directory)||is_blank(file)) return false;
	string filePath(directory);
	return fileExists(filePath.append("/").append(file).c_str());
}
bool fileExists(const char* const fullPath, const bool verifyReadable) {
	if (is_blank(fullPath)) return false;
	if (verifyReadable) {
		// Verify that the file can be opened for reading.
		ifstream in_test(fullPath);
		return in_test.good();
	} else {
		// Verify only that the file exists (and is not a directory)
		struct stat info;
		return stat(fullPath, &info)==0 && (info.st_mode & S_IFDIR)==0;
	}
}
bool dirExists(const char* const fullDirPath) {
	struct stat info;
	return !is_blank(fullDirPath) && stat(fullDirPath, &info)==0 && (info.st_mode&S_IFDIR)==S_IFDIR;
}
// Determine whether the specified file path signifies that STDIO should
// be used instead of a regular file. Currently this means the file path must be
// identically "-" (as commonly used in POSIX programs for the same purpose).
bool isStdIoFile(const char* const fullPath) {
	return fullPath!=NULL && fullPath[0]=='-' && fullPath[1]=='\0';
}

// Gets the portion of a path that represents the filename.
// I.e. removes the directory portion of a path. 
// If removeExtension is true, the file extension is also removed.
string getFileName(const char * const path, bool removeExtension) {
	string filename(path);
	std::size_t pos = filename.find_last_of("/\\");
	if (pos != string::npos) filename.erase(0, pos + 1);
	if (removeExtension) {
		pos = filename.rfind('.');
		if (pos != string::npos)
			filename.erase(pos);
	}
	return filename;
}

// Creates string and fills it by sprintf-formatting the arguments.
string sfmt(const char* const format, ...) {
    va_list args;
	int size = strlen(format) + 256;
	char* buf = new char[size];
    va_start(args, format);
    int req_size = vsnprintf(buf, size, format, args);
    va_end(args);
	if (req_size < 0)
		// negative return values indicate an error. 
		req_size = sprintf(buf, "Error formatting arguments: %d", req_size);
	else if (req_size >= size) {  // indicates required buffer size (not including terminating \0)
		// the formatted arguments could not all fit inside the buffer. So resize the buffer to the exact required amount and repeat.
		delete[] buf;
		size = req_size+1; // +1 for \0
		buf = new char[size];
		va_start(args, format);
		vsnprintf(buf, size, format, args);
		va_end(args);
	}
	string ret = buf; // copy the buffer and release it.
	delete[] buf;
	return ret;
}

// Find a char in a c-string
std::size_t findchr(const char* const subject, const char find) {
	const char* found = strchr(subject, find); // strchr returns a pointer to the char that was found or NULL (if not found)
	return found==NULL ? string::npos : found - subject;
}

string& trimLeft(string &s) {
	string::iterator it;
	for (it = s.begin(); it!=s.end(); it++) if (!::isspace(*it)) break; // we are guaranteed to find a non-whitespace character (the final '\0' at s.end() if nothing else)
	s.erase(s.begin(), it);
	return s;
}
string& trimRight(string &s) {
	string::iterator it;
	for (it = s.end()-1; it>=s.begin(); it--) if (!::isspace(*it)) { it++; break; } // start with s.end()-1 to skip the final '\0' at s.end()
	if (it < s.begin()) it++; // increment so we don't delete the NON-whitespace character that was found.
	s.erase(it, s.end());
	return s;
}
string& trim(string &s) {
	trimLeft(s);
	if (!s.empty())
		trimRight(s);
	return s;
}
// Converts a string to lower-case (operates on subject in-place).
string& toLower(string &s) {
	if (!s.empty())
		transform(s.begin(), s.end(), s.begin(), ::tolower);
	return s;
}
// Converts a string to upper-case (operates on subject in-place).
string& toUpper(string &s) {
	if (!s.empty())
		transform(s.begin(), s.end(), s.begin(), ::toupper);
	return s;
}

// The following string functions operate on const strings and 
// return a copy of the string with the desired modifiations.
string trimLeft(const string &subject) { string copy(subject); return trimLeft(copy);  }
string trimRight(const string &subject) { string copy(subject); return trimRight(copy);  }
string trim(const string &subject) { string copy(subject); return trim(copy);  }
string toLower(const string &subject) { string copy(subject); return toLower(copy);  }
string toUpper(const string &subject) { string copy(subject); return toUpper(copy);  }

// Parse a c-string to get an integer.
// Returns true on success or false if the string could not be converted to an interger.
// Note that intResult is NOT modified if conversion fails.
bool parseInt(const char* const input, int& intResult, const bool readToEnd) {
	long lnum; char *endptr;
	errno = 0; // reset global error indicator
	// long int strtol (const char* input, char** endptrptr, int base);
	lnum = strtol(input, &endptr, 0); // 0 specifies parsing of standard radix options (including 0x0 077, etc)
	if ( endptr==input //if no characters were converted, these pointers are equal.
		  // explicitly check for overflows. note that long *might* be the same size as int on some systems.
		  || errno == ERANGE || (lnum > INT_MAX) || (lnum < INT_MIN)  // number out of range
		  || (readToEnd&&(*nextNonWhitespaceChar(endptr)!='\0')) // extra (non-whitepsace) text after number.
		) return false;
	 //Success. Convert the result to a plain int.
	intResult = (int)lnum;
	return true;
}

// Parse a c-string to get a double.
// Return true on success or false if the string could not be converted to a double.
// Note that dblResult is NOT modified if conversion fails.
bool parseDbl(const char* input, double& dblResult, const bool readToEnd) {
	double num; char *endptr;
	errno = 0; // reset global error indicator
	//double strtod(const char *input, char **endptr);
	num = strtod(input, &endptr);
	if ( endptr==input //if no characters were converted these pointers would be equal
		 || errno!=0 // check global error indicator for non-zero values (e.g. ERANGE)
		 || (readToEnd&&(*nextNonWhitespaceChar(endptr)!='\0'))
		 ) return false;
	//Success.
	dblResult = num;
	return true;
}

// split a string using a delimiter.
vector<string> split(const string& str, const string& delim, const bool includeEmptyValues) {
    vector<string> tokens;
    std::size_t prev = 0, pos = 0;
    while (pos < str.length()) {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        //cout<<"\t"<<includeEmptyValues<<"\t"<<(pos>prev)<<endl;
        if (includeEmptyValues||pos>prev) tokens.push_back(token);
        prev = pos + delim.length();
    }
    return tokens;
}

NullStream NullStream::Default; // construct default member.
