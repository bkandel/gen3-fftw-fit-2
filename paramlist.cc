/**
 * ParamList.cpp
 * Author: SauravP
 * 
 */
#include "paramlist.h"

using namespace std;

int ParamList::print(const string& filename) {
    ofstream file;
    file.open(filename.c_str());
    if (!file) {
        cerr << "File "
            << filename 
            << " could not be opened."
            << endl;
        return 1;
    }
    for (unordered_map<string, string>::iterator it=defmap.begin(); it!=defmap.end(); ++it) {
        file << it->first << " = " << it->second << endl;
    }
    file.close();
    return 0;
}

int ParamList::parse(const string& filename) {

    ifstream file;
    file.open(filename.c_str());
    if (!file) {
        cerr << "File " 
            << filename 
            << " could not be opened." 
            << endl;
        return 1;
    }
    string line;
    while(getline(file, line)) {
        trim3(line);
        if (!line.empty()) {
            vector<string> tokens;
            tokenize(line, tokens, "=");
            trim3(tokens[0]);
            trim3(tokens[1]);
            defmap[tokens[0]] = tokens[1];
        }
    }
    file.close();

    return 0;
}

void ParamList::trim3(string& str) {
    string::size_type pos = str.find_first_of(param_file_comments_);
    if (pos != string::npos)  
        str.erase(pos); 
    pos = str.find_first_of(";");
    if (pos != string::npos) 
        str.erase(pos);
    pos = str.find_last_not_of(" \t\r");
    if (pos != string::npos) 
        str.erase(pos+1);
    pos = str.find_first_not_of(" \t\r"); 
    if (pos != 0) 
        str.erase(0, pos);
}

void ParamList::tokenize(const string& str,
        vector<string>& tokens,
        const string& delimiters) {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
