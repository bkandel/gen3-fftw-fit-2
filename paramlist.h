/**
 * ParamList.h
 * Author: SauravP
 *
 * ParamList is used to read and write parameters
 */
#ifndef _LIBAGYP_PARAMLIST_H_
#define _LIBAGYP_PARAMLIST_H_

#include <string>
#include <iostream>
#include <vector>
#include <complex>          
#include <fstream>
#include <sstream>
#include <unordered_map>

class ParamList {

    public:
        ParamList() : param_file_comments_("#%") {}
        int parse(const std::string& filename);
        int print(const std::string& filename);

        template<typename T>
        int getValue(const std::string s, T& x) {

            std::unordered_map<std::string, std::string>::iterator it;
            it = defmap.find(s);
            if (it == defmap.end())
                return 1;
            std::stringstream ss;
            ss << defmap[s];
            ss >> x;

            return 0;
        }

        template<typename T>
        int getValue(const std::string s, std::vector<T>& x) {

            std::unordered_map<std::string, std::string>::iterator it;
            it = defmap.find(s);
            if (it == defmap.end())
                return 1;
            std::stringstream ss;
            std::string value = defmap[s];

            //
            // erase the [ and ] from the value string
            std::string::size_type pos = value.find_first_of('[');
            if (pos != std::string::npos)  
                value.erase(0, pos+1);
            pos = value.find_last_of(']');
            if (pos != std::string::npos)  
                value.erase(pos);
            ss.str(value);
            T xx; 
            x.clear();
            while (ss >> xx) {
                x.push_back(xx);
            }
            return 0;
        }

        template<typename T> 
        void setValue(const std::string s, T& x) {
            std::stringstream ss;
            ss << x;
            defmap[s] = ss.str(); 
        }

        template<typename T> 
        void setValue(const std::string s, std::vector<T>& x) {
            std::stringstream ss;
            ss << "[";
            if (x.size() > 0) {
                ss << x[0];
                for (unsigned int i=1; i<x.size(); ++i) {
                    ss << " " << x[i];
                }
            }
            ss << "]" << std::endl;
            defmap[s] = ss.str();
        }

    private:
        std::unordered_map<std::string, std::string> defmap;
        const std::string param_file_comments_;
        void trim3(std::string& str);
        void tokenize(const std::string& str,
                std::vector<std::string>& tokens,
                const std::string& delimiters);
};

#endif
