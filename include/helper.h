#ifndef HELPER_H
#define HELPER_H

#include <vector>
#include <string>
#include <cmath>


inline std::vector<double> split_line(const std::string& line) {
    std::vector<double> res;
    int start_index = 0;
    bool in_number = false;

    for (int i = 0; i < line.length(); i++) {
        if (std::isspace(line[i])) {
            if (in_number) {
                res.push_back(std::stod(line.substr(start_index, i - start_index)));
                in_number = false;
            }
        } else {
            if (!in_number) { 
                start_index = i; 
                in_number = true; 
            }
        }
    }

    if (in_number) {
        res.push_back(std::stod(line.substr(start_index)));
    }

    return res;
}

inline bool compare(const double num1, const double num2, const double epsilon){
    if(fabs(num1 - num2) <= epsilon){
        return true;
    }
    return false;
}

#endif
