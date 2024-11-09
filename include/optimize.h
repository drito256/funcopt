#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>

namespace optimize{
    constexpr double epsilon = 10e-6;

    // TODO: function optimization methods 
    std::pair<double, double> golden_search(std::function<double(double)> func,
                                            const double point,
                                            const double e = epsilon);    
    std::pair<double, double> golden_search(std::function<double(double)> func,
                                            std::pair<double, double> &interval,
                                            const double e = epsilon);

    std::vector<double> coord_search(const std::vector<double> &starting_point, 
                                     const std::vector<double> e);

    
    

}

#endif



