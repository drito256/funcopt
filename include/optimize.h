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

    std::vector<double> coord_search(std::function<double(std::vector<double>)> func,
                                     const std::vector<double> &starting_point, 
                                     const std::vector<double> &e);
    
    std::vector<std::vector<double>> nm_simplex(std::function<double(std::vector<double>)> func,
                                   const std::vector<double> &starting_point,
                                   const double alpha,
                                   const double beta,
                                   const double gamma,
                                   const double e);






    // ========================   STATIC FUNCTIONS ===============================
    // helper functions for coord_search
    static bool compare_points(std::vector<double> p1,
                               std::vector<double> p2,
                               std::vector<double> e);
    static std::vector<double> axis_min(
                                    std::function<double(std::vector<double>)> func,
                                    std::vector<double> &point,
                                    int axis_index,
                                    double epsilon);


    // helper functions for nm_simplex
    static std::vector<double> calculate_centroid(
                                std::vector<std::vector<double>> simplex);
    static void print_centroid(std::vector<double> centroid);


}

#endif



