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
                                   const double alpha = 1,
                                   const double beta = 0.5,
                                   const double gamma = 2,
                                   const double sigma = 0.5,
                                   const double e = 10e-6);






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
                                std::vector<std::vector<double>> &simplex);

    static bool nm_exit_condition(std::function<double(std::vector<double>)> func,
                                  std::vector<std::vector<double>> &simplex,
                                  const double epsilon);
    static std::pair<int, int> get_best_worst_index_simplex(
                                    std::function<double(std::vector<double>)> func,
                                    std::vector<std::vector<double>> &simplex);

    static std::vector<double> simplex_reflexion(
                                    std::vector<double> &centroid, 
                                    std::vector<double> &worst_point, 
                                    const double alpha);
    static std::vector<double> simplex_expansion(
                                    std::vector<double> &centroid, 
                                    std::vector<double> &reflex_point, 
                                    const double gamma);
    static std::vector<double> simplex_contraction(
                                    std::vector<double> &centroid, 
                                    std::vector<double> &worst_point, 
                                    const double beta);
 

    static void print_centroid(std::vector<double> &centroid);
    static void print_simplex(std::vector<std::vector<double>> &simplex);


}

#endif



