#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>

#include "matrix.h"

namespace optimize{
    constexpr double epsilon = 10e-6;

    // -------------------------------------------------------------------
    // -------------- METHODS THAT DONT REQUIRE DERIVATIVES --------------
    // -------------------------------------------------------------------
    std::pair<double, double> golden_search (std::function<double(double)> func,
                                            const double point,
                                            const double e = epsilon);    
    std::pair<double, double> golden_search(std::function<double(double)> func,
                                            std::pair<double, double> &interval,
                                            const double e = epsilon);

    std::vector<double> coord_search(std::function<double(std::vector<double>)> func,
                                     const std::vector<double> &starting_point, 
                                     const std::vector<double> &e);
    
    std::vector<std::vector<double>> nm_simplex(
                                   std::function<double(std::vector<double>)> func,
                                   const std::vector<double> &starting_point,
                                   const double alpha = 1,
                                   const double beta = 0.5,
                                   const double gamma = 2,
                                   const double sigma = 0.5,
                                   const double e = epsilon);

    std::vector<double> hooke_jeeves(
                                     std::function<double(std::vector<double>)> func,
                                     const std::vector<double> &starting_point,
                                     std::vector<double> &e);

    
    // -------------------------------------------------------------------
    // ---------------- METHODS THAT REQUIRE DERIVATIVES -----------------
    // -------------------------------------------------------------------
    std::vector<double> gradient_desc(std::function<double(std::vector<double>)> func,
                                      std::vector<std::function<double(std::vector<double>)>> &part_derivative,
                                      const std::vector<double> &starting_point,
                                      const double e = epsilon,
                                      const bool golden_ratio_used = true); 

    std::vector<double> newton_raphson(std::function<double(std::vector<double>)> func,
                                      std::vector<std::function<double(std::vector<double>)>> &part_derivs,
                                      std::vector<std::vector<std::function<double(std::vector<double>)>>> &hesse_matrix,
                                      const std::vector<double> &starting_point,
                                      const double e = epsilon,
                                      const bool golden_ratio_used = true); 
;
    std::vector<double> gauss_newton(std::vector<std::function<double(std::vector<double>)>> &funcs,
                                      std::vector<std::vector<std::function<double(std::vector<double>)>>> &jacobian_matrix,
                                      const std::vector<double> &starting_point,
                                      const double e = epsilon,
                                      const bool golden_ratio_used = true);



    // ========================   PRIVATE FUNCTIONS ===============================
    // helper functions for coord_search
    inline static bool compare_points(std::vector<double> p1,
                               std::vector<double> p2,
                               std::vector<double> e);
    inline static std::vector<double> axis_min(
                                    std::function<double(std::vector<double>)> func,
                                    std::vector<double> &point,
                                    int axis_index,
                                    double epsilon);


    // helper functions for nm_simplex
    inline static std::vector<double> calculate_centroid(
                                std::vector<std::vector<double>> &simplex,
                                int worst_index);

    inline static bool nm_exit_condition(std::function<double(std::vector<double>)> func,
                                  std::vector<std::vector<double>> &simplex,
                                  const double epsilon);
    inline static std::pair<int, int> get_best_worst_index_simplex(
                                    std::function<double(std::vector<double>)> func,
                                    std::vector<std::vector<double>> &simplex);

    inline static std::vector<double> simplex_reflexion(
                                    std::vector<double> &centroid, 
                                    std::vector<double> &worst_point, 
                                    const double alpha);
    inline static std::vector<double> simplex_expansion(
                                    std::vector<double> &centroid, 
                                    std::vector<double> &reflex_point, 
                                    const double gamma);
    inline static std::vector<double> simplex_contraction(
                                    std::vector<double> &centroid, 
                                    std::vector<double> &worst_point, 
                                    const double beta);
    inline static void simplex_shift(std::vector<std::vector<double>> &simplex, 
                              std::vector<double> &best_point, 
                              const double sigma);
 

    inline static void print_centroid(std::vector<double> &centroid);
    inline static void print_simplex(std::vector<std::vector<double>> &simplex);

    inline static std::vector<double> discover(
                                     std::function<double(std::vector<double>)> func,
                                     const std::vector<double> &xp,
                                     const std::vector<double> &delta);

    inline static bool hj_exit_condition(const std::vector<double> &delta,
                                         const std::vector<double> &e);

    inline static bool grad_desc_exit_condition(std::vector<double> &grad,
                                                const double e = epsilon);
    inline static bool new_rap_exit_condition(Matrix delta,
                                              const double e = epsilon);
    inline static double vector_norm(std::vector<double> &vec);
}

#endif



