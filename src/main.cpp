#include "../include/optimize.h"

namespace functions{
    std::function<double(double)> parabolic = [](double x){
        return (x - 3) * (x - 3);
    };

}

void print_interval(std::pair<double, double> &interval){
    std::cout << "Found interval: [" << interval.first << ", " 
                                     << interval.second << "]\n";
}

int main(){


    constexpr double starting_point = 10; 
    constexpr double epsilon = 10e-6;
    std::pair<double, double> interval = optimize::golden_search(
                                                   functions::parabolic,
                                                   starting_point,
                                                   epsilon);
    print_interval(interval);
    return 0;
}

