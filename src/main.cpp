#include "../include/optimize.h"

namespace functions{
    std::function<double(double)> parabolic = [](double x){
        return (x - 3) * (x - 3);
    };
    std::function<double(std::vector<double>)> parabolic2 = [](std::vector<double> x){
        return (x[0] - 3) * (x[0] - 3);
    };

    std::function<double(std::vector<double>)> rosenbrock = [](std::vector<double> x){
        return 100 * ((x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]))  
           + ((1 - x[0]) * ( 1 - x[0]));
    };

    std::function<double(std::vector<double>)> f3 = [](std::vector<double> x){
        double sum = 0;
        for(int i=0;i<x.size();i++){
            sum += (x[i] - i) * (x[i] - i);
        }
        return sum;
    };

    std::function<double(std::vector<double>)> f2 = [](std::vector<double> x){
        return (x[0] - 4) * (x[0] - 4) + 4 * (x[1] - 2) * (x[1] - 2);
    };
    
    std::function<double(std::vector<double>)> f4 = [](std::vector<double> x){
        return fabs((x[0] - x[1]) * (x[0] + x[1])) + 
               sqrtf(x[0] * x[0] + x[1] * x[1]); 
    };
    
    std::function<double(std::vector<double>)> f6 = [](std::vector<double> x){
        double sum = 0;
        for(int i=0;i<x.size();i++){
            sum += (x[i] * x[i]);
        }

        return 0.5 + ((sin(sqrt(sum)) * sin(sqrt(sum)))  - 0.5) / 
              ((1 + 0.001 * sum) * (1 + 0.001 * sum));
    };



}




void print_interval(std::pair<double, double> &interval){
    std::cout << "Found interval: [ " << interval.first << " , " 
                                     << interval.second << " ]\n";
}

void print_point(std::vector<double> point){
    std::cout << "Found point : [ ";
    for(int i=0;i<point.size();i++){
        std::cout << point[i];
        if(i != point.size() - 1)
            std::cout << " , ";
    }
    std::cout << " ]\n";
}

int main(){
    



    std::cout << "==================================ZAD1=================================" << std::endl;
    constexpr double starting_point = 2; 
    constexpr double epsilon = 10e-6;
    std::pair<double, double> interval = optimize::golden_search(
                                                   functions::parabolic,
                                                   starting_point,
                                                   epsilon);
    std::cout << "Golden ratio search >>> ";
    print_interval(interval);
    
    std::vector<double> starting_point2 = {starting_point};
    std::vector<double> epsilon2 = {epsilon};
    std::vector<double> point = optimize::coord_search(functions::parabolic2,
                                                       starting_point2,
                                                       epsilon2);
    std::cout << "Coordinate system axis search >>> ";
    print_point(point);






    std::vector<double> starting_point3 = {-1.9, 2};
    std::vector<double> epsilon3 = {epsilon, epsilon};
    std::vector<double> point2 = optimize::coord_search(functions::rosenbrock,
                                                        starting_point3,
                                                        epsilon3);
    print_point(point2);

                                                        
    return 0;
}

