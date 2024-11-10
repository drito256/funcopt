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
            sum += (x[i] - (i+1)) * (x[i] - (i+1));
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

void print_simplex(std::vector<std::vector<double>> &simplex){
    std::cout << "Simplex points:\n";
    for(int i = 0; i < simplex.size(); i++){
        std::cout << "Point " << i << ": [ ";
        for(int j = 0; j < simplex[i].size(); j++){
            std::cout << simplex[i][j];
            if(j != simplex[i].size() - 1)
                std::cout << " , ";
        }
        std::cout <<" ]\n";
    }
}

int main(){
    



    std::cout << "==================================ZAD1=================================" << std::endl;
    constexpr double starting_point = 10; 
    constexpr double epsilon = 10e-6;
    std::pair<double, double> interval = optimize::golden_search(
                                                   functions::parabolic,
                                                   starting_point,
                                                   epsilon);
    std::cout << "Golden ratio search >>> ";
    print_interval(interval);
    std::cout << "----------------------------------------------------------------------\n";



    std::vector<double> starting_point2 = {starting_point};
    std::vector<double> epsilon2 = {epsilon};
    std::vector<double> point = optimize::coord_search(functions::parabolic2,
                                                       starting_point2,
                                                       epsilon2);
    std::cout << "Coordinate search result >>> ";
    print_point(point);
    std::cout << "----------------------------------------------------------------------\n";

    std::vector<std::vector<double>> simplex = optimize::nm_simplex(functions::parabolic2,
                                 starting_point2);

    std::cout << "Nelder-Mead simplex result >>> \n";
    print_simplex(simplex);
    std::cout << "----------------------------------------------------------------------\n";
    
    point = optimize::hooke_jeeves(functions::parabolic2,
                                   starting_point2,
                                   epsilon2);
    std::cout << "Hooke-Jeeves search >>> ";
    print_point(point);

    std::cout << "\n===============================ZAD2=============================\n\n";
    starting_point2 = {-1.9, 2};
    epsilon2 = {epsilon, epsilon};
    point = optimize::coord_search(functions::rosenbrock,
                                   starting_point2,
                                   epsilon2);
    std::cout << "Coord search result >>> ";
    print_point(point);
    std::cout << "----------------------------------------------------------------------\n";
    

    simplex = optimize::nm_simplex(functions::rosenbrock,
                                  starting_point2);
    std::cout << "Nelder-Mead simplex result >>> \n";
    print_simplex(simplex);
    std::cout << "----------------------------------------------------------------------\n";
    

    point = optimize::hooke_jeeves(
                           functions::rosenbrock,
                           starting_point2,
                           epsilon2);
    std::cout << "Hooke-Jeeves search >>> ";
    print_point(point);

    // Second function
    std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
    
    starting_point2 = {0.1, 0.3};
    point = optimize::coord_search(functions::f2,
                                   starting_point2,
                                   epsilon2);
    std::cout << "Coord search result >>> ";
    print_point(point);
    std::cout << "----------------------------------------------------------------------\n";
    

    simplex = optimize::nm_simplex(functions::f2,
                                  starting_point2);
    std::cout << "Nelder-Mead simplex result >>> \n";
    print_simplex(simplex);
    std::cout << "----------------------------------------------------------------------\n";
    

    point = optimize::hooke_jeeves(
                           functions::f2,
                           starting_point2,
                           epsilon2);
    std::cout << "Hooke-Jeeves search >>> ";
    print_point(point);
    std::cout << "----------------------------------------------------------------------\n";

// Third function
    std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
    
    starting_point2 = {0, 0, 0, 0, 0};
    epsilon2 = {epsilon, epsilon, epsilon, epsilon, epsilon};
    point = optimize::coord_search(functions::f3,
                                   starting_point2,
                                   epsilon2);
    std::cout << "Coord search result >>> ";
    print_point(point);
    std::cout << "----------------------------------------------------------------------\n";
    

    simplex = optimize::nm_simplex(functions::f3,
                                  starting_point2);
    std::cout << "Nelder-Mead simplex result >>> \n";
    print_simplex(simplex);
    std::cout << "----------------------------------------------------------------------\n";
    

    point = optimize::hooke_jeeves(
                           functions::f3,
                           starting_point2,
                           epsilon2);
    std::cout << "Hooke-Jeeves search >>> ";
    print_point(point);
    std::cout << "----------------------------------------------------------------------\n";

    // Fourth function
    std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
    
    starting_point2 = {5.1, 1.1};
    epsilon2 = {epsilon, epsilon};
    point = optimize::coord_search(functions::f4,
                                   starting_point2,
                                   epsilon2);
    std::cout << "Coord search result >>> ";
    print_point(point);
    std::cout << "----------------------------------------------------------------------\n";
    

    simplex = optimize::nm_simplex(functions::f4,
                                  starting_point2);
    std::cout << "Nelder-Mead simplex result >>> \n";
    print_simplex(simplex);
    std::cout << "----------------------------------------------------------------------\n";
    

    point = optimize::hooke_jeeves(
                           functions::f4,
                           starting_point2,
                           epsilon2);
    std::cout << "Hooke-Jeeves search >>> ";
    print_point(point);
    std::cout << "----------------------------------------------------------------------\n";

    std::cout << "\n===============================ZAD3=============================\n\n";
    starting_point2 = {5, 5};
    simplex = optimize::nm_simplex(functions::f4,
                                  starting_point2);
    std::cout << "Nelder-Mead simplex result >>> \n";
    print_simplex(simplex);
    std::cout << "----------------------------------------------------------------------\n";
    

    point = optimize::hooke_jeeves(
                           functions::f4,
                           starting_point2,
                           epsilon2);
    std::cout << "Hooke-Jeeves search >>> ";
    print_point(point);
    std::cout << "----------------------------------------------------------------------\n";

    std::cout << "\n===============================ZAD4=============================\n\n";
    starting_point2 = {0.5, 0.5};
    simplex = optimize::nm_simplex(functions::rosenbrock,
                                  starting_point2);
    std::cout << "Nelder-Mead simplex result >>> \n";
    print_simplex(simplex);
    std::cout << "----------------------------------------------------------------------\n";
    

    point = optimize::hooke_jeeves(
                           functions::rosenbrock,
                           starting_point2,
                           epsilon2);
    std::cout << "Hooke-Jeeves search >>> ";
    print_point(point);
    std::cout << "----------------------------------------------------------------------\n";

std::cout << "\n===============================ZAD5=============================\n\n";
    starting_point2 = {0.5, 0.5};
    simplex = optimize::nm_simplex(functions::f6,
                                  starting_point2);
    std::cout << "Nelder-Mead simplex result >>> \n";
    print_simplex(simplex);

                                                        
    return 0;
}

