#include "../include/optimize.h"

namespace functions{
    std::function<double(double)> parabolic = [](double x){
        return (x - 3) * (x - 3);
    };
    std::function<double(std::vector<double>)> parabolic2 = [](std::vector<double> x){
        return (x[0] - 2) * (x[0] + 3);
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

    std::function<double(std::vector<double>)> partial1 = [](std::vector<double> x){
        return 2 * x[0] - 4;
    };
    std::function<double(std::vector<double>)> partial2 = [](std::vector<double> x){
        return 2 * x[1] + 6;
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
    


    std::vector<std::function<double(std::vector<double>)>> partial;
    partial.push_back(functions::partial1);
    partial.push_back(functions::partial2);
    std::cout << "okkk" << std::endl;
    std::vector<double> stp = std::vector<double>{0,0};
    std::vector<double> point = optimize::gradient_desc(
                           functions::parabolic2,
                           partial,
                           stp);
    std::cout << "Gradient descent >>> ";
    print_point(point);
    std::cout << "----------------------------------------------------------------------\n";

                                                   
    return 0;
}

