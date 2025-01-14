#include "../include/optimize.h"

namespace functions{

    // implicit limits
    std::function<bool(std::vector<double>)> ogr1 = [](std::vector<double> x){
        return x[1] - x[0] >= 0;
    };

    std::function<bool(std::vector<double>)> ogr2 = [](std::vector<double> x){
        return 2 - x[0] >= 0;
    };

    std::function<bool(std::vector<double>)> ogr3 = [](std::vector<double> x){
        return 3 - x[0] - x[1] >= 0;
    };
    std::function<bool(std::vector<double>)> ogr4 = [](std::vector<double> x){
        return 3 + 1.5 * x[0] - x[1] >= 0;
    };
    std::function<bool(std::vector<double>)> ogr5 = [](std::vector<double> x){
        return x[1] - 1 == 0;
    };



    std::function<double(double)> parabolic = [](double x){
        return (x - 3) * (x - 3);
    };
    std::function<double(std::vector<double>)> parabolic2 = [](std::vector<double> x){
        return (x[0] - 2) * (x[0] - 2) + (x[1] + 3) * (x[1] + 3);
    };
    std::function<double(std::vector<double>)> parabolic3 = [](std::vector<double> x){
        return (x[0] - 3) * (x[0] - 3) + (x[1]) * (x[1]);
    };

    std::function<double(std::vector<double>)> f2 = [](std::vector<double> x){
        return (x[0] - 4) * (x[0] - 4) + 4 * (x[1] - 2) * (x[1] - 2);
    };

    // Rosenbrock function / Rosenbrock seperated into 2 funcs / partial derivatives
    std::function<double(std::vector<double>)> rosenbrock = [](std::vector<double> x){
        return 100 * ((x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]))  
           + ((1 - x[0]) * ( 1 - x[0]));
    };
    std::function<double(std::vector<double>)> rosenbrock1 = [](std::vector<double> x){
        return 10 * (x[1] - x[0] * x[0]);
    };
    std::function<double(std::vector<double>)> rosenbrock2 = [](std::vector<double> x){
        return 1 - x[0];
    };
    std::function<double(std::vector<double>)> rosenbrock1_x = [](std::vector<double> x){
        return -20 * x[0];
    };
    std::function<double(std::vector<double>)> rosenbrock1_y = [](std::vector<double> x){
        return 10;
    };
    std::function<double(std::vector<double>)> rosenbrock2_x = [](std::vector<double> x){
        return -1;
    };
    std::function<double(std::vector<double>)> rosenbrock2_y = [](std::vector<double> x){
        return 0;
    };


    std::function<double(std::vector<double>)> rosenbrock_x = [](std::vector<double> x){
        return -400 * x[0] * (-x[0] * x[0] + x[1]) + 2 * x[0] - 2;
    };
    std::function<double(std::vector<double>)> rosenbrock_y = [](std::vector<double> x){
        return -200 * x[0] * x[0] + 200 * x[1];
    };
    std::function<double(std::vector<double>)> rosenbrock_xx = [](std::vector<double> x){
        return 1200 * x[0] * x[0] - 400 * x[1] + 2;
    };
    std::function<double(std::vector<double>)> rosenbrock_yy = [](std::vector<double> x){
        return 200;
    };
    std::function<double(std::vector<double>)> rosenbrock_xy = [](std::vector<double> x){
        return -400 * x[0];
    };


    std::function<double(std::vector<double>)> f3 = [](std::vector<double> x){
        double sum = 0;
        for(int i=0;i<x.size();i++){
            sum += (x[i] - (i+1)) * (x[i] - (i+1));
        }
        return sum;
    };

    std::function<double(std::vector<double>)> f4 = [](std::vector<double> x){
        return pow(x[0], 4) / 4 - x[0] * x[0] + 2 * x[0] + pow(x[1] - 1, 2);
    };

    std::function<double(std::vector<double>)> f4_x = [](std::vector<double> x){
        return pow(x[0], 3) - 2 * x[0] + 2;
    };
    std::function<double(std::vector<double>)> f4_xx = [](std::vector<double> x){
        return 3 * x[0] * x[0] - 2;
    };
    std::function<double(std::vector<double>)> f4_y = [](std::vector<double> x){
        return 2 * x[1] - 2;
    };
    std::function<double(std::vector<double>)> f4_yy = [](std::vector<double> x){
        return 2;
    };
    std::function<double(std::vector<double>)> f4_xy = [](std::vector<double> x){
        return 0;
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

    std::function<double(std::vector<double>)> hesse1 = [](std::vector<double> x){
        return 2;
    };
    std::function<double(std::vector<double>)> hesse2 = [](std::vector<double> x){
        return 0;
    };
    std::function<double(std::vector<double>)> hesse3 = [](std::vector<double> x){
        return 0;
    };
    std::function<double(std::vector<double>)> hesse4 = [](std::vector<double> x){
        return 8;
    };
    std::function<double(std::vector<double>)> f2_partial1 = [](std::vector<double> x){
        return 2 * x[0] - 8;
    };
    std::function<double(std::vector<double>)> f2_partial2 = [](std::vector<double> x){
        return 8 * x[1] - 16;
    };

    std::function<double(std::vector<double>)> gn1 = [](std::vector<double> x){
        return x[0] * x[0] + x[1] * x[1] - 1;
    };
    std::function<double(std::vector<double>)> gn2 = [](std::vector<double> x){
        return x[1] - x[0] * x[0];
    };

    std::function<double(std::vector<double>)> jac1 = [](std::vector<double> x){
        return 2 * x[0];
    };
    std::function<double(std::vector<double>)> jac2 = [](std::vector<double> x){
        return 2 * x[1];
    };
    std::function<double(std::vector<double>)> jac3 = [](std::vector<double> x){
        return -2 * x[0];
    };
    std::function<double(std::vector<double>)> jac4 = [](std::vector<double> x){
        return 1;
    };

    std::function<double(std::vector<double>)> gn_large1 = [](std::vector<double> x){
        return x[0] * exp(x[1]) + x[2] - 3;
    };
    std::function<double(std::vector<double>)> gn_large2 = [](std::vector<double> x){
        return x[0] * exp(x[1] * 2) + x[2] - 4;
    };
    std::function<double(std::vector<double>)> gn_large3 = [](std::vector<double> x){
        return x[0] * exp(x[1] * 3) + x[2] - 4;
    };
    std::function<double(std::vector<double>)> gn_large4 = [](std::vector<double> x){
        return x[0] * exp(x[1] * 5) + x[2] - 5;
    };
    std::function<double(std::vector<double>)> gn_large5 = [](std::vector<double> x){
        return x[0] * exp(x[1] * 6) + x[2] - 6;
    };
    std::function<double(std::vector<double>)> gn_large6 = [](std::vector<double> x){
        return x[0] * exp(x[1] * 7) + x[2] - 8;
    };

     std::function<double(std::vector<double>)> jac_large11 = [](std::vector<double> x){
    return exp(x[1]);
    };
    std::function<double(std::vector<double>)> jac_large21 = [](std::vector<double> x){
    return exp(2*x[1]);
    };
    std::function<double(std::vector<double>)> jac_large31 = [](std::vector<double> x){
    return exp(3*x[1]);
    };
    std::function<double(std::vector<double>)> jac_large41 = [](std::vector<double> x){
    return exp(5*x[1]);
    };
    std::function<double(std::vector<double>)> jac_large51 = [](std::vector<double> x){
    return exp(6*x[1]);
    };
    std::function<double(std::vector<double>)> jac_large61 = [](std::vector<double> x){
    return exp(7*x[1]);
    };
    std::function<double(std::vector<double>)> jac_large12 = [](std::vector<double> x){
    return x[0] * exp(x[1]);
    };
    std::function<double(std::vector<double>)> jac_large22 = [](std::vector<double> x){
    return 2 * x[0] * exp(x[1] * 2);
    };
    std::function<double(std::vector<double>)> jac_large32 = [](std::vector<double> x){
    return 3 * x[0] * exp(x[1] * 3);
    };
    std::function<double(std::vector<double>)> jac_large42 = [](std::vector<double> x){
    return 5 * x[0] * exp(x[1] * 5);
    };
    std::function<double(std::vector<double>)> jac_large52 = [](std::vector<double> x){
    return 6 * x[0] * exp(x[1] * 6);
    };
    std::function<double(std::vector<double>)> jac_large62 = [](std::vector<double> x){
    return 7 * x[0] * exp(x[1] * 7);
    };

     std::function<double(std::vector<double>)> jac_last = [](std::vector<double> x){
        return 1;
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
    

    std::cout << "===================ZAD1===================" << std::endl;
    std::vector<std::function<bool(std::vector<double>)>> ogr;
    std::vector<double> stp{-1.9, 2};
    ogr.push_back(functions::ogr1);
    ogr.push_back(functions::ogr2);

    std::vector<std::vector<double>> simpl = optimize::box(functions::rosenbrock,
                                                           ogr,
                                                           std::make_pair(-100,100),
                                                           stp
                                                           );
    print_simplex(simpl);
    stp = {0.1, 0.3};
    std::vector<std::vector<double>> simpl2 = optimize::box(functions::f2,
                                                           ogr,
                                                           std::make_pair(-100,100),
                                                           stp
                                                           );


    print_simplex(simpl2);
    std::cout << "===================ZAD2===================" << std::endl;
    std::vector<std::function<double(std::vector<double>)>> ogr_eq;
    stp = {-1.9, 2};
    std::vector<double> p1 = optimize::penaltyBarrier(functions::rosenbrock,
                                                           ogr,
                                                           ogr_eq,
                                                           1,
                                                           stp,
                                                           1e-6);


    print_point(p1);

    stp = {0.1, 0.3};
    std::vector<double> p2 = optimize::penaltyBarrier(functions::f2,
                                                           ogr,
                                                           ogr_eq,
                                                           1,
                                                           stp,
                                                           1e-6);

    print_point(p2);

    stp = {0, 2}; //tocka iz koje nalazi rj
    p2 = optimize::penaltyBarrier(functions::f2,
                                                           ogr,
                                                           ogr_eq,
                                                           1,
                                                           stp,
                                                           1e-6);

    print_point(p2);
    std::cout << "===================ZAD3===================" << std::endl;
    ogr.clear();
    ogr.push_back(functions::ogr3);
    ogr.push_back(functions::ogr4);

    ogr_eq.clear();
    ogr_eq.push_back(functions::ogr5);
    stp = {0, 0};
    std::vector<double> p3 = optimize::penaltyBarrier(functions::parabolic3,
                                                           ogr,
                                                           ogr_eq,
                                                           1,
                                                           stp,
                                                           1e-6);
    print_point(p3); 

    stp = {5,5};
    std::vector<double> p_interior = optimize::findInteriorPoint(ogr, stp, 1e-6);
    p3 = optimize::penaltyBarrier(functions::parabolic3,
                                                           ogr,
                                                           ogr_eq,
                                                           1,
                                                           p_interior,
                                                           1e-6);
   print_point(p3); 
                                                   
    return 0;
}

