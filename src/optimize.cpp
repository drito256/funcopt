#include "../include/optimize.h"

std::pair<double, double> optimize::golden_search(
                                    std::function<double(double)> func,
                                    std::pair<double, double> &interval,
                                    const double e){
    
    const double k = 0.5 * (sqrt(5) - 1);
    double interval_length = interval.second - interval.first;
    double c = interval.second - k * (interval_length);
    double d = interval.first + k * (interval_length);
    double fc = func(c);
    double fd = func(d);

    while((interval.second - interval.first) > e){
        if(fc < fd){
            interval.second = d;
            d = c;
            c = interval.second - k * (interval.second - interval.first);
            fd = fc;
            fc = func(c);
        }
        else{
            interval.first = c;
            c = d;
            d = interval.first + k * (interval.second - interval.first);
            fc = fd;
            fd = func(d);
        }
    }

    return interval;
}

std::pair<double, double> optimize::golden_search(
                                    std::function<double(double)> func,
                                    const double point,
                                    const double e){
    
    std::pair<double, double> interval;

    // case 1
    if(func(point - e) < func(point) && func(point + e) < func(point)){
        throw "Given function is not unimodal";
    }
    // case 2
    else if(func(point - e) > func(point) && func(point + e) > func(point)){
        interval = std::make_pair(point - e, point + e);
    }
    // case 3
    else if(func(point - e) < func(point) && func(point + e) > func(point)){
        double step = e;
        while(func(point - step) < func(point - step / 2)) {
            step *= 2;
        }
        interval = std::make_pair(point - step, point - step / 4);
    }
    // case 4
    else {
        double step = e;
        while(func(point + step) < func(point + step / 2)) {
            step *= 2;
        }
        interval = std::make_pair(point + step / 4, point + step);
    }
    return golden_search(func, interval, e);
}


std::vector<double> optimize::coord_search(const std::vector<double>& point,
                                           const std::vector<double>& e){


}
