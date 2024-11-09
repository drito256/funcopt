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


std::vector<double> optimize::coord_search(std::function<double(std::vector<double>)> func,
                                           const std::vector<double>& point,
                                           const std::vector<double>& e){

    if(point.size() != e.size()){
       throw "Point dimension does not match epsilon dimension";
    }

    std::vector<double> minimum = point;
    std::vector<double> prev_min;
    do{
        prev_min = minimum;
        for(int i = 0; i < e.size(); i++){
            minimum = axis_min(func, minimum, i, e[i]);
        }
    }while(!compare_points(prev_min, minimum, e));


    return minimum;
}

static bool optimize::compare_points(std::vector<double> p1,
                           std::vector<double> p2,
                           std::vector<double> e){
    
    /*double sum_p1 = 0, avg_p1 = 0;
    double sum_p2 = 0, avg_p2 = 0;
    double sum_e = 0, avg_e = 0;

    // assuming p1 and p2 dimensions are same
    for(int i = 0; i < p1.size(); i++){
        sum_p1 += p1[i];
        sum_p2 += p2[i];
        sum_e += e[i];
    }

    avg_p1 = sum_p1 / p1.size();
    avg_p2 = sum_p2 / p2.size();
    avg_e = sum_e / e.size();

    if(fabs(avg_p1 - avg_p2) < avg_e)
        return true;
    return false;*/

    for(int i = 0;i<e.size();i++){
        if(fabs(p1[i] - p2[i]) > e[i]){
            return false;
        }
    }
    return true;
}

static std::vector<double> optimize::axis_min(
                               std::function<double(std::vector<double>)> func,
                               std::vector<double> &point,
                               int axis_index,
                               double epsilon){

    int direction = 1;
    std::vector<double> test_point = point;
    test_point[axis_index] += direction * epsilon;

    if(func(test_point) > func(point)){
        direction *= -1;
        test_point = point;
        test_point[axis_index] += direction * epsilon;
    }

    std::vector<double> previous_point = point;
    std::vector<double> new_point = previous_point;
    new_point[axis_index] += direction * epsilon;
    
    while(func(new_point) < func(previous_point)){
        previous_point = new_point;
        new_point[axis_index] += direction * epsilon;
    }
    
    return previous_point;
}

