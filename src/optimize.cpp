#include "../include/optimize.h"


int eval_counter = 0;
int grad_counter = 0;
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

    eval_counter = 2;

    while((interval.second - interval.first) > e){
        if(fc < fd){
            interval.second = d;
            d = c;
            c = interval.second - k * (interval.second - interval.first);
            fd = fc;
            fc = func(c);
            eval_counter++;
        }
        else{
            interval.first = c;
            c = d;
            d = interval.first + k * (interval.second - interval.first);
            fc = fd;
            fd = func(d);
            eval_counter++;
        }
    }
    // std::cout << "Number of evaluations: " << eval_counter << std::endl;
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
    eval_counter = 0;
    std::vector<double> minimum = point;
    std::vector<double> prev_min;
    do{
        prev_min = minimum;
        for(int i = 0; i < e.size(); i++){
            minimum = axis_min(func, minimum, i, e[i]);
        }
    }while(!compare_points(prev_min, minimum, e));


    // std::cout << "Number of evaluations: " << eval_counter << std::endl;
    return minimum;
}

inline static bool optimize::compare_points(std::vector<double> p1,
                           std::vector<double> p2,
                           std::vector<double> e){
    for(int i = 0;i<e.size();i++){
        if(fabs(p1[i] - p2[i]) > e[i]){
            return false;
        }
    }
    return true;
}

inline static std::vector<double> optimize::axis_min(
                               std::function<double(std::vector<double>)> func,
                               std::vector<double> &point,
                               int axis_index,
                               double epsilon){

    int direction = 1;
    std::vector<double> new_point = point;
    new_point[axis_index] += direction * epsilon;
    
    eval_counter +=2;
    if(func(new_point) > func(point)){
        direction *= -1;
        new_point = point;
        new_point[axis_index] += direction * epsilon;
    }

    std::vector<double> previous_point = point;
    new_point = previous_point;
    new_point[axis_index] += direction * epsilon;
    
    while(func(new_point) < func(previous_point)){
        previous_point = new_point;
        new_point[axis_index] += direction * epsilon;
        eval_counter +=2;
    }
    
    return previous_point;
}

std::vector<std::vector<double>> optimize::nm_simplex(
                              std::function<double(std::vector<double>)> func,
                              const std::vector<double> &starting_point,
                              const double alpha,
                              const double beta,
                              const double gamma,
                              const double sigma,
                              const double e){

    // generate starting simplex
    const double offset = 1;
    std::vector<std::vector<double>> simplex;
    simplex.push_back(starting_point);
    for(int i = 0; i < starting_point.size();i++){
        std::vector<double> new_point = starting_point;
        new_point[i] += offset;
        simplex.push_back(new_point);
    }
    eval_counter = 0;
    do {
        auto best_worst_index = get_best_worst_index_simplex(func, simplex);
        auto centroid = calculate_centroid(simplex, best_worst_index.second);
        auto r_point = simplex_reflexion(centroid, 
                                         simplex[best_worst_index.second], alpha);
       
        eval_counter += 3; //average 2 and 4
        if (func(r_point) < func(simplex[best_worst_index.first])) {
            auto e_point = simplex_expansion(centroid, r_point, gamma);

            eval_counter +=2;
            if (func(e_point) < func(simplex[best_worst_index.first])) {
                simplex[best_worst_index.second] = e_point;
            } 
            else {
                simplex[best_worst_index.second] = r_point;
            }
        } 
        else if (func(r_point) < func(simplex[best_worst_index.second])) {
            simplex[best_worst_index.second] = r_point;
        } 
        else {
            auto c_point = simplex_contraction(centroid, 
                                               simplex[best_worst_index.second], beta);
            
            eval_counter += 2;
            if (func(c_point) < func(simplex[best_worst_index.second])) {
                simplex[best_worst_index.second] = c_point;
            }
            else {
                simplex_shift(simplex, simplex[best_worst_index.first], sigma);
            }
        }
    } while (!nm_exit_condition(func, simplex, e));

    // std::cout << "Number of evaluations: " << eval_counter << std::endl;
    return simplex;
}

inline static std::vector<double> optimize::calculate_centroid(
                                        std::vector<std::vector<double>> &simplex,
                                        int worst_index){

    std::vector<double> centroid(simplex[0].size());

    for(int i = 0; i < simplex.size(); i++){
        if(i == worst_index)
            continue;
        for(int j = 0; j < simplex[i].size(); j++){
            centroid[j] += simplex[i][j];
        }
    }

    for(int i = 0; i < centroid.size(); i++){
        centroid[i] /= (simplex.size() - 1);
    }
    return centroid;
}

inline static bool optimize::nm_exit_condition(
                              std::function<double(std::vector<double>)> func,
                              std::vector<std::vector<double>> &simplex,
                              const double epsilon){
    
    int worst_index = get_best_worst_index_simplex(func, simplex).second;
    std::vector<double> centroid = calculate_centroid(simplex, worst_index);
    double sum = 0;
    for(int i = 0; i < simplex.size(); i++){
        sum += ((func(simplex[i]) - func(centroid)) * 
               (func(simplex[i]) - func(centroid)));
    }

    double root_avg = sqrtf(sum / simplex.size());
    
    return root_avg <= epsilon;
}

inline static std::pair<int, int> optimize::get_best_worst_index_simplex(
                                    std::function<double(std::vector<double>)> func,
                                    std::vector<std::vector<double>> &simplex){
    int best = 0, worst = 0;

    for(int i = 0; i < simplex.size(); i++){
        eval_counter += 3;
        double value = func(simplex[i]);
        if(value > func(simplex[worst])){
            worst = i;
        }
        if(value < func(simplex[best])){
            best = i;
        }
    }

    return std::make_pair(best, worst);
}

inline static std::vector<double> optimize::simplex_reflexion(
                                        std::vector<double> &centroid,
                                        std::vector<double> &worst_point,
                                        const double alpha){
    std::vector<double> r_point(centroid.size());
    for(int i = 0; i < centroid.size(); i++){
        r_point[i] = (1 + alpha) * centroid[i] - alpha * worst_point[i]; 
    }

    return r_point;
}

inline static std::vector<double> optimize::simplex_expansion(
                                                std::vector<double> &centroid,
                                                std::vector<double> &reflex_point,
                                                const double gamma){
    std::vector<double> expan_point(centroid.size());
    for(int i = 0; i < centroid.size(); i++){
        expan_point[i] = (1 - gamma) * centroid[i] + gamma * reflex_point[i]; 
    }

    return expan_point;
}

inline static std::vector<double> optimize::simplex_contraction(
                                                std::vector<double> &centroid,
                                                std::vector<double> &worst_point,
                                                const double beta){
    std::vector<double> contr_point(centroid.size());
    for(int i = 0; i < centroid.size(); i++){
        contr_point[i] = (1 - beta) * centroid[i] + beta * worst_point[i]; 
    }

    return contr_point;
}

inline static void optimize::simplex_shift(std::vector<std::vector<double>> &simplex,
                                    std::vector<double> &best_point,
                                    const double sigma){

    for(int i = 0; i < simplex.size(); i++){
        for(int j = 0; j < simplex[i].size(); j++){
            simplex[i][j] = 0.5 * (simplex[i][j] + best_point[j]);
        }
    }
}

inline static void optimize::print_simplex(std::vector<std::vector<double>> &simplex){
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

inline static void optimize::print_centroid(std::vector<double> &centroid){
    std::cout << "Centroid : [ ";
    for(int i = 0; i < centroid.size(); i++){
        std::cout << centroid[i];
        if(i != centroid.size() - 1)
            std::cout << centroid[i] << " , ";
    }
    std::cout << "]\n";
}


std::vector<double> optimize::hooke_jeeves(
                                    std::function<double(std::vector<double>)> func,
                                    const std::vector<double> &starting_point,
                                    std::vector<double> &e){
    std::vector<double> xb;
    std::vector<double> xp;
    xb = xp = starting_point;
    std::vector<double> xn;

    std::vector<double> delta(starting_point.size());
    for(int i = 0;i<delta.size();i++){
        delta[i] = 0.5;
    }
    eval_counter = 0;
    do{
        xn = discover(func, xp, delta); 
        eval_counter += 2;
        if(func(xn) < func(xb)){
            for(int i = 0; i < xp.size(); i++){
                xp[i] = 2 * xn[i] - xb[i];
            }
            xb = xn;
        }
        else{
            for(int i = 0; i < delta.size(); i++){
                delta[i] = delta[i] / 2;
            }
            xp = xb;
        }
    }while(!hj_exit_condition(delta, e));

    // std::cout << "Number of evaluations: " << eval_counter << std::endl;
    return xb;
}

static inline std::vector<double> optimize::discover(
                                    std::function<double(std::vector<double>)> func,
                                    const std::vector<double> &xp, 
                                    const std::vector<double> &delta){

    std::vector<double> x = xp;
    for(int i = 0; i < x.size(); i++){
        eval_counter +=2;
        double p = func(x);
        x[i] = x[i] + delta[i];
        double n  = func(x);

        if(n > p){
            x[i] = x[i] - 2 * delta[i];
            eval_counter++;
            n = func(x);
            if(n > p){
                x[i] = x[i] + delta[i];
            }
        }
    }
    return x;
}

static inline bool optimize::hj_exit_condition(const std::vector<double> &delta,
                                               const std::vector<double> &e){
    for(int i = 0; i < e.size(); i++){
        if(delta[i] >= e[i]){
            return false;
        }
    }
    return true;
}


std::vector<double> optimize::gradient_desc(std::function<double(std::vector<double>)> func,
                                  std::vector<std::function<double(std::vector<double>)>> &part_deriv,
                                  const std::vector<double> &starting_point,
                                  const double e,
                                  const bool golden_ratio_used){

    std::vector<double> x = starting_point;
    std::vector<double> grad(x.size());
    std::vector<double> new_x(x.size());
    std::vector<double> test_convergence_point(x.size());
    test_convergence_point = starting_point;
    int i = 1;
    eval_counter = 0;
    do {
        for (int i = 0; i < x.size(); i++) {
            grad[i] = part_deriv[i](x);
        }
        grad_counter++;

        if (golden_ratio_used) {
            auto line_search_func = [&](double lambda) {
                std::vector<double> temp_x = x;
                for (int i = 0; i < x.size(); i++) {
                    temp_x[i] -= lambda * grad[i];
                }
                return func(temp_x);
            };
            std::pair<double, double> interval = golden_search(line_search_func, 0.0, e);
            double l = (interval.first + interval.second) / 2;
            
            for (int i = 0; i < x.size(); i++) {
                new_x[i] = x[i] - l * grad[i];
            }
        } 
        else {
            for (int i = 0; i < x.size(); i++) {
                new_x[i] = x[i] - grad[i];
            }
        }
        if(i % 11 == 0){
            if(func(new_x) >= func(test_convergence_point)){
                std::cout << "Num of func evaluations:" << eval_counter << "\n";
                eval_counter = 0;
                std::cout << "Num of gradient calculations:" << grad_counter << "\n";
                grad_counter = 0;

                std::cout << "Gradient descent doesnt converge, exiting function\n";
                std::cout << "Calculated point might not be function minimum!\n";
                return new_x;
            }
            test_convergence_point = new_x;
        }
        i++;
        x = new_x;

    } while (!grad_desc_exit_condition(grad, e));
    std::cout << "Num of func evaluations:" << eval_counter << "\n";
    std::cout << "Num of gradient calculations:" << grad_counter << "\n";
    eval_counter = 0;
    return x;
};

static inline bool optimize::grad_desc_exit_condition(std::vector<double> &grad,
                                                      const double e){
    return vector_norm(grad) < e;
}

static inline double optimize::vector_norm(std::vector<double> &vec){
    double res = 0;
    for(int i = 0; i < vec.size(); i++)
        res += (vec[i] * vec[i]);
    res = sqrtf(res);

    return res;
}

std::vector<double> optimize::newton_raphson(std::function<double(std::vector<double>)> func,
                                  std::vector<std::function<double(std::vector<double>)>> &part_deriv,
                                  std::vector<std::vector<std::function<double(std::vector<double>)>>> &hesse_matrix,
                                  const std::vector<double> &starting_point,
                                  const double e,
                                  const bool golden_ratio_used){

    std::vector<double> x = starting_point;
    std::vector<double> test_convergence_point = starting_point;
    int i = 1;
    Matrix grad(x.size(), 1);
    Matrix hesse{x.size(), x.size()};
    Matrix delta{x.size(), 1};
    std::vector<double> new_x(x.size());
    eval_counter = 0;
    grad_counter = 0;
    int hesse_counter = 0; 
    do {
        for (int i = 0; i < x.size(); i++) {
            grad(i, 0) = part_deriv[i](x);
            for(int j=0; j < x.size(); j++){
                hesse(i,j) = hesse_matrix[i][j](x);
            }
        }
        grad_counter++;
        hesse_counter++;
        delta = hesse.inverse() * grad; 
        delta *= -1;

        if (golden_ratio_used) {
            auto line_search_func = [&](double lambda) {
                std::vector<double> temp_x = x;
                for (int i = 0; i < x.size(); i++) {
                    temp_x[i] -= lambda * delta(i,0);
                }
                return func(temp_x);
            };
            std::pair<double, double> interval = golden_search(line_search_func, 0.0, e);
            double l = (interval.first + interval.second) / 2;
            
            for (int i = 0; i < x.size(); i++) {
                new_x[i] = x[i] - l * delta(i, 0);
            }
        } 
        else {
            for (int i = 0; i < x.size(); i++) {
                new_x[i] = x[i] + delta(i, 0);
            }
        }
        if(i % 11 == 0){
            if(func(new_x) >= func(test_convergence_point)){
                std::cout << "Num of func evaluations:" << eval_counter << "\n";
                eval_counter = 0;
                std::cout << "Num of gradient calculations:" << grad_counter << "\n";
                grad_counter = 0;
                std::cout << "Num of hesse calculations:" << hesse_counter << "\n";
                hesse_counter = 0;



                std::cout << "Newton - Raphson doesnt converge, exiting function\n";
                std::cout << "Calculated point might not be function minimum!\n";
                return new_x;
            }
            test_convergence_point = new_x;
        }
        i++;

        x = new_x;

    } while (!new_rap_exit_condition(delta, e));
    std::cout << "Num of func evaluations:" << eval_counter << "\n";
    eval_counter = 0;
    std::cout << "Num of gradient calculations:" << grad_counter << "\n";
    grad_counter = 0;
    std::cout << "Num of hesse calculations:" << hesse_counter << "\n";
    hesse_counter = 0;

    return x;
};

static inline bool optimize::new_rap_exit_condition(Matrix delta, const double e){
    std::vector<double> v(delta.getColumns() * delta.getRows());

    for(int i = 0;i<v.size();i++){
        v[i] = delta(i, 0);
    }

    return vector_norm(v) < e;
}

std::vector<double> optimize::gauss_newton(
                                  std::vector<std::function<double(std::vector<double>)>> &funcs,
                                  std::vector<std::vector<std::function<double(std::vector<double>)>>> &jacobian_matrix,
                                  const std::vector<double> &starting_point,
                                  const double e,
                                  const bool golden_ratio_used){
    std::vector<double> x = starting_point;
    std::vector<double> test_convergence_point = starting_point;
    int i = 1;
    Matrix jacobian{x.size(), x.size()};
    Matrix delta{x.size(), 1};
    Matrix G{x.size(), 1};
    std::vector<double> new_x(x.size());
    eval_counter = 0;
    grad_counter = 0;

    do {
        for (int i = 0; i < x.size(); i++) {
            for(int j=0; j < x.size(); j++){
                jacobian(i,j) = jacobian_matrix[i][j](x);
            }
            G(i,0) = funcs[i](x);
        }
        grad_counter++;

        delta = jacobian.inverse() * G; 
        delta *= -1;

        if (golden_ratio_used) {  
            auto line_search_func = [&](double lambda) {
            std::vector<double> temp_x = x;
            for (int i = 0; i < x.size(); i++) {
                temp_x[i] -= lambda * delta(i, 0);
            }

            double func_value = 0.0;
            for (int i = 0; i < funcs.size(); i++) {
                func_value += std::pow(funcs[i](temp_x), 2);
            }
            return func_value;
        };
            std::pair<double, double> interval = golden_search(line_search_func, 0.0, 1.0);
            double lambda = (interval.first + interval.second) / 2;

            for (int i = 0; i < x.size(); i++) {
                new_x[i] = x[i] - lambda * delta(i, 0);
            }
        } 
        else {
            for (int i = 0; i < x.size(); i++) {
                new_x[i] = x[i] + delta(i, 0);
            }
        }
        if(i % 11 == 0){
            if(pow(funcs[1](new_x),2) + pow(funcs[0](new_x), 2) >= pow(funcs[1](test_convergence_point),2) + pow(funcs[0](test_convergence_point),2)){
                std::cout << "Num of func evaluations:" << eval_counter << "\n";
                eval_counter = 0;
                std::cout << "Num of gradient calculations:" << grad_counter << "\n";
                grad_counter = 0;
        
                std::cout << "Gauss newton doesnt converge, exiting function\n";
                std::cout << "Calculated point might not be function minimum!\n";
                return new_x;
            }
            test_convergence_point = new_x;
        }
        i++;

        x = new_x;
    } while (!new_rap_exit_condition(delta, e)); // its the same exit condition as previous algo
                
    std::cout << "Num of func evaluations:" << eval_counter << "\n";
    eval_counter = 0;
    std::cout << "Num of gradient calculations:" << grad_counter << "\n";
    grad_counter = 0;
        
    return x;

};
