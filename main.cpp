#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>

void fcn(double t, double lambda, double *X, double *f){
    f[1] = (1 - lambda * t * t - X[0] * X[0]);
    f[0] = X[1];
    // f[0] = X[1];
    // f[1] = (1 - X[0] * X[0]) * X[1] - X[0];

}


bool runge(int N, double t_0, double *X, double t_end, double tol, double lambda){
    std::ofstream fout1, fout2, fout3, fout4;
    fout1.open("a1.txt");
    fout2.open("second.txt");
    // fout3.open("third.txt");
    // fout4.open("fourth.txt");
    bool flag = true;
    double h = 0.01;
    int const size = N * 6;
    double K[size];
    int i = 0;
    double err_glob = 0;
    double tmp = 0;
    double fac = 0;
    double h_new;
    const double facmax = 2.5;   
    const double safety_factor = 0.8;  
    const double facmin = 0.3;
    const double min_step = 1e-6;  
    int N_accepted = 0;
    int N_rejected = 0;
    double X_tmp[N];
    double er_loc[N];
    double error = 0;
    while(t_0 < t_end){
        if (t_0 + h > t_end) {
            h = t_end - t_0;  
        }
        error = 0;
        fcn(t_0, lambda, X, K);
        for(i = 0; i < N; i++){
            X_tmp[i] = X[i] + K[i] * h / 2;
        }
        fcn(t_0 + 0.5 * h, lambda, X_tmp, K + N);
        for(i = 0; i < N; i++){
            X_tmp[i] = X[i] + h * (K[i] + K[i + N]) / 4;
        }
        fcn(t_0 + 0.5 * h, lambda, X_tmp, K + 2 * N);
        for(i = 0; i < N; i++){
            X_tmp[i] = X[i] + h * (- K[i + N] + 2 * K[i + 2 * N]);
        }
        fcn(t_0 + h, lambda, X_tmp, K + 3 * N);
        for(i = 0; i < N; i++){
            X_tmp[i] = X[i] + h * (7 * K[i] + 10 * K[i + N] + K[i + 3 * N]) / 27;
        }
        fcn(t_0 + 2 * h / 3, lambda, X_tmp, K + 4 * N);
        for(i = 0; i < N; i++){
            X_tmp[i] = X[i] + h * (28 * K[i] - 125 * K[i + N] + 546 * K[i + 2 * N] + 54 * K[i + 3 * N] - 378 * K[i + 4 * N]) / 625;
        }
        fcn(t_0 + h/5, lambda, X_tmp, K + 5 * N);
        for(i = 0; i < N; i++){
            X_tmp[i] = X[i] + h * (K[i] + 4 * K[i + 2 * N] + K[i + 3 * N]) / 6;
        }
        error = 0;
	for (size_t i = 0; i < N; i++) {
		er_loc[i] = - (1. / 336.) * h * (42. * K[i] + 224. * K[i + 2 * N] + 21. * K[i + 3 * N]
					- 162. * K[i + 4 * N] - 125. * K[i + 5 * N]);
                    error += er_loc[i] * er_loc[i];
	}
        error = sqrt(error);
        double l = 1 > 4 * X_tmp[0] * X_tmp[0] ? 1 : 4 * X_tmp[0] * X_tmp[0];
        tmp = error + tmp * exp(l * h);
        fac = fmax(fmin(pow(error / tol, (1. / 5.)) / safety_factor, facmax), facmin);
        h_new = h / fac;
        if(error <= tol){
            t_0 += h;
            for(i = 0; i < N; i++){
                X[i] = X_tmp[i];
            }
            err_glob = tmp;
            // fout1 << t_0 << " " << X[0] << "\n";
            // fout2 << t_0 << " " << X[1] << "\n";
            // fout3 << t_0 << " " << X[2] << "\n";
            // fout4 << t_0 << " " << X[3] << "\n";
            // std::cout << "t = " << t_0 << ", ";
            // for(i = 0; i < N; i++){
            //     std::cout << "X" << i << " = " << X[i] << ", ";
            // }
            // std::cout << "h = " << h << "\n";
            // std::cout << "t = " << t_0 << ", ";
            // for(i = 0; i < N; i++){
            //     std::cout << "X" << i << " = " << X[i] << ", ";
            // }
            // std::cout << "h = " << h << "\n";
            N_accepted++;
        }
        else{
            N_rejected++;
        }
        if (h < min_step) {
            printf("Шаг слишком мал, возможно ошибка!\n");
            flag = false;
	    break;
        }
        h = h_new;
    }   
    std::cout << "t = " << t_0 << ", ";
            for(i = 0; i < N; i++){
                std::cout << "X" << i << " = " << X[i] << ", ";
            }
            std::cout << "h = " << h << ", ";
            std::cout << "errglob = " << err_glob << "\n";
    std::cout << "N accepted = " << N_accepted << "\n";
    std::cout << "N rejected = " << N_rejected << "\n";
    fout1.close();
    fout2.close();
    // fout3.close();
    // fout4.close();
    return flag;
}

int main(){
    std::ofstream fout;
    fout.open("new1.txt");
    long double time;
    double t_0 = 0;
    int N = 2;
    double X[N];
    X[0] = 0;
    X[1] = 1;
    // X[2] = 0;
    // X[3] = sqrt(3);
    double t_end = 1;
    double tol = 2e-16;
    double lambda = 2.795;
    time = clock();
    int count = 0;
    // for(lambda = -100; lambda <= 100; lambda += 0.0001){
    //     std::cout << "lambda = " << lambda << ", ";
    //     if(!runge(N, t_0, X, t_end, tol, lambda)){
    //         count++;
    //         if(count >= 100){
    //             break;
    //         }
    //     }
    //     else{
    //         count = 0;
    //     }
    //     fout << lambda << " " << X[0] << "\n";
    //     X[0] = 0;
    //     X[1] = 1;
    // }

    runge(N, t_0, X, t_end, tol, 17.1878);
    time = (clock() - time) / CLOCKS_PER_SEC;
    std::cout << "time = " << time;
    fout.close();
    return 0;
}