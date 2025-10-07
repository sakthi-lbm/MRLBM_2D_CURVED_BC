#ifndef POST_STATS_H
#define POST_STATS_H

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<cmath>

const double PI = 4.0 * atan(1.0);

struct CFDdata{
    double Str;
    double Cd_avg;
    double Cl_avg;
    double Cps;
    double Cpb;

    std::vector<int> idx;
    std::vector<int> time;
    std::vector<double> ntime;
    std::vector<double> theta;
    std::vector<double> F_drag;
    std::vector<double> F_lift;
    std::vector<double> Cd;
    std::vector<double> Cl;
    std::vector<double> Cp;
    std::vector<double> Cd_smooth;
    std::vector<double> Cl_smooth;
    std::vector<double> Cp_smooth;
    std::vector<double> ps_avg;
    std::vector<std::vector<double>> ps;
    std::vector<double> theta_sort;
    std::vector<double> Cp_sort;
    std::vector<double> Cp_smooth_sort;

};

[[nodiscard]] inline CFDdata allocate_vectors(const int nsteps, const int npoints){
    CFDdata data;
    data.time = std::vector<int>(nsteps);
    data.ntime = std::vector<double>(nsteps);
    data.F_drag = std::vector<double>(nsteps); 
    data.F_lift = std::vector<double>(nsteps);
    data.Cd = std::vector<double>(nsteps); 
    data.Cl = std::vector<double>(nsteps); 
    data.Cd_smooth = std::vector<double>(nsteps); 
    data.Cl_smooth = std::vector<double>(nsteps); 
    
    data.ps = std::vector<std::vector<double>>(nsteps, std::vector<double>(npoints));
    data.theta = std::vector<double>(npoints);
    data.Cp = std::vector<double>(npoints); 
    data.Cp_smooth = std::vector<double>(npoints);  
    data.ps_avg = std::vector<double>(npoints); 
    data.idx = std::vector<int>(npoints); 
    data.theta_sort = std::vector<double>(npoints);
    data.Cp_sort = std::vector<double>(npoints);  
    data.Cp_smooth_sort = std::vector<double>(npoints);  

    return data;

} 

struct cycleAndPeriod
{
    const int start;
    const int period;
};

[[nodiscard]] inline cycleAndPeriod find_peaks(const std::vector<double>& signal){
    int peak1 = -1;
    int peak2 = -1;
    double signal_avg = 0.0;

    for(int i = 0; i<signal.size(); i++){
        signal_avg += signal[i];
    }
    signal_avg /= signal.size();
    std:: cout << "Average of the signal is: " << signal_avg << std::endl;


    for (int i = 2; i < signal.size() - 2; i++) {
        if (signal[i] >= signal[i - 2] && signal[i] >= signal[i - 1] && signal[i] >= signal[i + 1] && signal[i] >= signal[i + 2] && signal[i] > signal_avg) {
            if (peak1 == -1) {
                peak1 = i;
            } else if (peak2 == -1) {
                peak2 = i;
                break;  // exit loop after finding second peak
            }
        }
    }

    if (peak1 == -1 || peak2 == -1) {
        throw std::runtime_error("Error: Could not detect two peaks.");
    }

   
    return {peak1, peak2 - peak1};


}


inline void gaussian_smoothing(const std::vector<double>& sig_x,
                                             const std::vector<double>& var, 
                                             std::vector<double>& var_smooth, const double sigma)                                              
{
    const int nl = sig_x.size();
    for(int i = 0; i < nl; i++){
        double sum_weights = 0.0;
        double sum_values = 0.0;

        for(int j = 0; j < nl; j++){
            double sigx2 = (sig_x[i] - sig_x[j]) * (sig_x[i] - sig_x[j]);
            double weight = exp(-sigx2/(2.0*sigma*sigma));
            sum_weights += weight;
            sum_values += weight * var[j];
        }

        var_smooth[i] = sum_values/sum_weights;
    }
}

inline void moving_average_smoothing(const std::vector<double>& sig_x,
                                             const std::vector<double>& var, 
                                             std::vector<double>& var_smooth, const int window_size)                                              
{
    const int nl = sig_x.size();
     int half = window_size / 2;
    for (int i = 0; i < nl; i++) {
        double sum = 0.0;
        int count = 0;
        for (int j = -half; j <= half; j++) {
            int idx = i + j;
            if (idx >= 0 && idx < nl) {
                sum += var[idx];
                count++;
            }
        }
        var_smooth[i] = sum / count;
    }
}






#endif
