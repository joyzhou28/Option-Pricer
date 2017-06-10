#include <cmath>

typedef struct Conf{
	double value;
	double confs[2];
}Conf;

// sample mean
inline double Mean(double Arr[], int M){
	double mean = .0;
	for(int i=0; i<M; i++){ mean+= Arr[i]; }
	return mean/M;
}
// sample std
inline double Std(double Arr[], int M){
	double mean = Mean(Arr, M);
	double std = .0;
	for(int i=0; i<M; ++i){ std+= pow(Arr[i]-mean,2); }
	return sqrt(std/(M-1));
}

// for standard normal distribution - pdf
inline double norm_pdf(double x){
	double const Pi = acos(double(-1)); 
	return 1.0/sqrt(2*Pi)*exp(-pow(x,2)/2);
}

// for standard normal distribution - cdf
inline double norm_cdf(double x){ return 0.5*erfc(-x*sqrt(0.5)); }

inline double geoMean(double Arr[], int M){
	double mul = 1.0;
	for(int i=0; i<M; i++){ mul*= Arr[i]; }
	return pow(mul, 1.0/M);
}