#ifndef __Method_CPP
#define __Method_CPP

#include "method.h"
#include <cmath>
#include <random>

// black scholes
BlackScholes::BlackScholes(double _S, double _K, 
                            double _r, double _q, 
                            double _T, std::string _type):
	S(_S), K(_K), r(_r), q(_q), T(_T), type(_type){}


double BlackScholes::option_price(double sigma) const {
	double d1 = (log(S/K)+(r-q)*T)/(sigma*sqrt(T))+0.5*sigma*sqrt(T);
	double d2 = d1-sigma*sqrt(T);
	double value;
	if(this->type=="call"){
		value = S*exp(-q*T)*norm_cdf(d1)-K*exp(-r*T)*norm_cdf(d2);
	}else{
		value = -S*exp(-q*T)*norm_cdf(-d1)+K*exp(-r*T)*norm_cdf(-d2);
	}
	return value;
}

double BlackScholes::option_vega(double sigma) const {
	double d1 = (log(S/K)+(r-q)*T)/(sigma*sqrt(T))+0.5*sigma*sqrt(T);
	return S*exp(-q*T)*sqrt(T)*norm_pdf(d1);
}

// binomial
Binomial::Binomial(double _S, double _sigma, 
	double _r, double _T, double _K, 
	std::string _type, int _N):
	S(_S), sigma(_sigma), r(_r), 
	T(_T), K(_K), type(_type), N(_N){
		this->dt = T/N;
		this->u = exp(sigma*sqrt(dt)); 
		this->d = 1/u;
		this->p = (exp(r*T)-d)/(u-d);
		// std::cout << u << std::endl << d << std::endl << p << std::endl;
	}

Binomial::~Binomial(){}

double Binomial::OptionPrice(){
	//initial values at time T
	double* sp = new double[N+1];
	for (int i = 0; i <= N; ++i) {
		if(type == "call"){
			sp[i] = fmax(S*pow(u,2*i-N)-K, 0);
		}else{
			sp[i] = fmax(K-S*pow(u,2*i-N), 0);
		}
	}

	// move to earlier times
	for (int j = N - 1; j >= 0; --j) {
		for (int i = 0; i <= j; ++i) {
			// binomial value for European option
			// sp[i] = (1-p)*sp[i]+p*sp[i+1];    
			// for American option
			if(type=="call"){
				sp[i] = fmax((1-p)*sp[i]+p*sp[i+1],S*pow(u,2*i-j)-K);
			}else{
				sp[i] = fmax((1-p)*sp[i]+p*sp[i+1],K-S*pow(u,2*i-j));
			}
			
		}
	}

	double ret = sp[0];
	delete[] sp;
	// ret *= exp(-r * T);
	return ret;
}

// asian
double Asian::Geo(){
	double sigmasqT = pow(sigma,2)*T*(n+1)*(2*n+1)/(6*n*n);
	double muT = 0.5*sigmasqT + (r-0.5*pow(sigma,2))*T*(n+1)/(2*n);
	double d1 = (log(S/K) + (muT+0.5*sigmasqT))/(sqrt(sigmasqT));
	double d2 = d1 - sqrt(sigmasqT);

	if(type=="call"){
		return exp(-r*T)*(S*exp(muT)*norm_cdf(d1) - K*norm_cdf(d2));
	}else{
		return exp(-r*T)*(-S*exp(muT)*norm_cdf(-d1) + K*norm_cdf(-d2));
	}
}

void Asian::MonteCarlo(int _M){
	this->M = _M;
	double* arithPayoff = new double[M];
	double* geoPayoff = new double[M];
	double dt = T/n;

	std::default_random_engine gen;
  	// std::random_device rd;
    // std::mt19937 gen(rd());
	std::normal_distribution<> ND(.0,1.0); // declare the distribution
	ND(gen);

	double drift = exp(dt*(r-0.5*sigma*sigma));
	double phi;
	double Spath;
	double arithMean, geoMean;

	for(int i=0; i<M; i++){ // i-th simulation
		Spath = S; // start from S0
		arithMean = .0; geoMean = .0;
		for(int j=0; j<n; j++){ // n is obs times
			phi = ND(gen);
			Spath*= drift*exp(sigma*sqrt(dt)*phi);
			arithMean+= Spath;
			geoMean+= log(Spath);
		}
		arithMean/= n; 
		geoMean = exp(geoMean/n);
		if(this->type == "call"){
			arithPayoff[i] = exp(-r*T)*fmax(arithMean-K, 0);
			geoPayoff[i] = exp(-r*T)*fmax(geoMean-K, 0);
		}else{
			arithPayoff[i] = exp(-r*T)*fmax(K-arithMean, 0);
			geoPayoff[i] = exp(-r*T)*fmax(K-geoMean, 0);
		}
	}

	double std;
	if(isArith){
		MC.value = Mean(arithPayoff,M);
		std = Std(arithPayoff,M);
		MC.confs[0] = MC.value-1.96*std/sqrt(M);
		MC.confs[1] = MC.value+1.96*std/sqrt(M);
		// control variables
		double covXY = .0;
		for(int i=0; i<M; i++){
			covXY+= arithPayoff[i]*geoPayoff[i];
		}
		double theta = (covXY/M-MC.value*Mean(geoPayoff,M))/pow(Std(geoPayoff,M),2);
		double* Z = new double[M];
		double geo = this->Geo();
		for(int i=0; i<M; i++){
			Z[i] = arithPayoff[i]+theta*(geo-geoPayoff[i]);
		}
		CV.value = Mean(Z,M);
		double stdz = Std(Z,M);
		delete[] Z;
		CV.confs[0] = CV.value-1.96*stdz/sqrt(M);
		CV.confs[1] = CV.value+1.96*stdz/sqrt(M);
	}else{
		MC.value = Mean(geoPayoff,M);
		std = Std(geoPayoff,M);
		MC.confs[0] = MC.value-1.96*std/sqrt(M);
		MC.confs[1] = MC.value+1.96*std/sqrt(M);
	}

	delete[] arithPayoff;
	delete[] geoPayoff;
}

// basket
double Basket::sigmaBg(){
	double sig = .0;
	for(int i=0; i<n_option; i++){
		for (int j=0; j<n_option; j++){
			if(j == i){
				sig+= sigma[i]*sigma[j];
			}else{
				sig+= sigma[i]*sigma[j]*p[i][j];
			}
		}
	}
	return sqrt(sig)/n_option;
	// return (sqrt(pow(sigma[0],2)+pow(sigma[1],2)+2*sigma[0]*sigma[1]*0.5))/this->n_option;
}

double Basket::muBg(){
	double ssigsq = .0;
	for(int i=0; i<n_option; i++){
		ssigsq+= pow(sigma[i],2);
	}
	return r-0.5*ssigsq/n_option;
}

double Basket::Geo(){
	double sigBg = this->sigmaBg();
	double muBg = this->muBg()+0.5*pow(sigBg,2);
	double bg = geoMean(S,n_option);
	double d1 = (log(bg/K)+T*(muBg+0.5*pow(sigBg,2)))/(sqrt(T)*sigBg);
	double d2 = d1 - (sqrt(T)*sigBg);

	if(type=="call"){
		return exp(-r*T)*(bg*exp(muBg*T)*norm_cdf(d1) - K*norm_cdf(d2));
	}else{
		return exp(-r*T)*(-bg*exp(muBg*T)*norm_cdf(-d1) + K*norm_cdf(-d2));
	}
}

void Basket::MonteCarlo(int _M){
	std::default_random_engine gen;
	// std::random_device rd;
    // std::mt19937 gen(rd()); // different serial each time
	std::normal_distribution<> ND(.0,1.0); // declare the distribution
	ND(gen);

	this->M = _M;
	double* arithPayoff = new double[M];
	double* geoPayoff = new double[M];
	double sqrT = sqrt(T); double df = exp(-r*T);
	double* Spath = new double[n_option];
	double* drift = new double[n_option];
	double* rn = new double[n_option];
	for(int i=0; i<n_option; i++){
		drift[i] = exp(T*(r-0.5*pow(sigma[i],2)));
	}
	for(int i=0; i<M; i++){
		rn[0] = ND(gen);
		Spath[0] = S[0]*drift[0]*exp(sigma[0]*sqrT*rn[0]);
		for(int j=1; j<n_option; j++){
			rn[j] = .0;
			double spsq = .0;
			for(int k=0; k<j; k++){
				rn[j]+= p[k][j]*rn[k];
				spsq+= pow(p[k][j],2);
			}
			rn[j]+= sqrt(1-spsq)*ND(gen);
			Spath[j] = S[j]*drift[j]*exp(sigma[j]*sqrT*rn[j]);
		}

		if(type == "call"){
			arithPayoff[i] = df*fmax(0,Mean(Spath,n_option)-K);
			geoPayoff[i] = df*fmax(0,geoMean(Spath,n_option)-K);
		}else{
			arithPayoff[i] = df*fmax(0,K-Mean(Spath,n_option));
			geoPayoff[i] = df*fmax(0,K-geoMean(Spath,n_option));
		}

	}
	double std;
	if(isArith){
		MC.value = Mean(arithPayoff,M);
		std = Std(arithPayoff,M);
		MC.confs[0] = MC.value - 1.96*std/sqrt(M);
		MC.confs[1] = MC.value + 1.96*std/sqrt(M);

		// control variables
		double covXY = .0;
		for(int i=0; i<M; i++){
			covXY+= arithPayoff[i]*geoPayoff[i];
		}
		double theta = (covXY/M-MC.value*Mean(geoPayoff,M))/pow(Std(geoPayoff,M),2);
		double* Z = new double[M];
		double geo = this->Geo();
		for(int i=0; i<M; i++){
			Z[i] = arithPayoff[i]+theta*(geo-geoPayoff[i]);
		}
		CV.value = Mean(Z,M);
		double stdz = Std(Z,M);
		// std::cout<< geo << ", " << geoPayoff[0] << std::endl;
		delete[] Z;
		CV.confs[0] = CV.value-1.96*stdz/sqrt(M);
		CV.confs[1] = CV.value+1.96*stdz/sqrt(M);
	}else{
		MC.value = Mean(geoPayoff,M);
		std = Std(arithPayoff,M);
		MC.confs[0] = MC.value - 1.96*std/sqrt(M);
		MC.confs[1] = MC.value + 1.96*std/sqrt(M);
	}

	delete[] arithPayoff;
	delete[] geoPayoff;
	delete[] rn; 
	delete[] drift;
	delete[] Spath;
}
#endif