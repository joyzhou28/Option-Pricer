#ifndef __Asian_CPP
#define __Asian_CPP

#include "method.h"
#include <random>

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

  	std::random_device rd;
    std::mt19937 gen(rd());
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
#endif
