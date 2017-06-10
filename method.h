#ifndef __Method_H
#define __Method_H

#include <string>
#include <iostream>
#include "common.h"

// balck_scholes
class BlackScholes{
	private:
		double S; //S(0)
		double K; //strike price
		double r; // risk-free interest rate
		double q; // repo rate
		double T; // time to maturity
		std::string type;
	public:
		BlackScholes(double _S, double _K, 
                   double _r, double _q,
                   double _T, std::string _type);
		double option_price(double sigma) const;
		double option_vega(double sigma) const;
		~BlackScholes(){}
};

// binomial
class Binomial{
	private:
		double S;
		double sigma;
		double r;
		double T;
		double K;
		std::string type;
		int N;
		double dt, u, d, p;
	public:
		Binomial(double _S, double _sigma, double _r, double _T, double _K, std::string _type, int _N);
		~Binomial();
		double OptionPrice();
		void Display(){
			std::cout << u << "," << d << "," << p << std::endl;
		}
};

//asian
class Asian{
protected:
	double S,sigma, r,T,K; 
	int n; // obs times
	std::string type; // "call" or "put"
	int M; // path#
	bool isArith;
	Conf MC, CV;
public:
	Asian(double _S, double _sigma, double _r, double _T, double _K, int _n, std::string _type):S(_S), sigma(_sigma), r(_r), T(_T), K(_K), n(_n), type(_type){ isArith = 0; }
	double Geo();
	void MonteCarlo(int _M);
	Conf getMC(){ return this->MC; }
	~Asian(){};
};

class Arith_Asian : public Asian{
public:
	Arith_Asian(double _S, double _sigma, double _r, double _T, double _K, int _n, std::string _type):Asian(_S, _sigma, _r, _T, _K, _n, _type){ isArith = 1;}
	Conf getCV(){ return this->CV; }
	~Arith_Asian(){}
};

// basket
class Basket{
protected:
	double* S;
	double* sigma; 
	double** p;
	double r,T,K; 
	int n_option; // # of options
	std::string type; // "call" or "put"
	int M; // path#
	bool isArith;
	Conf MC, CV;
public:
	Basket(double* _S, double* _sigma, double** _p, double _r, double _T, double _K, int _n_option, std::string _type): S(_S), sigma(_sigma), p(_p), r(_r), T(_T), K(_K), n_option(_n_option), type(_type){ isArith = 0;}

	double sigmaBg();
	double muBg();
	double Geo();
	void MonteCarlo(int _M);
	Conf getMC(){ return this->MC; }
	~Basket(){}
};

class Arith_Basket : public Basket{
public:
	Arith_Basket(double* _S, double* _sigma, double** _p, double _r, double _T, double _K, int _n_option, std::string _type): Basket(_S, _sigma, _p, _r, _T, _K, _n_option, _type){ isArith = 1;}
	Conf getCV(){ return this->CV; }
	~Arith_Basket(){}
};

#endif