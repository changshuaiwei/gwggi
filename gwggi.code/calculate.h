#ifndef _CALCULATE_H_
#define _CALCULATE_H_

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;



namespace Stat_fuc
{
	double auc_frm_LR(vector<double> & LR, vector<int> & n_patient, vector<int> & n_health);
	double auc_frm_LR(vector<double> & LR, vector<int> & n_patient, vector<int> & n_health, double & var);
	double auc_frm_LR(vector<double> & indi_LR, vector<bool> & if_disease);
	double auc_frm_LR(vector<double> & indi_LR, vector<bool> & if_disease, double & var);
	double auc_frm_LR(vector<double> & indi_LR, vector< vector<int> > & compare_struct);
	double auc_frm_LR(vector<double> & indi_LR, vector< vector<int> > & compare_struct, double & var);

	void score_frm_LR(vector<double> & indi_LR, vector< vector<int> > & compare_struct, vector<double> & fam_score);


	void score_frm_LR(vector<double> & LR, vector< vector<int> > & patient_idx, vector< vector<int> > & health_idx, int tnp, int tnh, vector<double> & score);
	void score_frm_LR(vector<double> & LR, vector<bool> & if_disease, vector<double> & score);
	double auc_frm_score(vector<double> & score, vector<bool> & if_disease);
	double auc_frm_score(vector<double> & score, vector<bool> & if_disease, double & var);
	void aucVarCov_frm_score(vector<double> & pre_score, double pre_auc, vector<double> & score, vector<bool> & if_disease, double & auc, double & var, double & cov);

	int max_index(vector<double> &arr);
	int first_peak(vector<double> &arr);
	void btstrp_sampling(int n, vector<int> & rst, int & seed);

	double chi_sqrq(double k_sqr, double df);//return q_value, tail value(lower tail)
	double chi_sqrp(double k_sqr, double df);//return p_value
	double std_norm_p(double z);//return p_value, when z is positive, it's Cumulative distribution function
	double std_norm_p1(double z);//return p_value, higher one tail
	double std_norm_q(double p);//return q value
	double gammq(const double a, const double x);
	double gammp(const double a, const double x);
	void gser(double &gamser, const double a, const double x, double &gln);
	double gammln(const double xx);
	void gcf(double &gammcf, const double a, const double x, double &gln);

	void sampling(vector<int> & ori, int size, vector<int> & rst, int & seed);
	void sampling_wr(vector<int> & ori, int size, vector<int> & rst, int & seed, vector<int> & rst_residual);
	void ran_devide(vector<int> & ori, int fold, vector< vector<int> > & rst, int & seed, vector< vector<int> > & rst2);//random devide to equal size
	double ran1(int &idum);//seed should be negative value for the first time
	double ran2(int &idum);//seed should be negative value for the first time

	bool emprical_CI(vector<double> & arr, double & Upper, double & Lower, double pct);
	double median(vector<double> & arr);
	double mean(vector<double> & arr);
	void indexx(vector<double> &arr, vector<int> & indx);


	inline void SWAP(int &a, int &b)
	{
		int temp=a;
		a=b;
		b=temp;
	}

	struct Erf {
		static const int ncof=28;
		static const double cof[28];

		inline double SQR(double x){
			return x*x;
		}

		inline double erf(double x) {
			if (x >=0.) return 1.0 - erfccheb(x);
			else return erfccheb(-x) - 1.0;
		}

		inline double erfc(double x) {
			if (x >= 0.) return erfccheb(x);
			else return 2.0 - erfccheb(-x);
		}

		double erfccheb(double z){
			int j;
			double t,ty,tmp,d=0.,dd=0.;
			if (z < 0.) cerr<<"erfccheb requires nonnegative argument";
			t = 2./(2.+z);
			ty = 4.*t - 2.;
			for (j=ncof-1;j>0;j--) {
				tmp = d;
				d = ty*d - dd + cof[j];
				dd = tmp;
			}
			return t*exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
		}

		double inverfc(double p) {
			double x,err,t,pp;
			if (p >= 2.0) return -100.;
			if (p <= 0.0) return 100.;
			pp = (p < 1.0)? p : 2. - p;
			t = sqrt(-2.*log(pp/2.));
			x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
			for (int j=0;j<2;j++) {
				err = erfc(x) - pp;
				x += err/(1.12837916709551257*exp(-(x*x))-x*err);
			}
			return (p < 1.0? x : -x);
		}

		inline double inverf(double p) {return inverfc(1.-p);}

		double erfcc(const double x)
		{
			double t,z=fabs(x),ans;
			t=2./(2.+z);
			ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
				t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
				t*(-0.82215223+t*0.17087277)))))))));
			return (x >= 0.0 ? ans : 2.0-ans);
		}


	};


	struct Normaldist : Erf {
		double mu, sig;
		Normaldist(double mmu = 0., double ssig = 1.) : mu(mmu), sig(ssig) {
			if (sig <= 0.) cerr<<"bad sig in Normaldist";
		}
		double p(double x) {
			return (0.398942280401432678/sig)*exp(-0.5*SQR((x-mu)/sig));
		}
		double cdf(double x) {
			return 0.5*erfc(-0.707106781186547524*(x-mu)/sig);
		}
		double invcdf(double p) {
			if (p <= 0. || p >= 1.) cerr<<"bad p in Normaldist";
			return -1.41421356237309505*sig*inverfc(2.*p)+mu;
		}
	};





}


namespace Mat_fuc
{
	void svdcmp(vector< vector<double> > &a, vector<double> &w, vector< vector<double> > &v);//a=u*w*v(T),a is changed to u
	double SVD_inverse1(vector< vector<double> > & mat, vector< vector<double> > & inverse);//return rank, normal index(same as used in math)
	double SVD_inverse2(vector< vector<double> > & mat, vector< vector<double> > & inverse);//transfered index
	double pythag(const double a, const double b);
	void svbksb(vector< vector<double> > &u, vector<double> &w, vector< vector<double> > &v, vector<double> &b, vector<double> &x);
	void solve(vector< vector<double> > &X, vector<double> &b, vector<double> &Y);

	inline double SIGN(double a, double b);
	inline double MAX(double a, double b);
	inline double MIN(double a, double b);
	inline double SQR(double a);

	inline int MAX(int a, int b);
	inline int MIN(int a, int b);

}

inline double Mat_fuc::SIGN(double a, double b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline double Mat_fuc::MAX(double a, double b)
{
	return b > a ? (b) : (a);
}

inline double Mat_fuc::MIN(double a, double b) 
{
	return b < a ? (b) : (a);
}

inline double Mat_fuc::SQR(double a) 
{
	return a*a;
}

inline int Mat_fuc::MAX(int a, int b)
{
	return b > a ? (b) : (a);
}

inline int Mat_fuc::MIN(int a, int b) 
{
	return b < a ? (b) : (a);
}



#endif