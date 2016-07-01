#include "calculate.h"

double Stat_fuc::gammq(const double a, const double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0)
		cerr<<"Invalid arguments in routine gammq";
	if (x < a+1.0) {
		gser(gamser,a,x,gln);
		return 1.0-gamser;
	} else {
		gcf(gammcf,a,x,gln);
		return gammcf;
	}
}

void Stat_fuc::gser(double &gamser, const double a, const double x, double &gln)
{
	const int ITMAX=100;
	const double EPS=numeric_limits<double>::epsilon();
	int n;
	double sum,del,ap;

	gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) cerr<<"x less than 0 in routine gser";
		gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=0;n<ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				gamser=sum*exp(-x+a*log(x)-gln);
				return;
			}
		}
		cerr<<"a too large, ITMAX too small in routine gser";
		return;
	}
}

double Stat_fuc::gammln(const double xx)
{
	int j;
	double x,y,tmp,ser;
	static const double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
		-0.5395239384953e-5};

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<6;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

void Stat_fuc::gcf(double &gammcf, const double a, const double x, double &gln)
{
	const int ITMAX=100;
	const double EPS=numeric_limits<double>::epsilon();
	const double FPMIN=numeric_limits<double>::min()/EPS;
	int i;
	double an,b,c,d,del,h;

	gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) <= EPS) break;
	}
	if (i > ITMAX) cerr<<"a too large, ITMAX too small in gcf";
	gammcf=exp(-x+a*log(x)-gln)*h;
}

double Stat_fuc::chi_sqrq(double k_sqr, double df)
{
	return (gammq(df*0.5, k_sqr*0.5));
}

double Stat_fuc::std_norm_p(double z)
{
	if(z>=0) return (0.5+0.5*gammp(0.5, 0.5*z*z));
	else return (0.5-0.5*gammp(0.5, 0.5*z*z));
}

double Stat_fuc::std_norm_p1(double z)
{
	if(z>=0) return (0.5-0.5*gammp(0.5, 0.5*z*z));
	else return (0.5+0.5*gammp(0.5, 0.5*z*z));
}

double Stat_fuc::std_norm_q(double p)
{
	Normaldist norm;
	norm.mu=0; norm.sig=1;
	return norm.invcdf(p);
}


double Stat_fuc::chi_sqrp(double k_sqr, double df)
{
	return (gammp(df*0.5, k_sqr*0.5));
}

void Stat_fuc::indexx(vector<double> &arr, vector<int> &indx)
{
	const int M=7,NSTACK=50;
	int i,indxt,ir,j,k,jstack=-1,l=0;
	double a;
	vector<int> istack(NSTACK);

	int n=arr.size();
	ir=n-1;
	indx.resize(n);
	for (j=0;j<n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=l;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack < 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir]);
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir]);
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1]);
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j]);
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack >= NSTACK) cerr<<"NSTACK too small in indexx.";
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}



bool Stat_fuc::emprical_CI(vector<double> & arr, double & Upper, double & Lower, double pct)
{
	if(arr.size()==0) return false;
	else{
		vector<int> idx;
		indexx(arr, idx);

		double size=double(arr.size());
		int ipct=int((size-1)*(1-pct)*0.5);
		int ulimit=arr.size()-1-ipct; int llimit=ipct;

		bool U_mk=false; bool L_mk=false;
		for(int i=0; i<idx.size(); i++){
			if(i==llimit && (!L_mk)) { Lower=arr[idx[i]]; L_mk=true;}
			if(i==ulimit && (!U_mk)) { Upper=arr[idx[i]]; U_mk=true;}
			if(L_mk && U_mk) break;
		}

		if( (!L_mk) || (!U_mk) ) return false;
	}

	return true;
}

double Stat_fuc::median(vector<double> & arr)
{
	double med=0;
	if(arr.size()==0) return 0;
	else if(arr.size()==1) { med=arr[0];}
	else if(arr.size()==2) {med=(arr[0]+arr[1])/2; }
	else{
		vector<int> idx;
		indexx(arr,idx);
		int size=arr.size();
		if(size%2==0) med=(arr[(size/2)-1]+arr[(size/2)])*0.5;
		else med=arr[(size-1)/2];
	}
	return med;
}

double Stat_fuc::mean(vector<double> & arr)
{
	double tmp=0;
	for(int i=0; i<arr.size(); i++){
		tmp+=arr[i];
	}
	if(arr.size()>0) tmp/=double(arr.size());

	return tmp;
}

double Stat_fuc::gammp(const double a, const double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0)
		cerr<<"Invalid arguments in routine gammp";
	if (x < a+1.0) {
		gser(gamser,a,x,gln);
		return gamser;
	} else {
		gcf(gammcf,a,x,gln);
		return 1.0-gammcf;
	}
}

double Stat_fuc::ran1(int &idum)
{
	const int IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32;
	const int NDIV=(1+(IM-1)/NTAB);
	const double EPS=3.0e-16,AM=1.0/IM,RNMX=(1.0-EPS);
	static int iy=0;
	static vector<int> iv(NTAB);
	int j,k;
	double temp;

	if (idum <= 0 || !iy) {
		if (-idum < 1) idum=1;
		else idum = -idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ;
			idum=IA*(idum-k*IQ)-IR*k;
			if (idum < 0) idum += IM;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

double Stat_fuc::ran2(int &idum)
{
	const int IM1=2147483563,IM2=2147483399;
	const int IA1=40014,IA2=40692,IQ1=53668,IQ2=52774;
	const int IR1=12211,IR2=3791,NTAB=32,IMM1=IM1-1;
	const int NDIV=1+IMM1/NTAB;
	const double EPS=3.0e-16,RNMX=1.0-EPS,AM=1.0/double(IM1);
	static int idum2=123456789,iy=0;
	static vector<int> iv(NTAB);
	int j,k;
	double temp;

	if (idum <= 0) {
		idum=(idum==0 ? 1 : -idum);
		idum2=idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ1;
			idum=IA1*(idum-k*IQ1)-k*IR1;
			if (idum < 0) idum += IM1;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ1;
	idum=IA1*(idum-k*IQ1)-k*IR1;
	if (idum < 0) idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


void Stat_fuc::sampling(vector<int> & ori, int size, vector<int> & rst, int & seed)
{
	rst.resize(size);
	int total=ori.size();
	int i=0;
	int indx=0;
	double dtotal=total;
	for(i=0; i<size; i++){
		indx=(int)(ran1(seed)*(dtotal));
		while (indx<0 || indx>(total-1)) indx=(int)(ran1(seed)*(dtotal));
		rst[i]=ori[indx];

	}
}

void Stat_fuc::sampling_wr(vector<int> & ori, int size, vector<int> & rst, int & seed, vector<int> & rst_residual)
{
	int total=ori.size();
	int i=0;
	int indx=0;
	//rst.resize(size); rst_residual.resize(total-size);

	vector<int> tmp; tmp.resize(total,0);

	double dtotal=total;
	int count=0;
	while(1){
		indx=(int)(ran1(seed)*(dtotal));
		while (indx<0 || indx>(total-1)) indx=(int)(ran1(seed)*(dtotal));
		if(tmp[indx]==0){
			tmp[indx]=1;
			count++;
		}
		if( count==size ) break;
	}
	rst.clear(); rst_residual.clear();
	for(int i=0; i<tmp.size(); i++)	{
		if(tmp[i]==1) rst.push_back(ori[i]);
		else rst_residual.push_back(ori[i]);
	}
}


void Stat_fuc::ran_devide(vector<int> & ori, int fold, vector< vector<int> > & rst, int & seed, vector< vector<int> > & rst2)
{
	rst.clear(); rst2.clear();
	rst.resize(fold); rst2.resize(fold);
	int i=0, j=0;
	int indx=0;
	double dfold=fold;
	for(i=0; i<ori.size(); i++){
		indx=(int)(ran1(seed)*(dfold));
		while (indx<0 || indx>(fold-1)) indx=(int)(ran1(seed)*(dfold));
		rst[indx].push_back(ori[i]);

		for(j=0; j<fold; j++) if(j!=indx) rst2[j].push_back(ori[i]);
	}
}

int Stat_fuc::first_peak(vector<double> &arr)
{
	double temp=0, max=0;
	int i=0;
	int mark=0;
	for (i=0; i<arr.size(); i++)
	{
		temp=arr[i];
		if(temp>max){
			max=temp;
			mark=i;
		}
		if((i+1)!=arr.size()){
			if(max>arr[i+1]) break;
		}
	}

	return mark;

}

const double Stat_fuc::Erf::cof[28] = {-1.3026537197817094, 6.4196979235649026e-1,
1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
-1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
-1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17};


void Mat_fuc::svdcmp(vector< vector<double> > &a, vector<double> &w, vector< vector<double> > &v)
{
	bool flag;
	int i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z;

	int m=a.size();
	int n=a[0].size();
	vector<double> rv1(n);
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a[k][i]);
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i+1 != n) {
			for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l-1];
				g = -(SIGN(sqrt(s),f));
				h=f*g-s;
				a[i][l-1]=f-g;
				for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
					for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l-1;k<n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i][j]=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<m;j++) a[j][i] *= g;
		} else for (j=i;j<m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=true;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if (fabs(rv1[l])+anorm == anorm) {
					flag=false;
					break;
				}
				if (fabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if (fabs(f)+anorm == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 29) cerr<<"no convergence in 30 svdcmp iterations";
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
}

double Mat_fuc::pythag(const double a, const double b)
{
	double absa,absb;

	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

double Mat_fuc::SVD_inverse1(vector< vector<double> > & mat, vector< vector<double> > & inverse)
{
	vector< vector<double> > ori;
	ori=mat;
	int m=mat.size();
	int n=mat[0].size();
	vector< vector<double> > v;
	v.resize(n);
	int i=0, j=0, k=0;
	for(i=0; i<n; i++) v[i].resize(n);
	vector<double> w(n);
	svdcmp(ori, w, v);

	inverse.resize(n);
	for(i=0; i<n; i++) inverse[i].resize(m);
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			inverse[i][j]=0;
			for(k=0; k<n; k++){
				if(w[k]!=0) inverse[i][j]+=v[i][k]*(1.0/w[k])*ori[j][k];
			}
		}
	}

	double trace=0.0;
	for(i=0; i<n; i++)
	{
		if( w[i]!=0.0 ) trace+=1.0;
	}

	return trace;

}

double Mat_fuc::SVD_inverse2(vector< vector<double> > & mat, vector< vector<double> > & inverse)
{
	int i=0, j=0, k=0;

	//	ori=mat;
	int n=mat.size();
	int m=mat[0].size();
	vector< vector<double> > ori(m);

	for(i=0; i<m; i++) ori[i].resize(n);

	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			ori[i][j]=mat[j][i];
		}
	}

	vector< vector<double> > v;
	v.resize(n);

	for(i=0; i<n; i++) v[i].resize(n);
	vector<double> w(n);
	svdcmp(ori, w, v);

	inverse.resize(m);
	for(i=0; i<m; i++) inverse[i].resize(n);
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			inverse[i][j]=0;
			for(k=0; k<n; k++){
				if(w[k]!=0) inverse[i][j]+=v[j][k]*(1.0/w[k])*ori[i][k];
			}
		}
	}

	double trace=0.0;
	for(i=0; i<n; i++)
	{
		if( w[i]!=0.0 ) trace+=1.0;
	}

	return trace;

}

void Mat_fuc::svbksb(vector< vector<double> > &u, vector<double> &w, vector< vector<double> > &v, vector<double> &b, vector<double> &x)
{
	int jj,j,i;
	double s;

	int m=u.size();
	int n=u[0].size();
	x.resize(m);
	vector<double> tmp(n);
	for (j=0;j<n;j++) {
		s=0.0;
		if(w[j]!=0){
			for (i=0;i<n;i++) s += v[i][j]*b[i];
			s /= w[j];
		}

		tmp[j]=s;
	}
	for (j=0;j<m;j++) {
		s=0.0;
		for (jj=0;jj<n;jj++) s += u[j][jj]*tmp[jj];
		x[j]=s;
	}

	/*
	for (j=0;j<m;j++) {
	s=0.0;
	if (w[j] != 0.0) {
	for (i=0;i<m;i++) s += v[i][j]*b[i];
	s /= w[j];
	}
	tmp[j]=s;
	}
	for (j=0;j<n;j++) {
	s=0.0;
	for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
	x[j]=s;
	}
	*/
}

void Mat_fuc::solve(vector< vector<double> > &A, vector<double> &x, vector<double> &b)
{
	vector< vector<double> > ori;
	ori=A;
	int m=ori.size();
	int n=ori[0].size();
	vector< vector<double> > v;
	v.resize(n);
	int i=0, j=0, k=0;
	for(i=0; i<n; i++) v[i].resize(n);
	vector<double> w(n);
	svdcmp(ori, w, v);
	svbksb(ori,w,v,b,x);
}

double Stat_fuc::auc_frm_LR(vector<double> & LR, vector<int> & n_patient, vector<int> & n_health)
{
	vector<int> indx;
	indx.resize(LR.size());
	indexx(LR, indx);
	int np=0, nh=0;
	int i=0;
	double score=0, last_LR=LR[indx[0]];
	int cell_np=0, cell_nh=0;
	for(i=0; i<indx.size(); i++){
		int temp=indx[i];
		

		if(LR[indx[i]]!=last_LR){//in case that we have equal LR values
			if(cell_np!=0){
				score+=double(cell_np*nh);
				score+=0.5*(double(cell_np*cell_nh));
			}
			nh+=cell_nh; np+=cell_np;
			cell_nh=0; cell_np=0;
		}

		cell_nh+=n_health[temp]; cell_np+=n_patient[temp];
	}

	if(cell_np!=0){
		score+=double(cell_np*nh);
		score+=0.5*(double(cell_np*cell_nh));
	}
	nh+=cell_nh; np+=cell_np;
	cell_nh=0; cell_np=0;

	score/=((double)(nh*np));
	return score;
}

double Stat_fuc::auc_frm_LR(vector<double> & LR, vector<int> & n_patient, vector<int> & n_health, double & var)
{
	vector<int> indx;
	indx.resize(LR.size());
	indexx(LR, indx);
	double np=0, nh=0;
	int i=0, j=0;
	double score=0;

	vector<double> n_pa, n_he;
	double pa=0, he=0;
	double pre_lr=LR[indx[0]];

	for(i=0; i<indx.size(); i++){
		int temp=indx[i];
		if(LR[temp]>pre_lr){
			n_pa.push_back(pa);
			n_he.push_back(he);
			pa=n_patient[temp];
			he=n_health[temp];
			pre_lr=LR[temp];
		}else{
			pa+=n_patient[temp];
			he+=n_health[temp];
		}
	}
	n_pa.push_back(pa);
	n_he.push_back(he);

	for(i=0; i<n_pa.size(); i++){
		if(n_pa[i]!=0){
			score+=n_pa[i]*nh;
			score+=0.5*n_pa[i]*n_he[i];
		}
		nh+=n_he[i];
		np+=n_pa[i];
	}

	score/=((double)(nh*np));


	//cal auc var
	double t_p=double(np);
	double t_h=double(nh);
	vector<double> TPR, FPR;
	double temp_t=0, temp_f=0;

	for(i=n_pa.size()-1; i>=0; i--){

		temp_t+=(double(n_pa[i]))/t_p;
		TPR.push_back(temp_t);

		temp_f+=(double(n_he[i]))/t_h;
		FPR.push_back(temp_f);
	}

	double var_d=0, var_nd=0;

	for(i=0; i<TPR.size(); i++){
		if(i==0){
			var_d+=TPR[i]*TPR[i]*FPR[i]*0.5;
			var_nd+=TPR[i]*FPR[i]*FPR[i]*0.5;
		}else{
			var_d+=0.5*(TPR[i]*TPR[i]+TPR[i-1]*TPR[i-1])*(FPR[i]-FPR[i-1]);
			var_nd+=0.5*(FPR[i]*FPR[i]+FPR[i-1]*FPR[i-1])*(TPR[i]-TPR[i-1]);
		}
	}
	var_d-=score*score;
	var_nd-=(1-score)*(1-score);
	var=var_d/t_h+var_nd/t_p;

	return score;
}

int Stat_fuc::max_index(vector<double> &arr)
{
	double temp=0, max=0;
	int i=0;
	int mark=0;
	for (i=0; i<arr.size(); i++)
	{
		temp=arr[i];
		if(temp>max){
			max=temp;
			mark=i;
		}
	}

	return mark;

}


void Stat_fuc::btstrp_sampling(int n, vector<int> & rst, int & seed)
{
	int i, indx;
	rst.clear(); rst.resize(n,0);
	double dn=double(n);

	i=0;
	while(i<n){
		indx=int(dn*(Stat_fuc::ran1(seed)));
		if(indx>=0 && indx<n){
			rst[indx]++;
			i++;
		}
	}

}
void Stat_fuc::score_frm_LR(vector<double> & LR, vector< vector<int> > & patient_idx, vector< vector<int> > & health_idx, int tnp, int tnh, vector<double> & score)
{
	int i=0, j=0, k=0;
	vector<int> rindx;
	rindx.resize(LR.size()); 
	score.clear(); score.resize((tnh+tnp));

	indexx(LR, rindx);
	double cur_np=0, cur_nh=0;
	double pre_lr=LR[rindx[0]];
	vector<int> tmp_idx; vector< vector<int> > vec_idx;
	vector<double> vec_np, vec_nh;


	for(i=0; i<rindx.size(); i++){
		if(LR[rindx[i]]>pre_lr){
			pre_lr=LR[rindx[i]];

			vec_idx.push_back(tmp_idx);
			tmp_idx.clear(); tmp_idx.push_back(rindx[i]);
			vec_np.push_back(cur_np); vec_nh.push_back(cur_nh);
			cur_nh=health_idx[rindx[i]].size();
			cur_np=patient_idx[rindx[i]].size();
		}else{
			cur_nh+=health_idx[rindx[i]].size();
			cur_np+=patient_idx[rindx[i]].size();
			tmp_idx.push_back(rindx[i]);
		}
	}

	vec_idx.push_back(tmp_idx);
	vec_np.push_back(cur_np); vec_nh.push_back(cur_nh);

	double np=tnp, nh=0, dtnp=tnp, dtnh=tnh;

	double pscore=0.0, hscore=0.0;
	for(i=0; i<vec_idx.size(); i++){

		if(vec_np[i]!=0){ 
			np-=vec_np[i];
			pscore=nh+0.5*vec_nh[i];
			pscore/=dtnh;
			for(j=0; j<vec_idx[i].size(); j++){
				for(k=0;k<patient_idx[vec_idx[i][j]].size(); k++) score[patient_idx[vec_idx[i][j]][k]]=pscore;
			}
		}

		if(vec_nh[i]!=0){
			nh+=vec_nh[i];
			hscore=np+0.5*vec_np[i];
			hscore/=dtnp;
			for(j=0; j<vec_idx[i].size(); j++){
				for(k=0;k<health_idx[vec_idx[i][j]].size(); k++) score[health_idx[vec_idx[i][j]][k]]=hscore;
			}
		}
	}
}
void Stat_fuc::score_frm_LR(vector<double> & LR, vector<bool> & if_disease, vector<double> & score)
{
	int i=0, j=0, k=0;
	vector<int> rindx;
	score.clear(); score.resize(if_disease.size());

	Stat_fuc::indexx(LR, rindx);
	double cur_np=0, cur_nh=0;
	double pre_lr=LR[rindx[0]];
	vector<int> tmp_idx; vector< vector<int> > vec_idx;
	vector<double> vec_np, vec_nh;


	for(i=0; i<rindx.size(); i++){
		if(LR[rindx[i]]>pre_lr){
			pre_lr=LR[rindx[i]];

			vec_idx.push_back(tmp_idx);
			tmp_idx.clear(); 
			vec_np.push_back(cur_np); vec_nh.push_back(cur_nh);

			cur_nh=0; cur_np=0;
		}

		if(if_disease[rindx[i]]) cur_np+=1; else cur_nh+=1;
		tmp_idx.push_back(rindx[i]);
	}

	vec_idx.push_back(tmp_idx);
	vec_np.push_back(cur_np); vec_nh.push_back(cur_nh);

	double np=0, nh=0, dtnp=0, dtnh=0;

	for(i=0; i<if_disease.size(); i++){
		if(if_disease[i]) dtnp+=1; else dtnh+=1;
	}

	np=dtnp;

	double pscore=0.0, hscore=0.0;
	for(i=0; i<vec_idx.size(); i++){

		if(vec_np[i]!=0){ 
			np-=vec_np[i];
			pscore=nh+0.5*vec_nh[i];
			pscore/=dtnh;
			for(j=0; j<vec_idx[i].size(); j++){
				if(if_disease[ vec_idx[i][j] ]) score[vec_idx[i][j]]=pscore;
			}
		}

		if(vec_nh[i]!=0){
			nh+=vec_nh[i];
			hscore=np+0.5*vec_np[i];
			hscore/=dtnp;
			for(j=0; j<vec_idx[i].size(); j++){
				if(!if_disease[vec_idx[i][j]]) score[vec_idx[i][j]]=hscore;
			}
		}
	}
}
double Stat_fuc::auc_frm_score(vector<double> & score, vector<bool> & if_disease)
{
	double auc=0.0, nd=0;
	int i,j;
	for(i=0; i<score.size(); i++){
		if(if_disease[i]){
			auc+=score[i]; nd+=1;
		}
	}
	auc/=nd;
	return auc;
}
double Stat_fuc::auc_frm_score(vector<double> & score, vector<bool> & if_disease, double & var)
{
	double auc=0.0, nd=0, nh=0;
	int i,j;
	for(i=0; i<score.size(); i++){
		if(if_disease[i]){
			auc+=score[i]; nd+=1;
		}else{
			nh+=1;
		}
	}
	auc/=nd;

	vector<double> bias;
	bias=score;

	for(i=0; i<bias.size(); i++) bias[i]-=auc;

	double pSD=0, hSD=0;
	for(i=0; i<bias.size(); i++){
		if(if_disease[i]){
			pSD+=bias[i]*bias[i];
		}else{
			hSD+=bias[i]*bias[i];
		}
	}
	pSD/=(nd*nd);
	hSD/=(nh*nh);
	var=pSD+hSD;

	return auc;
}

void Stat_fuc::aucVarCov_frm_score(vector<double> & pre_score, double pre_auc, vector<double> & score, vector<bool> & if_disease, double & auc, double & var, double & cov)
{
	double nd=0, nh=0;
	int i,j;
	auc=0;
	for(i=0; i<score.size(); i++){
		if(if_disease[i]){
			auc+=score[i]; nd+=1;
		}else{
			nh+=1;
		}
	}
	auc/=nd;

	vector<double> bias, pre_bias;
	bias=score; pre_bias=pre_score;

	for(i=0; i<bias.size(); i++) bias[i]-=auc;
	for(i=0; i<pre_bias.size(); i++) pre_bias[i]-=pre_auc;

	double pSD=0, hSD=0;
	for(i=0; i<bias.size(); i++){
		if(if_disease[i]){
			pSD+=bias[i]*bias[i];
		}else{
			hSD+=bias[i]*bias[i];
		}
	}
	pSD/=(nd*nd);
	hSD/=(nh*nh);
	var=pSD+hSD;

	double pCD=0, hCD=0;
	for(i=0; i<if_disease.size(); i++){
		if(if_disease[i]){
			pCD+=pre_bias[i]*bias[i];
		}else{
			hCD+=pre_bias[i]*bias[i];
		}
	}
	pCD/=(nd*nd);
	hCD/=(nh*nh);

	cov=pCD+hCD;
}

double Stat_fuc::auc_frm_LR(vector<double> & indi_LR, vector<bool> & if_disease)
{
	vector<double> score;
	score_frm_LR(indi_LR,if_disease,score);
	double auc=auc_frm_score(score,if_disease);
	return auc;
}

double Stat_fuc::auc_frm_LR(vector<double> & indi_LR, vector<bool> & if_disease, double & var)
{
	vector<double> score;
	score_frm_LR(indi_LR,if_disease,score);
	double auc=auc_frm_score(score,if_disease,var);
	return auc;
}

double Stat_fuc::auc_frm_LR(vector<double> & indi_LR, vector< vector<int> > & compare_struct)
{
	int n_fam=compare_struct.size()/2;
	int i, j, k, idx_p, idx_h;
	double LR_p, LR_h, tmp_n_cmp=0, tmp_score=0, score=0;
	vector<double> fam_score;
	for(i=0; i<n_fam; i++){
		idx_p=i*2; idx_h=i*2+1;
		tmp_score=0; tmp_n_cmp=0;
		for(j=0; j<compare_struct[idx_p].size(); j++){
			LR_p=indi_LR[compare_struct[idx_p][j]];
			for(k=0; k<compare_struct[idx_h].size(); k++){
				LR_h=indi_LR[compare_struct[idx_h][k]];
				tmp_n_cmp+=1;
				if(LR_p>LR_h) tmp_score+=1;
				else if(LR_p==LR_h) tmp_score+=0.5;
			}
		}
		if(tmp_score!=0) tmp_score/=tmp_n_cmp;
		if(tmp_n_cmp!=0) fam_score.push_back(tmp_score);
		//fam_score.push_back(tmp_score);
	}

	for(i=0; i<fam_score.size(); i++) score+=fam_score[i];
	double fam_num=double(fam_score.size());
	score/=fam_num;

	return score;

}

double Stat_fuc::auc_frm_LR(vector<double> & indi_LR, vector< vector<int> > & compare_struct, double & var)
{
	int n_fam=compare_struct.size()/2;
	int i, j, k, idx_p, idx_h;
	double LR_p, LR_h, tmp_n_cmp=0, tmp_score=0, score=0;
	vector<double> fam_score;
	for(i=0; i<n_fam; i++){
		idx_p=i*2; idx_h=i*2+1;
		tmp_score=0; tmp_n_cmp=0;
		for(j=0; j<compare_struct[idx_p].size(); j++){
			LR_p=indi_LR[compare_struct[idx_p][j]];
			for(k=0; k<compare_struct[idx_h].size(); k++){
				LR_h=indi_LR[compare_struct[idx_h][k]];
				tmp_n_cmp+=1;
				if(LR_p>LR_h) tmp_score+=1;
				else if(LR_p==LR_h) tmp_score+=0.5;
			}
		}
		if(tmp_score!=0) tmp_score/=tmp_n_cmp;
		if(tmp_n_cmp!=0)fam_score.push_back(tmp_score);
		//fam_score.push_back(tmp_score);
	}

	for(i=0; i<fam_score.size(); i++) score+=fam_score[i];
	double fam_num=double(fam_score.size());
	score/=fam_num;

	var=0;
	for(i=0; i<fam_score.size(); i++){
		var+=fam_score[i]*fam_score[i];
	}
	var-=(fam_num*score*score);
	var/=(fam_num*(fam_num-1));

	return score;
}

void Stat_fuc::score_frm_LR(vector<double> & indi_LR, vector< vector<int> > & compare_struct, vector<double> & fam_score)
{
	int n_fam=compare_struct.size()/2;
	int i, j, k, idx_p, idx_h;
	double LR_p, LR_h, tmp_n_cmp=0, tmp_score=0;
	fam_score.clear();
	for(i=0; i<n_fam; i++){
		idx_p=i*2; idx_h=i*2+1;
		tmp_score=0; tmp_n_cmp=0;
		for(j=0; j<compare_struct[idx_p].size(); j++){
			LR_p=indi_LR[compare_struct[idx_p][j]];
			for(k=0; k<compare_struct[idx_h].size(); k++){
				LR_h=indi_LR[compare_struct[idx_h][k]];
				tmp_n_cmp+=1;
				if(LR_p>LR_h) tmp_score+=1;
				else if(LR_p==LR_h) tmp_score+=0.5;
			}
		}
		if(tmp_score!=0) tmp_score/=tmp_n_cmp;
		if(tmp_n_cmp!=0)fam_score.push_back(tmp_score);
		//fam_score.push_back(tmp_score);
	}
}


