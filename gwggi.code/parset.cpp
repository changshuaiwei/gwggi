#include "parset.h"

/*****globle function*******/

void gfun::NoMem()
{
	cerr << "*****************************************************\n"
		<< "* FATAL ERROR    Exhausted system memory            *\n"
		<< "*                                                   *\n"
		<< "* You need a smaller dataset or a bigger computer...*\n"
		<< "*                                                   *\n"
		<< "* Forced exit now...                                *\n"
		<< "*****************************************************\n\n";
	exit(1);
}

void gfun::shutdown()
{
	time_t curr=time(0);
	string tdstamp = ctime(&curr);

	if (!par::silent) cout << "\nAnalysis finished: " + tdstamp +"\n";
	LOG << "\nAnalysis finished: " + tdstamp +"\n";

	LOG.close();

	PP->clear();


	system("PAUSE");
	exit(0);
}

void gfun::error(string msg)
{
	cerr << "\nERROR: " << msg << "\n";
	LOG << "\nERROR: " << msg << "\n";  
	LOG.close();

	PP->clear();

	exit(1);
}

void gfun::getOutFileName(CArg &a)
{
	if (a.find("--out")) 
	{
		par::output_file_name = a.value("--out");
	}

	if (a.find("--silent")) 
	{
		par::silent = true;
	}

}

void gfun::setInitialValue()
{

}

void gfun::setPar(CArg &a)
{
	gfun::printLOG("Setting Parameters for command: "); a.showArg();

	if (a.find("--bfile")) 
	{ 
		par::read_bitfile = true;
		par::fileroot = a.value("--bfile");
		par::bitfilename = par::fileroot + ".bed";
		par::bitfilename_map= par::fileroot + ".bim";
		par::famfile = par::fileroot + ".fam";
	}

	if (a.find("--file")) 
	{ 
		par::read_ped = true;
		par::fileroot = a.value("--file");
		par::pedfile = par::fileroot + ".ped";
		par::mapfile = par::fileroot + ".map";
		par::genfile = par::fileroot + ".gen";
	}

	if(a.find("--split"))
	{
		par::split =true;
	}

	if(a.find("--lmw"))
	{
		par::lmw_run=true;
		par::lmw_scan_rst=par::output_file_name + ".scan";
		par::lmw_cv_rst=par::output_file_name + ".cv";
		par::lmw_subset=par::output_file_name + ".subset";
		par::lmw_lr=par::output_file_name + ".lmw.lr";
	}

	if(a.find("--lmw-apl")){
		par::lmw_apl=true;
		par::lmw_lr=a.value("--lmw-apl");
		par::lmw_apl_rst=par::output_file_name + ".apl";
		par::lmw_subset=par::output_file_name + ".subset";
	}

	if(a.find("--flmw"))
	{
		par::flmw_run=true;
		par::lmw_scan_rst=par::output_file_name + ".scan";
		par::lmw_cv_rst=par::output_file_name + ".cv";
		par::lmw_subset=par::output_file_name + ".subset";
		par::lmw_lr=par::output_file_name + ".flmw.lr";
	}

	if(a.find("--flmw-apl")){
		par::flmw_apl=true;
		par::lmw_lr=a.value("--flmw-apl");
		par::lmw_apl_rst=par::output_file_name + ".apl";
		par::lmw_subset=par::output_file_name + ".subset";
	}

	if(a.find("--tamw")){
		par::tamw_run=true;
		par::tamw_lr=par::output_file_name + ".tamw.lr";
		par::tamw_rst=par::output_file_name + ".tamw.rst";
		par::tamw_smr=par::output_file_name + ".tamw.smr";
		par::moniter_f=par::output_file_name + ".tamw.monitor";
	}

	if(a.find("--tamw-apl")){
		par::tamw_apl=true;
		par::tamw_lr=a.value("--tamw-apl");
		par::tamw_rst=par::output_file_name + ".tamw.rst";
		par::tamw_smr=par::output_file_name + ".tamw.smr";
		par::moniter_f=par::output_file_name + ".tamw.monitor";
	}

	if(a.find("--ftamw")){
		par::ftamw_run=true;
		par::tamw_lr=par::output_file_name + ".ftamw.lr";
		par::tamw_rst=par::output_file_name + ".ftamw.rst";
		par::tamw_smr=par::output_file_name + ".ftamw.smr";
		par::moniter_f=par::output_file_name + ".ftamw.monitor";
	}

	if(a.find("--ftamw-apl")){
		par::tamw_apl=true;
		par::tamw_lr=a.value("--ftamw-apl");
		par::tamw_rst=par::output_file_name + ".ftamw.rst";
		par::tamw_smr=par::output_file_name + ".ftamw.smr";
		par::moniter_f=par::output_file_name + ".ftamw.monitor";
	}

	if(a.find("--vote"))  par::vote=true;

	if(a.find("--converge")) par::moniter=true;

	if(a.find("--match-names")){
		par::match_name=true;
	}

	if(a.find("--maxauc")) par::largest_auc=str2double(a.value("--maxauc"));

	if(a.find("--no-cv")) par::cross_vali=false; else par::cross_vali=true;

	if(a.find("--lr-sub"))
	{
		par::out_nom_LR=true;
		par::out_nom_LR_f=par::output_file_name + ".nom.lr";
	}

	if(a.find("--seed")) par::seed=str2int(a.value("--seed"));

	if(a.find("--clsf")) {
		par::n_classifier=str2int(a.value("--clsf"));
		par::clf_ln=false; par::clf_sqrt=false;
	}

	if(a.find("--clsf-sqrt")){
		par::clf_sqrt=true; par::clf_ln=false;
	}

	if(a.find("--tree-depth")){
		par::tree_depth=str2int(a.value("--tree-depth"));
	}

	if(a.find("--ntree")){
		par::max_ntree=str2int(a.value("--ntree"));
	}

	if(a.find("--auto-depth")){
		par::auto_tree_depth=true;
	}

	if (a.find("--make-bed")) 
	{
		par::write_bitfile = true;
	}

	if (a.find("--recode")) 
	{
		par::recode = true;
	}

	if(a.find("--showcv"))
	{
		par::show_crossvali=true;
		par::show_iteration=true;
	}

	if(a.find("--showscan")){
		par::show_iteration=true;
	}

	if(a.find("--showntree")){
		par::show_ntree=true;
	}

	if(a.find("--maxcvauc")){
		par::choose_first_peak=false;
	}

	if(a.find("--hiorder")){
		par::most_nsnp=str2int(a.value("--hiorder"));
	}

	if(a.find("--ckmissing")){
		par::ckmissing=true;
		par::missing_rate=str2double(a.value("--ckmissing"));
	}

	if(a.find("--td-burnin")){
		par::burnin_tree_depth=true;
	}

	if(a.find("--nburnin")){
		par::ntree_burnin=str2int(a.value("--nburnin"));
	}

	if(a.find("--hz")){
		par::exclude_1=false;
	}

}

void gfun::checkFileExists(string f)
{

	ifstream inp;

	inp.open(f.c_str(), ifstream::in);
	if(inp.fail())
	{
		inp.clear(ios::failbit);
		inp.close();
		string msg = "No file [ " + f + " ] exists.";
		gfun::error(msg);
	}
	inp.close();
	return;

}

void gfun::printLOG(string s)
{
	LOG << s;
	LOG.flush();

	if (!par::silent)
	{
		cout << s;
		cout.flush();
	}
}

int gfun::split_string(const string &str, vector <string> &vec_str, string separator)
{
	if(str.empty()) return 0;
	vec_str.clear();

	int i=0;
	bool look=false;
	string str_buf;
	string symbol_pool="`1234567890-=~!@#$%^&*()_+qwertyuiop[]\\asdfghjkl;'zxcvbnm,./QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>? \t\n";
	string::size_type pos;

	for(i=0; i<separator.size(); i++) 
	{
		pos=symbol_pool.find(separator[i]);
		if( pos!=string::npos ) symbol_pool.erase(symbol_pool.begin()+pos);
	}

	for(i=0; i<str.size(); i++)
	{		
		if( symbol_pool.find(str[i])!=string::npos )
		{
			if(!look) look=true;
			str_buf += str[i];
		}
		else
		{
			if(look) 
			{
				look=false;
				vec_str.push_back(str_buf);
				str_buf.erase(str_buf.begin(), str_buf.end());
			}
		}
	}
	if(look) vec_str.push_back(str_buf);

	return vec_str.size();
}

/**********std**************/

int std::str2int (const string &str) {
	stringstream ss(str);
	int n;
	ss >> n;
	return n;
}

string std::int2str (int n) {
	stringstream ss;
	ss << n;
	return ss.str();
}

double std::str2double(const string &str)
{
	stringstream ss(str);
	double n;
	ss >> n;
	return n;
}



/*****CArg class*******/

CArg::CArg(int n, char *argv[])
{
	_n = n;
	a.push_back(argv[0]);
	_parsed.resize(_n,false);
	_option.resize(_n,false);
	for (int i=1 ; i < _n ; i++ )
		a.push_back(argv[i]);
	_original = a;
}

void CArg::showArg()
{
	for(int i=0; i<a.size(); i++){
		gfun::printLOG(a[i]+" ");
	}
	gfun::printLOG("\n");
}

bool CArg::find(string s)
{  
	for (int i=0 ; i < _n ; i++ )
		if (a[i] == s) { _parsed[i]=true; return true; }
		return false;
}

string CArg::value(string s)
{
	for (int i=0 ; i < _n ; i++ )
		if (a[i] == s && (i+1 < _n) ) { _parsed[i+1]=true; return a[i+1]; }
		gfun::error("Missing an argument for "+s);
} 

vector<string> CArg::value(string s, int c)
{
	vector<string> r(0);

	for (int i=0 ; i < _n ; i++ )
		if (a[i] == s && (i+1 < _n) ) 
		{
			for (int j=1;j<=c;j++) 
			{
				if ( (i+j) < a.size() )
				{
					_parsed[i+j]=true;
					r.push_back(a[i+j]);
				}
				else
					gfun::error("Not enough arguments given for option: "+s+" ");
			} 
		}

		if (r.size() != c) gfun::error("Not enough arguments given for option: "+s+" ");
		return r;
} 



///////////////////////

string par::output_file_name="pamuRst";
bool par::silent=false;
double par::replace0_cor=0.5;
double par::replace0_sml=0.00001;
bool par::read_bitfile=false;
string par::fileroot="";
string par::bitfilename="";
string par::bitfilename_map="";
string par::famfile="";
bool par::qt=false;//quantitative trait
bool par::bt=true;
bool par::coding01=false;
int par::missing_int=-9;
string par::missing_str=".";
string par::out_missing_phenotype = "-9";
bool par::SNP_major=true;
bool par::write_bitfile=false;
bool par::out_SNP_major = true;
string par::recode_delimit = "\t";
string par::recode_indelimit = "";
string par::out_missing_genotype = "0";

bool par::read_ped = false;
string par::pedfile = "";
string par::mapfile = "";
string par::genfile = "";
bool par:: recode=false;

bool par::split=false;

int par::most_nsnp=10;
double par::largest_auc=0.9;
bool par::out_nom_LR=false;
string par::out_nom_LR_f="pamuRst.nom.lr";
bool par::exclude_1=true;

bool par::lmw_run=false;
bool par::flmw_run=false;
string par::lmw_scan_rst="pamuRst.scan";
string par::lmw_cv_rst="pamuRst.cv";
int par::n_fold=10;
int par::seed=-2011;
bool par::choose_first_peak=true;
bool par::cross_vali=true;
string par::lmw_subset="pamuRst.subset";
string par::lmw_lr="pamuRst.lmw.lr";
bool par::match_name=false;
bool par::lmw_apl=false;
bool par::flmw_apl=false;
string par::lmw_apl_rst="pamuRst.apl";
int par::tree_depth=6;
int par::n_classifier=2;
bool par::clf_ln=true;
bool par::clf_sqrt=false;
string par::tamw_lr="pamuRst.tamw.lr";
bool par::auto_tree_depth=false;
int par::between_ntree=10;
int par::max_ntree=2000;
double par::thrh_Zscore=0;
int par::thrh_seltimes=5;
bool par::tamw_run=false;
string par::tamw_rst="pamuRst.tamw.rst";
string par::tamw_smr="pamuRst.tamw.smr";
bool par::moniter=false;
string par::moniter_f="pamuRst.tamw.monitor";
bool par::tamw_apl=false;
bool par::ftamw_run=false;
bool par::ftamw_apl=false;
bool par::show_iteration=false;
bool par::show_crossvali=false;
bool par::show_ntree=false;
bool par::ckmissing=false;
double par::missing_rate=0.1;
bool par::burnin_tree_depth=false;
int par::ntree_burnin=20;
bool par::vote=false;

