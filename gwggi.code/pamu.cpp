#include "pamu.h"

ofstream LOG;
Pamu * PP;

int main(int argc, char* argv[]) 
{

	set_new_handler(gfun::NoMem);

	snpdt SDT;//store all the data
	
	Pamu P; //perform all the function
	P.initialize(&SDT);
	PP=&P;



	//////set initial value for globe parameters
	gfun::setInitialValue();

	// Command line arguments

	CArg a(argc,argv);

	gfun::getOutFileName(a);
	LOG.open(string(par::output_file_name + ".log").c_str());
	P.printLOG("\n"
		"@----------------------------------------------------------@\n"
		"|        GWGGI!                                            |\n"
		"|----------------------------------------------------------|\n"
		"|        Genome-wide Gene-Gene Interaction Analysis        |\n"
		"|----------------------------------------------------------|\n"
		"|        Copyright (C) 2011                                |\n"
		"|----------------------------------------------------------|\n"
		"|  Changshuai Wei, Michigan Stat Univ.                     |\n"
		"|        https://www.msu.edu/~changs18/index.html          |\n"
		"@----------------------------------------------------------@\n"
		"\n");

	gfun::setPar(a);

	// Time stamp

	P.printLOG("Writing log file at [ "+
		par::output_file_name + ".log ]\n");

	time_t curr=time(0);
	string tdstamp = (string)ctime(&curr);
	P.printLOG("Analysis started: " + tdstamp +"\n");

	////input data
	if (par::read_bitfile) P.readBinData();
	else if(par::read_ped) P.readPedData();

	///ifneeded split the data
	if( par::split) P.split_3();

	////forward selection algorithm
	////sample command 1:pamu --bfile GVT2D_HPFS --lmw --out HPFS
	////sample command 2:pamu --bfile GVT2D_HPFS --lmw --no-cv --out HPFS
	if(par::lmw_run) {
		if(par::cross_vali) P.cvLMW();
		else P.assocLMW();
	}

	///Appling model
	///sample command:pamu --bfile ND_90 --lmw-apl ND.lmw.lr --out NDpdt
	if(par::lmw_apl) P.aplLMW();

	//forward for family data
	if(par::flmw_run) {
		if(par::cross_vali) P.cvFLMW();
		else P.assocFLMW();
	}

	if(par::flmw_apl) P.aplFLMW();

	//tree assembling method
	////sample command:pamu --bfile ND_90 --tamw --converge --tree-depth 8 --out NDta
	if(par::tamw_run) P.scanTAMW();
	////sample command:pamu --bfile ND_90 --tamw-apl NDta.tamw.lr --out NDta
	if(par::tamw_apl) P.aplTAMW();

	//family tree assembling method
	if(par::ftamw_run) P.scanFTAMW();
	if(par::ftamw_apl) P.aplFTAMW();

	/////output data
	if (par::write_bitfile) P.write_BITFILE();
	else if(par::recode) P.write_PEDFILE();

	gfun::shutdown();
}


void Pamu::initialize(snpdt * data)
{
	_data=data;
}

void Pamu::clear()
{
	_data->clear();
}

void Pamu::printLOG(string s)
{
	LOG << s;
	LOG.flush();

	if (!par::silent)
	{
		cout << s;
		cout.flush();
	}
}

void Pamu::readBinData()
{
	_data->readBinData();
}

void Pamu::readPedData()
{
	_data->readGenMapPed();
}

void Pamu::write_BITFILE()
{
	_data->write_BITFILE();
}

void Pamu::write_PEDFILE()
{
	_data->writeGenMapPed();
}

snpdt* Pamu::dataPointer()
{
	snpdt * tmp=_data;
	return tmp;
}

void Pamu::split_3()
{
	int total=_data->totalIndi();
	int size=double(total)*0.6666;

	vector<int> ori;
	for(int i=0; i<total; i++){
		ori.push_back(i);
	}
	vector<int> rst, rst2;
	Stat_fuc::sampling_wr(ori,size, rst, par::seed, rst2 );

	_data->sampling(rst);
	
	string rstfile=par::output_file_name + ".train";
	if(par::read_bitfile){
		_data->write_BITFILE(rstfile);
	}else{
		_data->writeGenMapPed(rstfile);
	}

	_data->reset();
	_data->sampling(rst2);

	rstfile=par::output_file_name + ".test";
	if(par::read_bitfile){
		_data->write_BITFILE(rstfile);
	}else{
		_data->writeGenMapPed(rstfile);
	}

	_data->reset();
}

void Pamu::assocLMW()
{
	gfun::printLOG("\n*****LMW without Cross-validation*****\n");

	gfun::printLOG("Forward Scanning...\n");
	LMW asso;
	asso.initialize(_data);
	asso.iteration();
	asso.wtSubset(par::lmw_subset);
	asso.wtResult(par::lmw_scan_rst);

	kfold cv;
	cv.initialize(_data, &asso);
	cv.noCvWtLr(par::lmw_lr);
}

void Pamu::assocFLMW()
{

	gfun::printLOG("\n*****FLMW without Cross-validation*****\n");

	gfun::printLOG("Forward Scanning...\n");
	FLMW asso;
	asso.initialize(_data);
	asso.iteration();
	asso.wtSubset(par::lmw_subset);
	asso.wtResult(par::lmw_scan_rst);

	fkfold cv;
	cv.initialize(_data, &asso);
	cv.noCvWtLr(par::lmw_lr);
}

void Pamu::cvLMW()
{
	gfun::printLOG("\n*****LMW with Cross-validation*****\n");

	gfun::printLOG("Forward Scanning...\n");

	LMW asso;
	asso.initialize(_data);
	asso.iteration();
	asso.wtSubset(par::lmw_subset);
	asso.wtResult(par::lmw_scan_rst);

	gfun::printLOG("Cross validating...\n");
	kfold cv;
	cv.initialize(_data, &asso);
	cv.crossVali(par::lmw_lr);
	cv.wtResult(par::lmw_cv_rst);
}

void Pamu::cvFLMW()
{

	gfun::printLOG("\n*****FLMW with Cross-validation*****\n");

	gfun::printLOG("Forward Scanning...\n");
	FLMW asso;
	asso.initialize(_data);
	asso.iteration();
	asso.wtSubset(par::lmw_subset);
	asso.wtResult(par::lmw_scan_rst);

	gfun::printLOG("Cross validating...\n");
	fkfold cv;
	cv.initialize(_data, &asso);
	cv.crossVali(par::lmw_lr);
	cv.wtResult(par::lmw_cv_rst);
}

void Pamu::aplLMW()
{
	gfun::printLOG("\n*****Apply LMW*****\n");

	gfun::printLOG("Applying model...\n");
	kfold apl;
	apl.applyModel(_data,par::lmw_lr);
	apl.wtResult(par::lmw_apl_rst);
	apl.wtSubset(par::lmw_subset);
}

void Pamu::aplFLMW()
{
	gfun::printLOG("\n*****Apply FLMW*****\n");

	gfun::printLOG("Applying model...\n");
	fkfold apl;
	apl.applyModel(_data,par::lmw_lr);
	apl.wtResult(par::lmw_apl_rst);
	apl.wtSubset(par::lmw_subset);
}

void Pamu::scanTAMW()
{
	gfun::printLOG("\n*****TAMW*****\n");

	gfun::printLOG("Building assembling model...\n");
	TAMW ta;
	ta.initialize(_data);
	ta.assembling();
	ta.wtResult(par::tamw_rst);
	ta.wtPredictorStat(par::tamw_smr);
	if(par::moniter) ta.wtConverge(par::moniter_f);
}

void Pamu::aplTAMW()
{
	gfun::printLOG("\n*****Apply TAMW*****\n");

	gfun::printLOG("Applying assembling model...\n");
	TAMW ta;
	ta.initialize(_data);
	ta.apply_model(par::tamw_lr);
	ta.wtResult(par::tamw_rst);
	ta.wtPredictorStat(par::tamw_smr);
	if(par::moniter) ta.wtConverge(par::moniter_f);
}

void Pamu::scanFTAMW()
{
	gfun::printLOG("\n*****FTAMW*****\n");

	gfun::printLOG("Building assembling model...\n");

	FTAMW ta;
	ta.initialize(_data);
	ta.assembling();
	ta.wtResult(par::tamw_rst);
	ta.wtPredictorStat(par::tamw_smr);
	if(par::moniter) ta.wtConverge(par::moniter_f);
}

void Pamu::aplFTAMW()
{
	gfun::printLOG("\n*****Apply FTAMW*****\n");

	gfun::printLOG("Applying assembling model...\n");
	FTAMW ta;
	ta.initialize(_data);
	ta.apply_model(par::tamw_lr);
	ta.wtResult(par::tamw_rst);
	ta.wtPredictorStat(par::tamw_smr);
	if(par::moniter) ta.wtConverge(par::moniter_f);
}









