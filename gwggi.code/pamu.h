#ifndef _PAMU_H
#define _PAMU_H


#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <list>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <new>

#include "parset.h"
#include "snpdt.h"
#include "lmw.h"
#include "kfold.h"
#include "tamw.h"


class Individual;
class Locus;
class CSNP;
class snpdt;

using namespace std;

class Pamu
{
public:
	Pamu(){
		_data=0;
	}
	~Pamu(){
		_data=0;
	}
	
	void initialize(snpdt *);
	void clear();
	void printLOG(string);
	void readBinData();
	void readPedData();
	void write_BITFILE();
	void write_PEDFILE();
	void assocLMW();
	void cvLMW();
	void aplLMW();
	void assocFLMW();
	void cvFLMW();
	void aplFLMW();
	void scanTAMW();
	void aplTAMW();
	void scanFTAMW();
	void aplFTAMW();
	snpdt* dataPointer();
	void split_3();

private:
	snpdt * _data;

};


#endif