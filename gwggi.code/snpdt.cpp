#include "snpdt.h"

snpdt::snpdt()
{
	_sample.resize(0);
	_locus.resize(0);
	_phenotype.resize(0);
	_SNP.resize(0);
	_resample=false;
	_relocus=false;
}

snpdt::~snpdt()
{
	if(_resample){
		for (int l=0; l<_SNP.size(); l++)
			delete _SNP[l];
		_sample=_sample_bkup; _sample_bkup.clear();
		_locus=_locus_bkup; _locus_bkup.clear();
		_SNP=_SNP_bkup;	_SNP_bkup.clear();

	}
	for (int i=0; i<_sample.size(); i++)
		delete _sample[i];

	for (int l=0; l<_locus.size(); l++)
		delete _locus[l];      

	for (int l=0; l<_SNP.size(); l++)
		delete _SNP[l];
}

void snpdt::test()
{
	cout<<"\ntest snpdt class\n";
	//add code here...
}

void snpdt::readBinData()
{
	// We cannot assume the file will be in order
	// Check files exist

	gfun::checkFileExists(par::famfile);
	gfun::checkFileExists(par::bitfilename_map);
	gfun::checkFileExists(par::bitfilename);

	gfun::printLOG("Reading map (extended format) from [ " 
		+ par::bitfilename_map + " ] \n");

	vector<Locus> ordered;

	ifstream MAP(par::bitfilename_map.c_str(), ios::in);
	MAP.clear();

	int c=0;
	while(!MAP.eof())
	{

		Locus * loc = new Locus;

		MAP >> loc->chr   // will automatically by numeric
			>> loc->name 
			>> loc->pos   
			>> loc->bp     
			>> loc->allele1
			>> loc->allele2;

		if ( MAP.eof() )
		{
			delete loc; 
			continue;
		}
		else if ( MAP.fail() )
		{
			delete loc;
			gfun::error("Problem reading BIM file, line " + int2str(c+1) + "\n");
		}

		// Use the frequency slot temporarily to 
		// store order information
		loc->freq = c++;


		// Always included, but not always in correct order
		if (loc->name!="") 
		{
			_locus.push_back(loc);
			ordered.push_back(*loc);
		}  
		else
			delete loc;
	}


	gfun::printLOG( int2str(_locus.size()) 
		+ " markers to be included from [ " +
		par::bitfilename_map + " ]\n");

	MAP.clear();
	MAP.close();

	if ( _locus.size() == 0 ) 
		gfun::shutdown();

	///////////////////////////////////////////////
	// Build ordered table, so that genotypes can 
	// be inserted in correct order; then swap locus 
	// file over

	// Sort vector of pointers to Locus
	stable_sort(_locus.begin(),_locus.end(),less<Locus*>());

	// Sort normal vector Locus
	stable_sort(ordered.begin(),ordered.end());

	c=0;
	for (int i=0; i<_locus.size(); i++)
	{
		// swap file order into locus position
		// and make all same chromosome
		ordered[i].bp = (int)ordered[i].freq;
		ordered[i].chr = 1;
		// keep track of genetic order, but only 
		ordered[i].freq = c++;
	}

	// resort to get lookup table
	stable_sort(ordered.begin(),ordered.end());  


	///////////////////////////////////////
	// Do we want to look at all the data?

	vector<int> include(0);
	int nl_actual = _locus.size();

	// If we do want to look at all the data
	for (int j=0; j<ordered.size(); j++)
		include.push_back((int)ordered[j].freq);



	//////////////////////////////
	// Read individual information

	readFamFile(par::famfile);


	gfun::printLOG("Reading genotype bitfile from [ " 
		+ par::bitfilename + " ] \n");

	ifstream BIT;
	bool bfile_SNP_major = openBinaryFile(par::bitfilename, BIT);

	if (!bfile_SNP_major) gfun::error("Sorry, but we don't support individual-major mode now\n");


	/////////////////////////////////////////
	// Read entire genotype data into a temp
	// array


	//////////////////////////////
	// Allocate space for SNPs


	for (int i=0; i<nl_actual; i++)
	{
		CSNP * newlocus = new CSNP;
		newlocus->one.resize( _sample.size() );
		newlocus->two.resize( _sample.size() );
		_SNP.push_back(newlocus);
	}



	///////////////////////////
	// SNP-major mode

	if (bfile_SNP_major)
	{

		CSNP * snp;

		// Outer loop for SNPs
		int s=0;
		while (s<_locus.size()) // for all SNPs
		{

			// Do we want to include this SNP?
			if ( include[s] > -1 )
				snp = _SNP[ include[s] ];	  
			else
				snp = NULL;

			// Inner loop for individuals

			int indx = 0;
			int ss = _sample.size();

			while ( indx < ss )
			{

				bitset<8> b;
				
				char ch[1];
				BIT.read(ch,1);
				if (!BIT) 
					gfun::error("Problem with the BED file...has the FAM/BIM file been changed?\n");		
				b = ch[0];

				int c=0;

				while (c<7 && indx < ss ) 
				{
					if (snp)
					{
						snp->one[indx] = b[c++];
						snp->two[indx] = b[c++];

					}
					else
					{
						c+=2;
					}	  
					++indx;
				}

			}	  

			s++;
		}

		// Set file mode
		par::SNP_major = true;

	}



	// Check that we got what we expected

	char ch[1];
	BIT.read(ch,1);
	if (BIT) 
		gfun::error("Problem with the BED file... has the FAM/BIM file been changed?\n");

	BIT.clear();
	BIT.close();


	////////////////////////////////////////
	// If need be, now prune the MAP file 
	// i.e. if --chr or --from/--to were used

}




bool snpdt::openBinaryFile(string s, ifstream & BIT)
{

	BIT.open(s.c_str(), ios::in | ios::binary);

	// 1) Check for magic number
	// 2) else check for 0.99 SNP/Ind coding
	// 3) else print warning that file is too old

	char ch[1];
	BIT.read(ch,1);
	bitset<8> b;
	b = ch[0];	  

	bool bfile_SNP_major = false;
	bool v1_bfile = true;

	// If v1.00 file format
	// Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
	if (   ( b[2] && b[3] && b[5] && b[6] ) && 
		! ( b[0] || b[1] || b[4] || b[7] )    )
	{

		// Next number
		BIT.read(ch,1);
		b = ch[0];	  
		if (   ( b[0] && b[1] && b[3] && b[4] ) && 
			! ( b[2] || b[5] || b[6] || b[7] )    )
		{
			// Read SNP/Ind major coding
			BIT.read(ch,1);
			b = ch[0];	  
			if ( b[0] ) bfile_SNP_major = true;
			else bfile_SNP_major = false;

			if (bfile_SNP_major) 
				gfun::printLOG("Detected that binary PED file is v1.00 SNP-major mode\n");
			else
				gfun::printLOG("Detected that binary PED file is v1.00 individual-major mode\n");

		} else v1_bfile = false;

	} else v1_bfile = false;


	// Reset file if < v1
	if ( ! v1_bfile ) 
	{
		gfun::error("Error, old BED file <v1.00 \n");
		  
	}

	return bfile_SNP_major;

}

void snpdt::readFamFile(string filename)
{

	//////////////////////////////
	// Read pedigree information

	// Initially, assume a binary trait
	par::qt = false;
	par::bt = true; 

	gfun::printLOG("Reading pedigree information from [ " 
		+ filename + " ] \n");

	gfun::checkFileExists( filename );

	ifstream PED;
	PED.open(filename.c_str());
	PED.clear();

	vector<Individual*> ambiguous;

	int nmale = 0;
	int nfemale = 0;
	int nambig = 0;

	int c=0;
	while(!PED.eof())
	{

		Individual * person = new Individual;

		// First 6 obligatory fields
		string phenotype;

		PED >> person->fid 
			>> person->iid 
			>> person->pat
			>> person->mat
			>> person->sexcode
			>> person->pheno_str;

		phenotype=person->pheno_str;

		// Are we using 0/1 coding?
		if (par::coding01) 
		{
			if ( phenotype == "1" ) 
				phenotype = "2";      
			else if ( phenotype == "0" )
				phenotype = "1";
			else 
				phenotype = "0";
		}

		///set the binary phenotype
		if(phenotype == "2") person->aff = 1;
		else if (phenotype == "1") person->aff = 0;

		//set the family number
		if(person->fid == "-9") person->nfid= -9;
		else istringstream(person->fid) >> person->nfid;

		// Skip last empty line that gets read
		if (person->fid=="") 
		{
			delete person;
			break;
		}

		// Check sex
		if (person->sexcode=="1")
		{
			person->sex = 1; // male
			nmale++;
		}
		else if (person->sexcode=="2")
		{
			person->sex = 0;  // female (default)
			nfemale++;
		}
		else 
		{
			person->sex = par::missing_int;//missing
			nambig++;
		}

		// Increase person counter
		c++;

		// Add individual to list
		_sample.push_back(person);

	}

	PED.clear();
	PED.close();





	gfun::printLOG(int2str(c)+" individuals read from [ " 
		+ filename + " ] \n");
	
	if (par::bt) 
	{

		if (par::coding01) 
			gfun::printLOG("Assuming a disease phenotype (0=unaff, 1=aff, other=miss)\n");
		else
		{
			gfun::printLOG("Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)\n");
		}
	}
	
}


void snpdt::clear()
{
	if(_resample){
		for (int l=0; l<_SNP.size(); l++)
			delete _SNP[l];
		_sample=_sample_bkup; _sample_bkup.clear();
		_locus=_locus_bkup; _locus_bkup.clear();
		_SNP=_SNP_bkup;	_SNP_bkup.clear();

	}
	for (int i=0; i<_sample.size(); i++)
		delete _sample[i];

	for (int l=0; l<_locus.size(); l++)
		delete _locus[l];      

	for (int l=0; l<_SNP.size(); l++)
		delete _SNP[l];
}

void snpdt::subset(vector<int> sample, vector<int> locus)
{
	reset();

	_sample_bkup=_sample; _sample.clear();
	_SNP_bkup=_SNP; _SNP.clear();
	_locus_bkup=_locus; _locus.clear();
	_resample=true;

	for(int i=0; i<sample.size(); i++){
		_sample.push_back(_sample_bkup[sample[i]]);
	}

	for(int i=0; i<locus.size(); i++){
		_locus.push_back(_locus_bkup[locus[i]]);
	}

	for (int i=0; i<_locus.size(); i++)
	{
		CSNP * newlocus = new CSNP;
		newlocus->one.resize( _sample.size() );
		newlocus->two.resize( _sample.size() );
		_SNP.push_back(newlocus);
	}

	for(int i=0; i<locus.size(); i++){
		for(int j=0; j<sample.size(); j++){
			_SNP[i]->one[j]=_SNP_bkup[locus[i]]->one[sample[j]];
			_SNP[i]->two[j]=_SNP_bkup[locus[i]]->two[sample[j]];
		}
	}

}

void snpdt::sampling(vector<int> sample)
{
	reset();

	_sample_bkup=_sample; _sample.clear();
	_SNP_bkup=_SNP; _SNP.clear();
	_locus_bkup=_locus;
	_resample=true; _relocus=true;

	for(int i=0; i<sample.size(); i++){
		_sample.push_back(_sample_bkup[sample[i]]);
	}

	for (int i=0; i<_locus.size(); i++)
	{
		CSNP * newlocus = new CSNP;
		newlocus->one.resize( _sample.size() );
		newlocus->two.resize( _sample.size() );
		_SNP.push_back(newlocus);
	}

	for(int i=0; i<_locus.size(); i++){
		for(int j=0; j<sample.size(); j++){
			_SNP[i]->one[j]=_SNP_bkup[i]->one[sample[j]];
			_SNP[i]->two[j]=_SNP_bkup[i]->two[sample[j]];
		}
	}
}

void snpdt::subset(vector<int> locus)
{
	reset();
	_sample_bkup=_sample;
	_SNP_bkup=_SNP; _SNP.clear();
	_locus_bkup=_locus; _locus.clear();
	_relocus=true;

	for(int i=0; i<locus.size(); i++){
		_locus.push_back(_locus_bkup[locus[i]]);
		_SNP.push_back(_SNP_bkup[locus[i]]);
	}
}

void snpdt::reset()
{
	if(_resample){
		_relocus=false; _resample=false;

		for (int l=0; l<_SNP.size(); l++)
			delete _SNP[l];

		_SNP=_SNP_bkup; _SNP_bkup.clear();
		_sample=_sample_bkup; _sample_bkup.clear();
		_locus=_locus_bkup; _locus_bkup.clear();

	}else if(_relocus){
		_relocus=false;
		_SNP=_SNP_bkup; _SNP_bkup.clear();
		_sample=_sample_bkup; _sample_bkup.clear();
		_locus=_locus_bkup; _locus_bkup.clear();
	}

}

void snpdt::write_BITFILE()
{

	gfun::printLOG( "Writing pedigree information to [ " 
		+ par::output_file_name + ".fam ] \n");
	ofstream BIT((par::output_file_name+".fam").c_str(), ios::out);
	// For each individual

	for (int i=0;i<_sample.size();i++)
	{
		Individual * person = _sample[i];
		BIT << person->fid << " "
			<< person->iid << " "
			<< person->pat << " "
			<< person->mat << " "
			<< person->sexcode << " "
			<< person->pheno_str << "\n";
	}
	BIT.clear();
	BIT.close();

	gfun::printLOG( "Writing map (extended format) information to [ " 
		+ par::output_file_name + ".bim ] \n");

	BIT.open((par::output_file_name+".bim").c_str(), ios::out);

	for (int l=0;l<_locus.size();l++)
	{

		if (_locus[l]->allele1=="")
		{
			if (_locus[l]->allele2!="0") _locus[l]->allele1="0";
			else _locus[l]->allele1="X";
		}

		if (_locus[l]->allele2=="")
			_locus[l]->allele2="0";

		BIT << _locus[l]->chr << "\t"
			<< _locus[l]->name << "\t"
			<< _locus[l]->pos << "\t"
			<< _locus[l]->bp << "\t" 
			<< _locus[l]->allele1 << "\t"
			<< _locus[l]->allele2 << "\n";
	}

	BIT.clear();
	BIT.close();


	//////////////////////////////////////
	// Save genotype data in BITFILE format

	gfun::printLOG("Writing genotype bitfile to [ " 
		+ par::output_file_name + ".bed ] \n");

	BIT.open((par::output_file_name+".bed").c_str(), ios::out | ios::binary);


	gfun::printLOG("Using (default) SNP-major mode\n");


	bitset<8> b;
	char ch[1];


	// Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file

	b.reset();
	b.set(2);  b.set(3);  b.set(5);  b.set(6);
	ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);

	b.reset();
	b.set(0);  b.set(1);  b.set(3);  b.set(4);
	ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);


	// BIT represents status of SNP-major (true) or Ind-major (false)

	b.reset();
	b.set(0);
	ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);


	//////////////////////////////////////////
	// Now consider genotypes: SNP-major mode

	if (par::out_SNP_major)
	{  

		vector<CSNP*>::iterator s = _SNP.begin();

		// Outer loop over SNPs
		while ( s != _SNP.end() )
		{

			vector<bool>::iterator i1 = (*s)->one.begin();
			vector<bool>::iterator i2 = (*s)->two.begin();
			//	  vector<Individual*>::iterator person = sample.begin();

			// Inner loop over individuals
			while ( i1 != (*s)->one.end() )
			{
				bitset<8> b;
				b.reset();
				int c=0;      

				while (c<8 && i1 != (*s)->one.end() )
				{
					if ( *(i1) ) b.set(c);
					i1++;
					c++;

					if ( *(i2) ) b.set(c);		  
					i2++;
					c++;

					//person++;
				}

				char ch[1];
				ch[0] = (char)b.to_ulong();
				BIT.write(ch,1);

			}

			// next SNP
			s++;
		}
	}


	BIT.close();

}


void snpdt::write_BITFILE(string file_root)
{

	gfun::printLOG( "Writing pedigree information to [ " 
		+ file_root + ".fam ] \n");
	ofstream BIT((file_root+".fam").c_str(), ios::out);
	// For each individual

	for (int i=0;i<_sample.size();i++)
	{
		Individual * person = _sample[i];
		BIT << person->fid << " "
			<< person->iid << " "
			<< person->pat << " "
			<< person->mat << " "
			<< person->sexcode << " "
			<< person->pheno_str << "\n";
	}
	BIT.clear();
	BIT.close();

	gfun::printLOG( "Writing map (extended format) information to [ " 
		+ file_root + ".bim ] \n");

	BIT.open((file_root+".bim").c_str(), ios::out);

	for (int l=0;l<_locus.size();l++)
	{

		if (_locus[l]->allele1=="")
		{
			if (_locus[l]->allele2!="0") _locus[l]->allele1="0";
			else _locus[l]->allele1="X";
		}

		if (_locus[l]->allele2=="")
			_locus[l]->allele2="0";

		BIT << _locus[l]->chr << "\t"
			<< _locus[l]->name << "\t"
			<< _locus[l]->pos << "\t"
			<< _locus[l]->bp << "\t" 
			<< _locus[l]->allele1 << "\t"
			<< _locus[l]->allele2 << "\n";
	}

	BIT.clear();
	BIT.close();


	//////////////////////////////////////
	// Save genotype data in BITFILE format

	gfun::printLOG("Writing genotype bitfile to [ " 
		+ file_root + ".bed ] \n");

	BIT.open((file_root+".bed").c_str(), ios::out | ios::binary);


	gfun::printLOG("Using (default) SNP-major mode\n");


	bitset<8> b;
	char ch[1];


	// Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file

	b.reset();
	b.set(2);  b.set(3);  b.set(5);  b.set(6);
	ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);

	b.reset();
	b.set(0);  b.set(1);  b.set(3);  b.set(4);
	ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);


	// BIT represents status of SNP-major (true) or Ind-major (false)

	b.reset();
	b.set(0);
	ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);


	//////////////////////////////////////////
	// Now consider genotypes: SNP-major mode

	if (par::out_SNP_major)
	{  

		vector<CSNP*>::iterator s = _SNP.begin();

		// Outer loop over SNPs
		while ( s != _SNP.end() )
		{

			vector<bool>::iterator i1 = (*s)->one.begin();
			vector<bool>::iterator i2 = (*s)->two.begin();
			//	  vector<Individual*>::iterator person = sample.begin();

			// Inner loop over individuals
			while ( i1 != (*s)->one.end() )
			{
				bitset<8> b;
				b.reset();
				int c=0;      

				while (c<8 && i1 != (*s)->one.end() )
				{
					if ( *(i1) ) b.set(c);
					i1++;
					c++;

					if ( *(i2) ) b.set(c);		  
					i2++;
					c++;

					//person++;
				}

				char ch[1];
				ch[0] = (char)b.to_ulong();
				BIT.write(ch,1);

			}

			// next SNP
			s++;
		}
	}


	BIT.close();

}



string snpdt::genotypeToStr(int i, int l)
{

	// Return a genotype in suitable text format to be 
	// written to most output files

	string a1 = _locus[l]->allele1;
	string a2 = _locus[l]->allele2;

	bool s1 = _SNP[l]->one[i];
	bool s2 = _SNP[l]->two[i];

	if      ( (!s1) && (!s2) )  
		return a1+par::recode_indelimit+a1;
	else if ( (!s1) && s2 ) 
		return a1+par::recode_indelimit+a2;
	else if (  s1   && s2 ) 
		return a2+par::recode_indelimit+a2;
	else 
		return par::out_missing_genotype+ par::recode_indelimit+par::out_missing_genotype;

}

string snpdt::intToGenotype(int snp, int code)
{
	string a1 = _locus[snp]->allele1;
	string a2 = _locus[snp]->allele2;

	if(code==0) return a1+par::recode_indelimit+a1;
	else if(code==1) return a1+par::recode_indelimit+a2;
	else if(code==2) return a2+par::recode_indelimit+a2;
	else if(code==-9) return par::out_missing_genotype + par::recode_indelimit+par::out_missing_genotype;
}

int snpdt::genoStrToInt(int snp, string A1, string A2)
{
	string a1 = _locus[snp]->allele1;
	string a2 = _locus[snp]->allele2;

	if(A1==A2){
		if(A1==a1) return 0;
		else if(A1==a2) return 2;
		else return par::missing_int;
	}else{
		if( (A1==a1 && A2==a2) || (A1==a2 && A2==a1) ) return 1;
		else return par::missing_int;
	}
}



void snpdt::writeIntToGenotype(int indi, int snp, int code)
{
	bool s1;
	bool s2;

	if(code==0) { s1=false; s2=false; }
	else if(code==1) { s1=false; s2=true; }
	else if(code==2) { s1=true; s2=true; }
	else if(code==-9) { s1=true; s2=false; }

	_SNP[snp]->one[indi] = s1;
	_SNP[snp]->two[indi] = s2;
}


void snpdt::writeGenMapPed()
{

	string f = par::output_file_name + ".ped";

	gfun::printLOG("Writing recoded ped file to [ " + f + " ] \n");
	ofstream PED(f.c_str(), ios::out);
	PED.clear();

	string missingCode = par::out_missing_phenotype;

	for (int i=0;i<_sample.size() ;i++)
	{
		Individual * person = _sample[i];
		PED << person->fid << par::recode_delimit
			<< person->iid << par::recode_delimit
			<< person->pat << par::recode_delimit
			<< person->mat << par::recode_delimit
			<< person->sexcode << par::recode_delimit
			<< person->pheno_str;

		PED << "\n";
	}

	PED.close();

	//then write the gen file

	f = par::output_file_name + ".gen";

	gfun::printLOG("Writing recoded gen file to [ " + f + " ] \n");
	gfun::printLOG("Assuming a coding rule (0=A1A1, 1=A1A2, 2=A2A2, -9=missing)\n");
	ofstream GEN(f.c_str(), ios::out);
	GEN.clear();

	for (int i=0;i<_sample.size() ;i++)
	{
		if(_locus.size()>0) GEN <<genotypeToInt(i, 0);
		for (int l=1;l<_locus.size();l++)
			GEN <<par::recode_delimit<<genotypeToInt(i, l);

		GEN << "\n";
	}

	GEN.close();


	// And we also need a new map file

	
	f = par::output_file_name + ".map";
	gfun::printLOG( "Writing new map file to [ " + f + " ] \n");


	ofstream MAP(f.c_str(), ios::out);
	MAP.clear();

	for (int l=0;l<_locus.size();l++)
	{
		MAP << _locus[l]->chr << "\t"
			<< _locus[l]->name << "\t"
			<< _locus[l]->pos << "\t"
			<< _locus[l]->bp << "\t"
			<< _locus[l]->allele1 <<"\t"
			<< _locus[l]->allele2 <<"\n";
	}


	MAP.close();

}

void snpdt::writeGenMapPed(string file_root)
{

	string f = file_root + ".ped";

	gfun::printLOG("Writing recoded ped file to [ " + f + " ] \n");
	ofstream PED(f.c_str(), ios::out);
	PED.clear();

	string missingCode = par::out_missing_phenotype;

	for (int i=0;i<_sample.size() ;i++)
	{
		Individual * person = _sample[i];
		PED << person->fid << par::recode_delimit
			<< person->iid << par::recode_delimit
			<< person->pat << par::recode_delimit
			<< person->mat << par::recode_delimit
			<< person->sexcode << par::recode_delimit
			<< person->pheno_str;

		PED << "\n";
	}

	PED.close();

	//then write the gen file

	f = file_root + ".gen";

	gfun::printLOG("Writing recoded gen file to [ " + f + " ] \n");
	gfun::printLOG("Assuming a coding rule (0=A1A1, 1=A1A2, 2=A2A2, -9=missing)\n");
	ofstream GEN(f.c_str(), ios::out);
	GEN.clear();

	for (int i=0;i<_sample.size() ;i++)
	{
		if(_locus.size()>0) GEN <<genotypeToInt(i, 0);
		for (int l=1;l<_locus.size();l++)
			GEN <<par::recode_delimit<<genotypeToInt(i, l);

		GEN << "\n";
	}

	GEN.close();


	// And we also need a new map file


	f = file_root + ".map";
	gfun::printLOG( "Writing new map file to [ " + f + " ] \n");


	ofstream MAP(f.c_str(), ios::out);
	MAP.clear();

	for (int l=0;l<_locus.size();l++)
	{
		MAP << _locus[l]->chr << "\t"
			<< _locus[l]->name << "\t"
			<< _locus[l]->pos << "\t"
			<< _locus[l]->bp << "\t"
			<< _locus[l]->allele1 <<"\t"
			<< _locus[l]->allele2 <<"\n";
	}


	MAP.close();

}

void snpdt::readMtxFile()
{

}

void snpdt::readGenMapPed()
{
	gfun::checkFileExists(par::pedfile);
	gfun::checkFileExists(par::mapfile);
	gfun::checkFileExists(par::genfile);

	gfun::printLOG("Reading map (extended format) from [ " 
		+ par::mapfile + " ] \n");

	vector<Locus> ordered;

	ifstream MAP(par::mapfile.c_str(), ios::in);
	MAP.clear();

	int c=0;
	while(!MAP.eof())
	{

		Locus * loc = new Locus;

		MAP >> loc->chr   // will automatically by numeric
			>> loc->name 
			>> loc->pos   
			>> loc->bp     
			>> loc->allele1
			>> loc->allele2;

		if ( MAP.eof() )
		{
			delete loc; 
			continue;
		}
		else if ( MAP.fail() )
		{
			delete loc;
			gfun::error("Problem reading MAP file, line " + int2str(c+1) + "\n");
		}

		// Use the frequency slot temporarily to 
		// store order information
		loc->freq = c++;


		// Always included, but not always in correct order
		if (loc->name!="") 
		{
			_locus.push_back(loc);
			ordered.push_back(*loc);
		}  
		else
			delete loc;
	}


	gfun::printLOG( int2str(_locus.size()) 
		+ " markers to be included from [ " +
		par::mapfile + " ]\n");

	MAP.clear();
	MAP.close();

	if ( _locus.size() == 0 ) 
		gfun::shutdown();

	///////////////////////////////////////////////
	// Build ordered table, so that genotypes can 
	// be inserted in correct order; then swap locus 
	// file over

	// Sort vector of pointers to Locus
	stable_sort(_locus.begin(),_locus.end(),less<Locus*>());

	// Sort normal vector Locus
	stable_sort(ordered.begin(),ordered.end());

	c=0;
	for (int i=0; i<_locus.size(); i++)
	{
		// swap file order into locus position
		// and make all same chromosome
		ordered[i].bp = (int)ordered[i].freq;
		ordered[i].chr = 1;
		// keep track of genetic order, but only 
		ordered[i].freq = c++;
	}

	// resort to get lookup table
	stable_sort(ordered.begin(),ordered.end());  


	///////////////////////////////////////
	// Do we want to look at all the data?

	vector<int> include(0);
	int nl_actual = _locus.size();

	// If we do want to look at all the data
	for (int j=0; j<ordered.size(); j++)
		include.push_back((int)ordered[j].freq);

	


	//////////////////////////////
	// Read pedigree information

	// Initially, assume a binary trait
	par::qt = false;
	par::bt = true; 

	gfun::printLOG("Reading pedigree information from [ " 
		+ par::pedfile + " ] \n");

	ifstream PED;
	PED.open(par::pedfile.c_str());
	PED.clear();

	vector<Individual*> ambiguous;

	int nmale = 0;
	int nfemale = 0;
	int nambig = 0;

	c=0;
	while(!PED.eof())
	{

		Individual * person = new Individual;

		// First 6 obligatory fields
		string phenotype;

		PED >> person->fid 
			>> person->iid 
			>> person->pat
			>> person->mat
			>> person->sexcode
			>> person->pheno_str;

		phenotype=person->pheno_str;

		// Are we using 0/1 coding?
		if (par::coding01) 
		{
			if ( phenotype == "1" ) 
				phenotype = "2";      
			else if ( phenotype == "0" )
				phenotype = "1";
			else 
				phenotype = "0";
		}

		///set the binary phenotype
		if(phenotype == "2") person->aff = 1;
		else if (phenotype == "1") person->aff = 0;

		//set the family number
		if(person->fid == "-9") person->nfid= -9;
		else istringstream(person->fid) >> person->nfid;


		// Skip last empty line that gets read
		if (person->fid=="") 
		{
			delete person;
			break;
		}

		// Check sex
		if (person->sexcode=="1")
		{
			person->sex = 1; // male
			nmale++;
		}
		else if (person->sexcode=="2")
		{
			person->sex = 0;  // female (default)
			nfemale++;
		}
		else 
		{
			person->sex = par::missing_int;//missing
			nambig++;
		}

		// Increase person counter
		c++;

		// Add individual to list
		_sample.push_back(person);

	}

	PED.clear();
	PED.close();





	gfun::printLOG(int2str(c)+" individuals read from [ " 
		+ par::pedfile + " ] \n");

	if (par::bt) 
	{

		if (par::coding01) 
			gfun::printLOG("Assuming a disease phenotype (0=unaff, 1=aff, other=miss)\n");
		else
		{
			gfun::printLOG("Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)\n");
		}
	}


	//////////////////////////////
	// Allocate space for SNPs


	for (int i=0; i<nl_actual; i++)
	{
		CSNP * newlocus = new CSNP;
		newlocus->one.resize( _sample.size() );
		newlocus->two.resize( _sample.size() );
		_SNP.push_back(newlocus);
	}

	gfun::printLOG("Reading genotype information from [ " 
		+ par::genfile + " ] \n");

	gfun::printLOG("Assuming a coding rule (0=A1A1, 1=A1A2, 2=A2A2, -9=missing)\n");

	ifstream GEN;
	GEN.open(par::genfile.c_str());
	GEN.clear();

	for (c=0; c<_sample.size(); c++)
	{
		////now read genotype information
		vector<int> indx_tmp; indx_tmp.resize(_locus.size(),par::missing_int);

		for (int i=0; i<indx_tmp.size(); i++) GEN>>indx_tmp[i];

		if (!GEN) 
			gfun::error("Problem with the GEN file "+int2str(c)+"th line...has the MAP file been changed?\n");	

		for(int i=0; i<include.size(); i++){
			if(include[i]>-1) writeIntToGenotype(c,include[i],indx_tmp[i]);
		}

	}

	PED.clear();
	PED.close();


}

vector<int> snpdt::matchName(vector<string> names)
{
	vector<int> colidx; colidx.resize(names.size(),par::missing_int);

	for(int i=0; i<names.size(); i++){
		for(int j=0; j<_locus.size(); j++){
			if(_locus[j]->name == names[i]) {
				colidx[i]=j; break;
			}
		}
	}

	return colidx;
}

void snpdt::matchNameAllel(vector<string> names, vector<string> allels, vector<int> & colidx, vector<int> & geno_type)
{
	colidx.clear(); geno_type.clear();
	colidx.resize(names.size(),par::missing_int); geno_type.resize(allels.size(),par::missing_int);

	for(int i=0; i<names.size(); i++){
		for(int j=0; j<_locus.size(); j++){
			if(_locus[j]->name == names[i]) {
				colidx[i]=j; break;
			}
		}
	}

	for(int i=0; i<colidx.size(); i++){
		if(colidx[i]>=0) {
			geno_type[i]=genoStrToInt(colidx[i],string(1,allels[i][0]),string(1,allels[i][1]));
		}
	}
}