# gwggi

PROGRAM: GWGGI

DESCRIPTION: Genome Wide Joint Association Software

AUTHOR: Changshuai Wei

CONTACT: weichangshuai@gmail.com

YEAR: 2011

Precompiled excutable for windows and unix, as well as example data can be found at
http://changshuaiwei.github.io/software.html#GWGGI


COMPILATION: You will need a standard C/C++ compiler such as GNU gcc
(version 3).

Using TAMW:

gwggi --bfile example --tamw --converge --ntree 2000 --clsf-sqrt --tree-depth 4 --showntree --out ex.rst

Using LRMW:

gwggi --bfile example --lmw --showcv --maxauc 0.99 --out ex.rst


##Data format

GWGGI can read plink binary format, using --bfile filename.

GWGGI can also read text format, using --file filename.

 

##Likelihood Ratio Forword Mann-Whitney

Fast forward searching higher order joint association, using --lmw.

Example: gwggi --bfile filename --lmw --out resultfile.

More complicated example: gwggi --bfile filename --lmw --hz --hiorder 10 --maxauc 0.99 --out resultfile

 

##Trees Assembling Mann-Whitney

Assembling multiple trees to search low-marginal effect, using --tamw.

Example: gwggi --bfile filename --tamw --out resultfile.

More complicated example: gwggi --bfile filename --tamw --converge --ntree 10000 --td-burnin --showntree --out resultfile

 

##Options&Functions
###Functions

--bfile string, read binary file

--file string, read text file

--out string, file for output

--lmw, Likelihood Ratio Mann-Whitney

--lmw-apl string, apply lmw to a new data

--tamw, Trees Assembling Mann-Whitney

--tamw-apl, apply tamw to a new data

###Options

--hz, allow heterozygosity ( Aa v.s. aa,AA )

--hiorder int, define the highest order for searching, default is 10

--maxauc float, define the highest auc for searching, default is 0.9

--converge, show the converging status of tamw in *.monitor file

--ntree int, number of trees in tamw, default is 2000

--tree-depth int, set tree depth, default is 6

--td-burnin, use first 50 bootstrap sample to determine the trees depth

--showntree, show the number of tree built as program is running

--clsf int, number of genetic variants randomly selected in tamw

--clsf-sqrt, using sqrt of total genetic variant for --clsf in tamw, default is ln

--showcv, show the progress of the program when runing lmw

--no-cv, scaning without cross-validation for lmw

--lr-sub string, output the subject specific LR value

