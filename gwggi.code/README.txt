PROGRAM: GWGGI

DESCRIPTION: Genome Wide Joint Association Software

AUTHOR: Changshuai Wei

CONTACT: weichangshuai@gmail.com

YEAR: 2011

LICENSE: Released under GNU General Public License, v2 (see
COPYING.txt)


COMPILATION: You will need a standard C/C++ compiler such as GNU gcc
(version 3).

Using TAMW:

gwggi --bfile example --tamw --converge --ntree 2000 --clsf-sqrt --tree-depth 4 --showntree --out ex.rst

Using LRMW:

gwggi --bfile example --lmw --showcv --maxauc 0.99 --out ex.rst


