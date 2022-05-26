=====
*This readme template is provided on the NTUEE Spring 2019 Algorithm Course.
=====
This is the README file for Physical Design Final Project
Author: B07901056 張凱鈞
Date: 2021.06.27
=====
SYNOPSIS:

This program tries to solve the macro legalization problem
=====
DIRECTORY: b07901056_final_project

bin/	  macroLegalizer (executable binary)
src/ 	  *.cpp, *.h (source C++ codes)
inputs/   input cases
outputs/  output results (the best results reported)
Makefile
readme.txt
*.plt     the initial and result placement for case 1
======
HOW TO COMPILE:

Under the b07901056_final_project/ directory, type: 
	make
======
HOW TO RUN:

Under the b07901056_final_project/ directory, type: 
	./bin/macroLegalizer <DEF file> <LEF file> <constraint file> <output file>

For example, ./bin/macroLegalizer inputs/case1/case1.def inputs/case1/case1.lef inputs/case1/case1.txt outputs/case1.out
 
======
HOW TO VISUALIZE:

After running the program, under the b07901056_final_project/ directory, type: 
	gnuplot result.plt

=====
NOTICE:

1. If the submitted binary cannot run, type "make clean" first and then compile again.
2. The submitted version is optimized for the larger case (case1), thus the result of caseSample might be a little bit different with the one in the report. 
