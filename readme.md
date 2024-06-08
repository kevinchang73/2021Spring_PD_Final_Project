# Macro Legalization

Macro legalization is an important step in today’s placement flow. In this work, a constraint-graph based macro legalization algorithm flow which combines iterative refinement and simulated annealing is applied to solve the macro legalization problem. Experimental result shows the flow could achieve a balance between cost and runtime.

Please visit https://kevinchang73.github.io/ for detailed information.

## HOW TO COMPILE:

Under the root directory, type: 
	`make`

## HOW TO RUN:

Under the root directory, type: 
	`./bin/macroLegalizer <DEF file> <LEF file> <constraint file> <output file>`

For example, `./bin/macroLegalizer inputs/case1/case1.def inputs/case1/case1.lef inputs/case1/case1.txt outputs/case1.out`
 
## HOW TO VISUALIZE:

After running the program, under the root directory, type: 
	`gnuplot result.plt`

## NOTICE:

1. If the submitted binary cannot run, type "make clean" first and then compile again.
2. The submitted version is optimized for the larger case (case1), thus the result of caseSample might be a little bit different with the one in the report. 
