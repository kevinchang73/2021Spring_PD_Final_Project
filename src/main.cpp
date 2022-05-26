#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include "legalizer.h"
using namespace std;

int main(int argc, char** argv)
{
    fstream def, lef, constraint, output;
    if (argc == 5) {
        def.open(argv[1], ios::in);
        lef.open(argv[2], ios::in);
        constraint.open(argv[3], ios::in);
        output.open(argv[4], ios::out);
        if (!def) {
            cerr << "Cannot open the input file \"" << argv[2]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if (!lef) {
            cerr << "Cannot open the input file \"" << argv[3]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if (!constraint) {
            cerr << "Cannot open the input file \"" << argv[3]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if (!output) {
            cerr << "Cannot open the output file \"" << argv[4]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
    }
    else {
        cerr << "Usage: ./macroLegalization <def file> <lef file> " <<
                "<constraint file> <output file>" << endl;
        exit(1);
    }

    Legalizer* legal = new Legalizer();
    legal->parseInput(def, lef, constraint);
    legal->plotResult("init.plt", false);
    legal->legalize();
    legal->writeOutput(output);
    legal->plotResult("result.plt", false);
    cout << "\n";
    return 0;
}
