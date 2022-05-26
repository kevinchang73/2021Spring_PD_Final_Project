#ifndef LEGALIZER_H
#define LEGALIZER_H

#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include "macro.h"
using namespace std;

class Legalizer
{
public:
    Legalizer() {}
    ~Legalizer() {}

    //input
    void parseInput(fstream& def, fstream& lef, fstream& constraint);
    //output
    void writeOutput(fstream& out);
    void plotResult(const string outfilename, bool isPrompt);
    //legalize
    void legalize();

private:
    //chip
    string _version;
    string _designName;
    int _unit;
    vector<pair<int, int> > _chipBoundary;
    int _chipLeft;
    int _chipBottom;
    int _chipRight;
    int _chipTop;
    double _powerplanWidth;
    double _minChannelSpacing;
    double _bufferAreaReservation;
    int _alpha, _beta;
    
    //macros / components
    size_t _numMacros;
    vector<Macro*> _macroList;
    vector<Macro*> _sortMacroList;
    unordered_map<string, pair<int, int> > _macroLib;
    
    //constraint graph
    vector<VertexXDP*> _graphXDPVertical;
    vector<VertexXDP*> _graphXDPHorizontal;
    vector<VertexXDP*> _topologicalSeqH;
    vector<VertexXDP*> _topologicalSeqV;
    vector< vector<bool> > _deTransitiveGraphVertical;
    vector< vector<bool> > _deTransitiveGraphHorizontal;

    //grids
    set<pair<int, int> > _noCellGrids;

    //cost
    double _bestCost;

    //initialization
    void _buildPsuedoMacro();
    void _initializeGraphXDP();
    EdgeXDP* _addEdge(VertexXDP* u, VertexXDP* v, int);
    void _delEdge(EdgeXDP* e);
    void _topologicalSortHorizontal(vector<VertexXDP*>&);
    void _topologicalSortVertical(vector<VertexXDP*>&);
    void _reportTopologicalSeq(bool);
    void _buildDeTransitiveGraph();
    //util
    void _plotBoxPLT( ofstream& stream, int x1, int y1, int x2, int y2 );
    void _plotHorizontalGraph(bool, string file_name);
    void _plotVerticalGraph(bool, string file_name);
    void _sortMacros();
    //iteration
    void _computeLongestPath(bool, int& length, vector<VertexXDP*>& path);
    void _reportLongestPath(const int& length, const vector<VertexXDP*>& path);
    void _updateLB(bool);
    void _updateRT(bool);
    void _determineMacroLocation();
    void _storeBestLocation();
    void _refineLongestPath(bool, const vector<VertexXDP*>& path, const int& limit = -1);
    void _resetBoundaryEdgeWeight();
    bool _checkBufferArea();
    void _buildCellArea();
    void _buildCellArea(Macro* m);
    bool _checkAround(Macro* m);
    //cost
    double _computeDisplacementCurrent() const;
    double _computeDisplacementBest() const;
    double _computeNoCellAreaCurrent() const;
    double _computeNoCellAreaBest() const;
    double _computeCostCurrent() const;
    double _computeCostBest() const;
    //simulated annealing
    void _simulatedAnnealing();
    void _perturb();
    void _swap(Macro* a, Macro* b);
    EdgeXDP* _reverse(EdgeXDP* e);
    EdgeXDP* _move(EdgeXDP* e);
    Macro* _randomMacro() const;
    EdgeXDP* _randomEdge() const;

    struct compare_macro
    {
        inline bool operator() (Macro* m1, Macro* m2)
        {
            return (m1->getSize() > m2->getSize());
        }
    };
    struct compare_macro_horizontal
    {
        inline bool operator() (Macro* m1, Macro* m2)
        {
            return (m1->_vertexXDPHorizontal->slack() < m2->_vertexXDPHorizontal->slack());
        }
    };
    struct compare_macro_vertical
    {
        inline bool operator() (Macro* m1, Macro* m2)
        {
            return (m1->_vertexXDPVertical->slack() < m2->_vertexXDPVertical->slack());
        }
    };
    struct compare_refine_horizontal
    {
        inline bool operator() (pair<VertexXDP*, int> p1, pair<VertexXDP*, int> p2)
        {
            return (p1.second > p2.second);
        }
    };

};

#endif  // LEGALIZER_H
