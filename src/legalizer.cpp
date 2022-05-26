#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include "macro.h"
#include "legalizer.h"
using namespace std;

#define HORIZONTAL true
#define VERTICAL false

void Legalizer::parseInput(fstream& def, fstream& lef, fstream& constraint){
    string temp;
    getline(def, temp);
    temp = temp.substr(temp.find(' ')+1);
    _version = temp.substr(0, temp.find(' '));
    getline(def, temp);
    temp = temp.substr(temp.find(' ')+1);
    _designName = temp.substr(0, temp.find(' '));
    getline(def, temp);
    temp = temp.substr(temp.find(' ')+1);
    temp = temp.substr(temp.find(' ')+1);
    temp = temp.substr(temp.find(' ')+1);
    _unit = stoi(temp.substr(0, temp.find(' ')));
    getline(def, temp);

    while(!lef.eof()){
        getline(lef, temp);
        if(temp == "\n") continue;
        if(temp == "") continue;
        if(temp.substr(0,3) == "END") continue;
        if(temp.substr(0,5) == "MACRO"){
            string name = temp.substr(6);
            getline(lef, temp);
            while(temp.substr(0,1) == " "){
                temp = temp.substr(1);
            }
            temp = temp.substr(5);
            string s_width = temp.substr(0, temp.find(' '));
            temp = temp.substr(temp.find(' ')+4);
            string s_height = temp.substr(0, temp.find(' '));
            _macroLib[name] = pair<int, int>(_unit*stod(s_width), _unit*stod(s_height));
        }
    }
    
    int x_coord, y_coord;
    while(true){
        def >> temp;
        if(temp == "DIEAREA") continue;
        if(temp == ";") break;
        if(temp == "("){
            def >> x_coord;
            def >> y_coord;
            _chipBoundary.push_back(pair<int, int>(x_coord, y_coord));
            def >> temp;
        }
    }
    getline(def, temp);
    def >> temp;
    def >> _numMacros;
    def >> temp;
    while(true){
        def >> temp;
        if(temp == "END") break;
        if(temp == "-"){
            string name, macroName;
            bool fixed;
            double x, y;
            def >> name;
            def >> macroName;
            def >> temp;
            def >> temp;
            if(temp == "FIXED") fixed = true;
            else fixed = false;
            def >> temp;
            def >> x >> y;
            Macro* m = new Macro(name, macroName, _macroLib.find(macroName)->second.first, _macroLib.find(macroName)->second.second, x, y, fixed);
            _macroList.push_back(m);
            
        }
    }
    assert(_numMacros == _macroList.size());

    constraint >> temp;
    constraint >> _powerplanWidth;
    _powerplanWidth = _powerplanWidth * _unit;
    constraint >> temp;
    constraint >> _minChannelSpacing;
    _minChannelSpacing = _minChannelSpacing * _unit;
    constraint >> temp;
    constraint >> _bufferAreaReservation;
    _bufferAreaReservation = _bufferAreaReservation * _unit;
    constraint >> temp;
    constraint >> _alpha;
    constraint >> temp;
    constraint >> _beta;

    std::cout << _designName << "\n";
}

void Legalizer::legalize(){
    cout << "Import psuedo macros...\n";
    _buildPsuedoMacro();
    cout << "Initialize constraint graphs...\n";
    _initializeGraphXDP();
    _buildDeTransitiveGraph();
    _sortMacros();
    //_plotHorizontalGraph(true, "horizontal.graph");
    //_plotVerticalGraph(true, "vertical.graph");
    cout << "Topological sort...\n";
    _topologicalSortHorizontal(_topologicalSeqH);
    _topologicalSortVertical(_topologicalSeqV);
    
    cout << "Compute longest path...\n";
    int length_h, length_v;
    vector<VertexXDP*> path_h, path_v;
    _computeLongestPath(HORIZONTAL, length_h, path_h);
    _computeLongestPath(VERTICAL, length_v, path_v);
    //_reportLongestPath(length_h, path_h);
    //_reportLongestPath(length_v, path_v);
    cout << "Iteratively refine...\n";
    int count = 0;
    while((length_h > (_chipRight - _chipLeft)) || (length_v > (_chipTop - _chipBottom))){
        ++count;
        if(length_h <= (_chipRight - _chipLeft)){
            cout << "Violate vertical chip length\n";
            _refineLongestPath(VERTICAL, path_v);
            _topologicalSortHorizontal(_topologicalSeqH);
            _topologicalSortVertical(_topologicalSeqV);
            _computeLongestPath(VERTICAL, length_v, path_v);
            //_reportLongestPath(length_v, path_v);
            if(length_v <= (_chipTop - _chipBottom)){
                break;
            }
        }
        else if(length_v <= (_chipTop - _chipBottom)){
            cout << "Violate horizontal chip length\n";
            _refineLongestPath(HORIZONTAL, path_h);
            _topologicalSortHorizontal(_topologicalSeqH);
            _topologicalSortVertical(_topologicalSeqV);
            _computeLongestPath(HORIZONTAL, length_h, path_h);
            //_reportLongestPath(length_h, path_h);
            if(length_h <= (_chipRight - _chipLeft)){
                break;
            }
        }
        else{
            if((length_h - (_chipRight - _chipLeft)) >= (length_v - (_chipTop - _chipBottom))){ //TODO: difference or ratio
                cout << "Violate both chip lengths. Refine horizontal length first.\n";
                _refineLongestPath(HORIZONTAL, path_h, length_v);
                _topologicalSortHorizontal(_topologicalSeqH);
                _topologicalSortVertical(_topologicalSeqV);
                _computeLongestPath(HORIZONTAL, length_h, path_h);
                //_reportLongestPath(length_h, path_h);
            }
            else{
                cout << "Violate both chip lengths. Refine vertical length first.\n";
                _refineLongestPath(VERTICAL, path_v, length_h);
                _topologicalSortHorizontal(_topologicalSeqH);
                _topologicalSortVertical(_topologicalSeqV);
                _computeLongestPath(VERTICAL, length_v, path_v);
                //_reportLongestPath(length_v, path_v);
            }
        }
    }
    cout << "Number of Iterations = " << count << "\n";
    cout << "Determine macro location...\n";
    _updateLB(HORIZONTAL);
    _updateLB(VERTICAL);
    _updateRT(HORIZONTAL);
    _updateRT(VERTICAL);
    _determineMacroLocation();
    _storeBestLocation();
    cout << "Estimate cost...\n";
    double cost;
    cost = _computeCostBest();
    _bestCost = cost;
    plotResult("iterative.plt", false);
    //TODO: SA
    cout << "\n";
    _simulatedAnnealing();
    //_computeCostCurrent();
    cout << "Best cost = " << _bestCost << "\n";
    
    //_plotHorizontalGraph(false);
    //_plotVerticalGraph(false);
}

void Legalizer::_buildPsuedoMacro(){
    if(_chipBoundary.size() == 2){
        _chipLeft = _chipBoundary[0].first;
        _chipBottom = _chipBoundary[0].second;
        _chipRight = _chipBoundary[1].first;
        _chipTop = _chipBoundary[1].second;
    }
    else{
        _chipLeft = _chipBoundary[0].first;
        _chipRight = _chipBoundary[0].first;
        _chipBottom = _chipBoundary[0].second;
        _chipTop = _chipBoundary[0].second;
        for(unsigned i=1; i<_chipBoundary.size(); ++i){
            if(_chipBoundary[i].first > _chipRight){
                _chipRight = _chipBoundary[i].first;
            }
            else if(_chipBoundary[i].first < _chipLeft){
                _chipLeft = _chipBoundary[i].first;
            }
            if(_chipBoundary[i].second > _chipTop){
                _chipTop = _chipBoundary[i].second;
            }
            else if(_chipBoundary[i].second < _chipBottom){
                _chipBottom = _chipBoundary[i].second;
            }
        }
        
        for(unsigned i=0; i<_chipBoundary.size(); ++i){
            if(_chipBoundary[i].first != _chipLeft && _chipBoundary[i].first != _chipRight && _chipBoundary[i].second != _chipBottom && _chipBoundary[i].second != _chipTop){
                if(_chipBoundary[i].first == _chipBoundary[i-1].first){
                    if(_chipBoundary[i].second > _chipBoundary[i-1].second){
                        string name = "psuedo_" + to_string(i);
                        Macro* m = new Macro(name, "psuedo", (_chipBoundary[i].first - _chipBoundary[i+1].first), (_chipBoundary[i].second - _chipBoundary[i-1].second), _chipBoundary[i+1].first, _chipBoundary[i-1].second, true, true);
                        _macroList.push_back(m);
                    }
                    else{
                        string name = "psuedo_" + to_string(i);
                        Macro* m = new Macro(name, "psuedo", (_chipBoundary[i+1].first - _chipBoundary[i].first), (_chipBoundary[i-1].second - _chipBoundary[i].second), _chipBoundary[i].first, _chipBoundary[i].second, true, true);
                        _macroList.push_back(m);
                    }
                }
                else{
                    if(_chipBoundary[i].first > _chipBoundary[i-1].first){
                        string name = "psuedo_" + to_string(i);
                        Macro* m = new Macro(name, "psuedo", (_chipBoundary[i].first - _chipBoundary[i-1].first), (_chipBoundary[i+1].second - _chipBoundary[i].second), _chipBoundary[i-1].first, _chipBoundary[i-1].second, true, true);
                        _macroList.push_back(m);
                    }
                    else{
                        string name = "psuedo_" + to_string(i);
                        Macro* m = new Macro(name, "psuedo", (_chipBoundary[i-1].first - _chipBoundary[i].first), (_chipBoundary[i].second - _chipBoundary[i+1].second), _chipBoundary[i+1].first, _chipBoundary[i+1].second, true, true);
                        _macroList.push_back(m);
                    }
                }
            }
        }
    }
}

void Legalizer::_initializeGraphXDP(){
    //Vertex initialization
    for(unsigned i=0; i<_macroList.size(); ++i){
        VertexXDP* v = new VertexXDP(_macroList[i], true, i);
        _graphXDPVertical.push_back(v);
        _macroList[i]->_vertexXDPVertical = v;
        VertexXDP* h = new VertexXDP(_macroList[i], false, i);
        _graphXDPHorizontal.push_back(h);
        _macroList[i]->_vertexXDPHorizontal = h;
    }
    Macro* ml = new Macro("left", "boundary", 0, _chipTop - _chipBottom, _chipLeft, _chipBottom, true, true);
    Macro* mr = new Macro("right", "boundary", 0, _chipTop - _chipBottom, _chipRight, _chipBottom, true, true);
    Macro* mb = new Macro("bottom", "boundary", _chipLeft - _chipRight, 0, _chipLeft, _chipBottom, true, true);
    Macro* mt = new Macro("top", "boundary", _chipLeft - _chipRight, 0, _chipLeft, _chipTop, true, true);
    _macroList.push_back(ml);
    _macroList.push_back(mr);
    _macroList.push_back(mb);
    _macroList.push_back(mt);
    VertexXDP* source_v = new VertexXDP(mb, true, _macroList.size()-4);
    _graphXDPVertical.push_back(source_v);
    mb->_vertexXDPVertical = source_v;
    VertexXDP* source_h = new VertexXDP(ml, false, _macroList.size()-4);
    _graphXDPHorizontal.push_back(source_h);
    ml->_vertexXDPHorizontal = source_h;
    VertexXDP* sink_v = new VertexXDP(mt, true, _macroList.size()-3);
    _graphXDPVertical.push_back(sink_v);
    mt->_vertexXDPVertical = sink_v;
    VertexXDP* sink_h = new VertexXDP(mr, false, _macroList.size()-3);
    _graphXDPHorizontal.push_back(sink_h);
    mr->_vertexXDPHorizontal = sink_h;

    //Edge initialization
    for(unsigned i=0; i<_macroList.size()-4; ++i){
        VertexXDP* m_h = _macroList[i]->_vertexXDPHorizontal;
        VertexXDP* m_v = _macroList[i]->_vertexXDPVertical;
        if(_macroList[i]->isFixed()){
            _addEdge(source_h, m_h, _macroList[i]->getX() - _chipLeft);
            _addEdge(m_h, sink_h, _chipRight - _macroList[i]->getX());
            _addEdge(source_v, m_v, _macroList[i]->getY() - _chipBottom);
            _addEdge(m_v, sink_v, _chipTop - _macroList[i]->getY());
        }
        else{
            _addEdge(source_h, m_h, 0);
            _addEdge(m_h, sink_h, _macroList[i]->getWidth());
            _addEdge(source_v, m_v, 0);
            _addEdge(m_v, sink_v, _macroList[i]->getHeight());
            //_addEdge(source_h, m_h, 0+_powerplanWidth);
            //_addEdge(m_h, sink_h, _macroList[i]->getWidth()+_powerplanWidth);
            //_addEdge(source_v, m_v, 0);
            //_addEdge(m_v, sink_v, _macroList[i]->getHeight());
        }
    }
    //_minChannelSpacing = _powerplanWidth;
    for(unsigned i=0; i<_macroList.size()-4; ++i){
        for(unsigned j=i+1; j<_macroList.size()-4; ++j){
            Macro* m1_v = _macroList[i];
            Macro* m1_h = _macroList[i];
            Macro* m2_v = _macroList[j];
            Macro* m2_h = _macroList[j];
            int x1 = m1_h->getX(), y1 = m1_v->getY();
            int x2 = m2_h->getX(), y2 = m2_v->getY();
            int w1 = m1_h->getWidth(), h1 = m1_v->getHeight();
            int w2 = m2_h->getWidth(), h2 = m2_v->getHeight();
            if(x1 > x2) { swap(x1, x2); swap(w1, w2); swap(m1_h, m2_h); }
            if(y1 > y2) { swap(y1, y2); swap(h1, h2); swap(m1_v, m2_v); }
            //case 1: no overlap
            if((x1 + w1 < x2) && (y1 + h1 < y2)){ //TODO: consider min. channel spacing
                if((x2 - (x1 + w1)) > (y2 - (y1 + h1))){
                    if(m1_h->getType() == "psuedo" || m2_h->getType() == "psuedo"){
                        _addEdge(m1_h->_vertexXDPHorizontal, m2_h->_vertexXDPHorizontal, w1);
                    }
                    else{
                        _addEdge(m1_h->_vertexXDPHorizontal, m2_h->_vertexXDPHorizontal, w1+_powerplanWidth/*/w1+_minChannelSpacing*/);
                    }
                }
                else{
                    if(m1_v->getType() == "psuedo" || m2_v->getType() == "psuedo"){
                        _addEdge(m1_v->_vertexXDPVertical, m2_v->_vertexXDPVertical, h1);
                    }
                    else{
                        _addEdge(m1_v->_vertexXDPVertical, m2_v->_vertexXDPVertical, /*h1+_powerplanWidth/*/h1+_minChannelSpacing);
                    }
                }
            }
            //case 2: overlap
            else if((x1 + w1 >= x2) && (y1 + h1 >= y2)){
                if((x1 + w1 - x2) < (y1 + h1 - y2)){
                    if(m1_h->getType() == "psuedo" || m2_h->getType() == "psuedo"){
                        _addEdge(m1_h->_vertexXDPHorizontal, m2_h->_vertexXDPHorizontal, w1);
                    }
                    else{
                        _addEdge(m1_h->_vertexXDPHorizontal, m2_h->_vertexXDPHorizontal, w1+_powerplanWidth/*/w1+_minChannelSpacing*/);
                    }
                }
                else{
                    if(m1_v->getType() == "psuedo" || m2_v->getType() == "psuedo"){
                        _addEdge(m1_v->_vertexXDPVertical, m2_v->_vertexXDPVertical, h1);
                    }
                    else{
                        _addEdge(m1_v->_vertexXDPVertical, m2_v->_vertexXDPVertical, /*h1+_powerplanWidth/*/h1+_minChannelSpacing);
                    }
                }
            }
            //case 3: only horizontal overlap
            else if((x1 + w1 >= x2)){
                if(m1_v->getType() == "psuedo" || m2_v->getType() == "psuedo"){
                    _addEdge(m1_v->_vertexXDPVertical, m2_v->_vertexXDPVertical, h1);
                }
                else{
                    _addEdge(m1_v->_vertexXDPVertical, m2_v->_vertexXDPVertical, /*h1+_powerplanWidth/*/h1+_minChannelSpacing);
                }
            }
            //case 4: only vertical overlap
            else{
                if(m1_h->getType() == "psuedo" || m2_h->getType() == "psuedo"){
                    _addEdge(m1_h->_vertexXDPHorizontal, m2_h->_vertexXDPHorizontal, w1);
                }
                else{
                    _addEdge(m1_h->_vertexXDPHorizontal, m2_h->_vertexXDPHorizontal, w1+_powerplanWidth/*/w1+_minChannelSpacing*/);
                }
            }
        }
    }
}

EdgeXDP* Legalizer::_addEdge(VertexXDP* u, VertexXDP* v, int weight){
    EdgeXDP* e = new EdgeXDP(u, v, weight);
    u->outList.push_back(e);
    v->inList.push_back(e);
    return e;
}

void Legalizer::_delEdge(EdgeXDP* e){
    VertexXDP* u = e->u;
    VertexXDP* v = e->v;
    for(auto ite = u->outList.begin(); ite != u->outList.end(); ++ite){
        if(*ite == e){
            u->outList.erase(ite);
            break;
        }
    }
    for(auto ite = v->inList.begin(); ite != v->inList.end(); ++ite){
        if(*ite == e){
            v->inList.erase(ite);
            break;
        }
    }
    delete e;
}

void Legalizer::_plotHorizontalGraph(bool svg, string file_name){
    ofstream out(file_name, ofstream::out);
    out << "digraph {\n";
    //out << "concentrate = true;\n";
    for(unsigned i=0; i<_graphXDPHorizontal.size(); ++i){
        out << "\"" << _graphXDPHorizontal[i]->m->getName() << "\"" << " [label = <" << _graphXDPHorizontal[i]->m->getName() << ">];\n";
    }
    for(unsigned i=0; i<_graphXDPHorizontal.size(); ++i){
        VertexXDP* u = _graphXDPHorizontal[i];
        for(unsigned j=0; j<u->outList.size(); ++j){
            VertexXDP* v = u->outList[j]->v;
            out << "\"" << u->m->getName() << "\" -> \"" << v->m->getName() << "\" [label = <" << u->outList[j]->weight << ">] \n";
        }
    }
    out << "}";
    out.close();
    if(svg){
        string cmd = "dot -Tsvg -O " + file_name;
        system(cmd.c_str());
    }
}

void Legalizer::_plotVerticalGraph(bool svg, string file_name){
    ofstream out(file_name, ofstream::out);
    out << "digraph {\n";
    out << "concentrate = true;\n";
    for(unsigned i=0; i<_graphXDPVertical.size(); ++i){
        out << "\"" << _graphXDPVertical[i]->m->getName() << "\"" << " [label = <" << _graphXDPVertical[i]->m->getName() << ">];\n";
    }
    for(unsigned i=0; i<_graphXDPVertical.size(); ++i){
        VertexXDP* u = _graphXDPVertical[i];
        for(unsigned j=0; j<u->outList.size(); ++j){
            VertexXDP* v = u->outList[j]->v;
            out << "\"" << u->m->getName() << "\" -> \"" << v->m->getName() << "\" [label = <" << u->outList[j]->weight << ">] \n";
        }
    }
    out << "}";
    out.close();
    if(svg){
        string cmd = "dot -Tsvg -O " + file_name;;
        system(cmd.c_str());
    }
}

void Legalizer::_topologicalSortHorizontal(vector<VertexXDP*>& seq){
    seq.resize(0);
    size_t i = 0;
    vector<bool> removed(_graphXDPHorizontal.size(), false);
    vector<int> inListNum;
    for(unsigned i=0; i<_graphXDPHorizontal.size(); ++i){
        inListNum.push_back(_graphXDPHorizontal[i]->inList.size());
    }
    while(seq.size() < _graphXDPHorizontal.size()){
        if(inListNum[i]==0 && !removed[i]){
            seq.push_back(_graphXDPHorizontal[i]);
            for(unsigned j=0; j<_graphXDPHorizontal[i]->outList.size(); ++j){
                VertexXDP* endpoint = _graphXDPHorizontal[i]->outList[j]->v;
                --inListNum[endpoint->idx];
            }
            removed[i] = true;
        }
        if(i == _graphXDPHorizontal.size() - 1){
            i = 0;
        }
        else{
            ++i;
        }
    }
}

void Legalizer::_topologicalSortVertical(vector<VertexXDP*>& seq){
    seq.resize(0);
    size_t i = 0;
    vector<bool> removed(_graphXDPVertical.size(), false);
    vector<int> inListNum;
    for(unsigned i=0; i<_graphXDPVertical.size(); ++i){
        inListNum.push_back(_graphXDPVertical[i]->inList.size());
    }
    int count = 0;
    while(seq.size() < _graphXDPVertical.size()){
        if(inListNum[i]==0 && !removed[i]){
            seq.push_back(_graphXDPVertical[i]);
            for(unsigned j=0; j<_graphXDPVertical[i]->outList.size(); ++j){
                VertexXDP* endpoint = _graphXDPVertical[i]->outList[j]->v;
                --inListNum[endpoint->idx];
            }
            removed[i] = true;
        }
        if(i == _graphXDPVertical.size() - 1){
            i = 0;
        }
        else{
            ++i;
        }
        ++count;
    }
}

void Legalizer::_reportTopologicalSeq(bool horizontal){
    if(horizontal){
        for(unsigned i=0; i<_topologicalSeqH.size(); ++i){
            cout << _topologicalSeqH[i]->m->getName() << " ";
        }
        cout << "\n";
    }
    else{
        for(unsigned i=0; i<_topologicalSeqV.size(); ++i){
            cout << _topologicalSeqV[i]->m->getName() << " ";
        }
        cout << "\n";
    }
}

void Legalizer::_buildDeTransitiveGraph(){
    _deTransitiveGraphHorizontal.resize(0);
    _deTransitiveGraphVertical.resize(0);
    for(size_t i=0; i<_macroList.size()-2; ++i){
        _deTransitiveGraphHorizontal.push_back(vector<bool>(_macroList.size()-2, false));
        _deTransitiveGraphVertical.push_back(vector<bool>(_macroList.size()-2, false));
    }
    for(size_t i=0; i<_macroList.size()-2; ++i){
        VertexXDP* u = _macroList[i]->_vertexXDPHorizontal;
        for(size_t j=0; j<u->outList.size(); ++j){
            VertexXDP* v = u->outList[j]->v;
            _deTransitiveGraphHorizontal[u->idx][v->idx] = true;
            for(size_t k=0; k<v->outList.size(); ++k){
                VertexXDP* w = v->outList[k]->v;
                _deTransitiveGraphHorizontal[u->idx][w->idx] = false;
            }
        }
    }
    for(size_t i=0; i<_macroList.size()-4; ++i){
        VertexXDP* u = _macroList[i]->_vertexXDPVertical;
        for(size_t j=0; j<u->outList.size(); ++j){
            VertexXDP* v = u->outList[j]->v;
            _deTransitiveGraphVertical[u->idx][v->idx] = true;
            for(size_t k=0; k<v->outList.size(); ++k){
                VertexXDP* w = v->outList[k]->v;
                _deTransitiveGraphVertical[u->idx][w->idx] = false;
            }
        }
    }
    for(size_t i=_macroList.size()-2; i<_macroList.size(); ++i){
        VertexXDP* u = _macroList[i]->_vertexXDPVertical;
        for(size_t j=0; j<u->outList.size(); ++j){
            VertexXDP* v = u->outList[j]->v;
            _deTransitiveGraphVertical[u->idx][v->idx] = true;
            for(size_t k=0; k<v->outList.size(); ++k){
                VertexXDP* w = v->outList[k]->v;
                _deTransitiveGraphVertical[u->idx][w->idx] = false;
            }
        }
    }
}

void Legalizer::_computeLongestPath(bool horizontal, int& length, vector<VertexXDP*>& path){
    if(horizontal){
        vector<pair<int, VertexXDP*> > dist(_graphXDPHorizontal.size(), pair<int, VertexXDP*>(-1, NULL));
        VertexXDP* source = _topologicalSeqH[0];
        VertexXDP* sink = _topologicalSeqH[_topologicalSeqH.size()-1];
        dist[source->idx] = pair<int, VertexXDP*>(0, NULL);
        for(size_t i=0; i<_topologicalSeqH.size(); ++i){
            VertexXDP* current = _topologicalSeqH[i];
            for(size_t j=0; j<current->outList.size(); ++j){
                EdgeXDP* edge = current->outList[j];
                VertexXDP* target = edge->v;
                if(dist[current->idx].first + edge->weight > dist[target->idx].first){
                    dist[target->idx].first = dist[current->idx].first + edge->weight;
                    dist[target->idx].second = current;
                }
            }
        }
        //for(size_t i=0; i<_graphXDPHorizontal.size(); ++i){
        //    VertexXDP* current = _graphXDPHorizontal[i];
        //    current->LB = dist[i].first + _chipLeft;
        //    //assert(current->idx == (int)i);
        //}
        length = dist[sink->idx].first;
        path.resize(0);
        VertexXDP* backtrace = sink;
        while(backtrace != source){
            path.insert(path.begin(), backtrace);
            backtrace = dist[backtrace->idx].second;
        }
        path.insert(path.begin(), source);
    }
    else{
        vector<pair<int, VertexXDP*> > dist(_graphXDPVertical.size(), pair<int, VertexXDP*>(-1, NULL));
        VertexXDP* source = _topologicalSeqV[0];
        VertexXDP* sink = _topologicalSeqV[_topologicalSeqV.size()-1];
        dist[source->idx] = pair<int, VertexXDP*>(0, NULL);
        for(size_t i=0; i<_topologicalSeqV.size(); ++i){
            VertexXDP* current = _topologicalSeqV[i];
            for(size_t j=0; j<current->outList.size(); ++j){
                EdgeXDP* edge = current->outList[j];
                VertexXDP* target = edge->v;
                if(dist[current->idx].first + edge->weight > dist[target->idx].first){
                    dist[target->idx].first = dist[current->idx].first + edge->weight;
                    dist[target->idx].second = current;
                }
            }
        }
        //for(size_t i=0; i<_graphXDPVertical.size(); ++i){
        //    VertexXDP* current = _graphXDPVertical[i];
        //    current->LB = dist[i].first + _chipBottom;
        //    //assert(current->idx == (int)i);
        //}
        length = dist[sink->idx].first;
        path.resize(0);
        VertexXDP* backtrace = sink;
        while(backtrace != source){
            path.insert(path.begin(), backtrace);
            backtrace = dist[backtrace->idx].second;
        }
        path.insert(path.begin(), source);
    }
}

void Legalizer::_reportLongestPath(const int& length, const vector<VertexXDP*>& path){
    cout << "length = " << length << "\n";
    cout << "path: (";
    for(size_t i=0; i<path.size(); ++i){
        cout << path[i]->m->getName() << ", ";
    }
    cout << ")\n";
}

void Legalizer::_updateLB(bool horizontal){
    if(horizontal){
        vector<int> dist(_graphXDPHorizontal.size(), INT32_MIN);
        VertexXDP* source = _topologicalSeqH[0];
        dist[source->idx] = _chipLeft;
        for(size_t i=0; i<_topologicalSeqH.size(); ++i){
            VertexXDP* current = _topologicalSeqH[i];
            for(size_t j=0; j<current->outList.size(); ++j){
                EdgeXDP* edge = current->outList[j];
                VertexXDP* target = edge->v;
                if(dist[current->idx] + edge->weight > dist[target->idx]){
                    dist[target->idx] = dist[current->idx] + edge->weight;
                }
            }
        }
        for(size_t i=0; i<_graphXDPHorizontal.size(); ++i){
            VertexXDP* current = _graphXDPHorizontal[i];
            current->LB = dist[i];
            //assert(current->idx == (int)i);
        }
    }
    else{
        vector<int> dist(_graphXDPVertical.size(), INT32_MIN);
        VertexXDP* source = _topologicalSeqV[0];
        dist[source->idx] = _chipBottom;
        for(size_t i=0; i<_topologicalSeqV.size(); ++i){
            VertexXDP* current = _topologicalSeqV[i];
            for(size_t j=0; j<current->outList.size(); ++j){
                EdgeXDP* edge = current->outList[j];
                VertexXDP* target = edge->v;
                if(dist[current->idx] + edge->weight > dist[target->idx]){
                    dist[target->idx] = dist[current->idx] + edge->weight;
                }
            }
        }
        for(size_t i=0; i<_graphXDPVertical.size(); ++i){
            VertexXDP* current = _graphXDPVertical[i];
            current->LB = dist[i];
            //assert(current->idx == (int)i);
        }
    }
}

void Legalizer::_updateRT(bool horizontal){
    if(horizontal){
        vector<int> dist(_graphXDPHorizontal.size(), INT32_MAX);
        VertexXDP* sink = _topologicalSeqH[_topologicalSeqH.size()-1];
        dist[sink->idx] = _chipRight;
        for(int i=_topologicalSeqH.size()-1; i>=0; --i){
            VertexXDP* current = _topologicalSeqH[i];
            for(size_t j=0; j<current->inList.size(); ++j){
                EdgeXDP* edge = current->inList[j];
                VertexXDP* target = edge->u;
                if(dist[current->idx] - edge->weight < dist[target->idx]){
                    dist[target->idx] = dist[current->idx] - edge->weight;
                }
            }
        }
        for(size_t i=0; i<_graphXDPHorizontal.size(); ++i){
            VertexXDP* current = _graphXDPHorizontal[i];
            current->RT = dist[i];
            //assert(current->idx == (int)i);
        }
    }
    else{
        vector<int> dist(_graphXDPVertical.size(), INT32_MAX);
        VertexXDP* sink = _topologicalSeqV[_topologicalSeqV.size()-1];
        dist[sink->idx] = _chipTop;
        for(int i=_topologicalSeqV.size()-1; i>=0; --i){
            VertexXDP* current = _topologicalSeqV[i];
            for(size_t j=0; j<current->inList.size(); ++j){
                EdgeXDP* edge = current->inList[j];
                VertexXDP* target = edge->u;
                if(dist[current->idx] - edge->weight < dist[target->idx]){
                    dist[target->idx] = dist[current->idx] - edge->weight;
                }
            }
        }
        for(size_t i=0; i<_graphXDPVertical.size(); ++i){
            VertexXDP* current = _graphXDPVertical[i];
            current->RT = dist[i];
            //assert(current->idx == (int)i);
        }
    }
}

void Legalizer::_determineMacroLocation(){
/*    for(size_t i=0; i<_numMacros; ++i){
        if(!_macroList[i]->isFixed()){
            VertexXDP* vertex_h = _macroList[i]->_vertexXDPHorizontal;
            VertexXDP* vertex_v = _macroList[i]->_vertexXDPVertical;
            int init_x = _macroList[i]->getInitX();
            int init_y = _macroList[i]->getInitY();
            int x, y;
            if(vertex_h->LB <= init_x && init_x <= vertex_h->RT){
                x = init_x;
            }
            else if(init_x < vertex_h->LB){
                x = vertex_h->LB;
            }
            else if(init_x > vertex_h->RT){
                x = vertex_h->RT;
            }
            if(vertex_v->LB <= init_y && init_y <= vertex_v->RT){
                y = init_y;
            }
            else if(init_y < vertex_v->LB){
                y = vertex_v->LB;
            }
            else if(init_y > vertex_v->RT){
                y = vertex_v->RT;
            }
            _macroList[i]->setCurrentPos(x, y);
            _macroList[i]->storeBestPos();
            for(unsigned j=0; j<_macroList[i]->_vertexXDPHorizontal->outList.size(); ++j){
                EdgeXDP* edge = _macroList[i]->_vertexXDPHorizontal->outList[j];
                if(edge->v->m->getName() == "right"){
                    edge->weight = _chipRight - x;
                }
            }
            for(unsigned j=0; j<_macroList[i]->_vertexXDPHorizontal->inList.size(); ++j){
                EdgeXDP* edge = _macroList[i]->_vertexXDPHorizontal->inList[j];
                if(edge->u->m->getName() == "left"){
                    edge->weight = x - _chipLeft;
                }
            }
            for(unsigned j=0; j<_macroList[i]->_vertexXDPVertical->outList.size(); ++j){
                EdgeXDP* edge = _macroList[i]->_vertexXDPVertical->outList[j];
                if(edge->v->m->getName() == "top"){
                    edge->weight = _chipTop - y;
                }
            }
            for(unsigned j=0; j<_macroList[i]->_vertexXDPVertical->inList.size(); ++j){
                EdgeXDP* edge = _macroList[i]->_vertexXDPVertical->inList[j];
                if(edge->u->m->getName() == "bottom"){
                    edge->weight = y - _chipBottom;
                }
            }
            _updateLB(HORIZONTAL);
            _updateLB(VERTICAL);
            _updateRT(HORIZONTAL);
            _updateRT(VERTICAL);
            cout << "Move macro " << _macroList[i]->getName() << " from (" << _macroList[i]->getInitX() << ", "
                 << _macroList[i]->getInitY() << ") to (" << _macroList[i]->getBestX() << ", "
                 << _macroList[i]->getBestY() << ")\n";
            string name1 = "horizontal" + to_string(i);
            string name2 = "vertical" + to_string(i);
            _plotHorizontalGraph(false, name1);
            _plotVerticalGraph(false, name2);
        }
    }
*/
    _sortMacroList.resize(0);
    for(size_t i=0; i<_numMacros; ++i){
        _sortMacroList.push_back(_macroList[i]);
    }
    std::sort(_sortMacroList.begin(), _sortMacroList.end(), compare_macro_horizontal());
    for(size_t i=0; i<_numMacros; ++i){
        if(!_sortMacroList[i]->isFixed()){
            VertexXDP* vertex_h = _sortMacroList[i]->_vertexXDPHorizontal;
            int init_x = _sortMacroList[i]->getInitX();
            int curr_y = _sortMacroList[i]->getInitY();
            int x;
            if(vertex_h->LB <= init_x && init_x <= vertex_h->RT){
                x = init_x;
            }
            else if(init_x < vertex_h->LB){
                x = vertex_h->LB;
            }
            else if(init_x > vertex_h->RT){
                x = vertex_h->RT;
            }
            _sortMacroList[i]->setCurrentPos(x, curr_y);
            for(unsigned j=0; j<_sortMacroList[i]->_vertexXDPHorizontal->outList.size(); ++j){
                EdgeXDP* edge = _sortMacroList[i]->_vertexXDPHorizontal->outList[j];
                if(edge->v->m->getName() == "right"){
                    edge->weight = _chipRight - x;
                }
            }
            for(unsigned j=0; j<_sortMacroList[i]->_vertexXDPHorizontal->inList.size(); ++j){
                EdgeXDP* edge = _sortMacroList[i]->_vertexXDPHorizontal->inList[j];
                if(edge->u->m->getName() == "left"){
                    edge->weight = x - _chipLeft;
                }
            }
            _updateLB(HORIZONTAL);
            _updateRT(HORIZONTAL);
            //cout << "Move macro " << _sortMacroList[i]->getName() << " from (" << _sortMacroList[i]->getInitX() << ", "
            //     << _sortMacroList[i]->getInitY() << ") to (" << _sortMacroList[i]->getBestX() << ", "
            //     << _sortMacroList[i]->getBestY() << ")\n";
            //string name = "horizontal" + to_string(i);
            //_plotHorizontalGraph(false, name);
        }
    }
    std::sort(_sortMacroList.begin(), _sortMacroList.end(), compare_macro_vertical());
    for(size_t i=0; i<_numMacros; ++i){
        if(!_sortMacroList[i]->isFixed()){
            VertexXDP* vertex_v = _sortMacroList[i]->_vertexXDPVertical;
            int curr_x = _sortMacroList[i]->getX();
            int init_y = _sortMacroList[i]->getInitY();
            int y;
            if(vertex_v->LB <= init_y && init_y <= vertex_v->RT){
                y = init_y;
            }
            else if(init_y < vertex_v->LB){
                y = vertex_v->LB;
            }
            else if(init_y > vertex_v->RT){
                y = vertex_v->RT;
            }
            _sortMacroList[i]->setCurrentPos(curr_x, y);
            for(unsigned j=0; j<_sortMacroList[i]->_vertexXDPVertical->outList.size(); ++j){
                EdgeXDP* edge = _sortMacroList[i]->_vertexXDPVertical->outList[j];
                if(edge->v->m->getName() == "top"){
                    edge->weight = _chipTop - y;
                }
            }
            for(unsigned j=0; j<_sortMacroList[i]->_vertexXDPVertical->inList.size(); ++j){
                EdgeXDP* edge = _sortMacroList[i]->_vertexXDPVertical->inList[j];
                if(edge->u->m->getName() == "bottom"){
                    edge->weight = y - _chipBottom;
                }
            }
            _updateLB(VERTICAL);
            _updateRT(VERTICAL);
            //cout << "Move macro " << _sortMacroList[i]->getName() << " from (" << _sortMacroList[i]->getX() << ", "
            //     << _sortMacroList[i]->getInitY() << ") to (" << _sortMacroList[i]->getBestX() << ", "
            //     << _sortMacroList[i]->getBestY() << ")\n";
            //string name = "vertical" + to_string(i);
            //_plotVerticalGraph(false, name);
        }
    }

}

void Legalizer::_storeBestLocation(){
    for(size_t i=0; i<_numMacros; ++i){
        if(!_macroList[i]->isFixed()){
            _macroList[i]->storeBestPos();
        }
    }
}

void Legalizer::_refineLongestPath(bool horizontal, const vector<VertexXDP*>& path, const int& l){
    int length_limit;
    if(horizontal){
        if(l == -1){
            length_limit = _chipTop - _chipBottom;
        }
        else{
            length_limit = l;
        }
    }
    else{
        if(l == -1){
            length_limit = _chipLeft - _chipRight;
        }
        else{
            length_limit = l;
        }
    }
    if(horizontal){
        pair<VertexXDP*, VertexXDP*> max_pair;
        int max_slack = INT32_MIN;
        for(size_t i=1; i<path.size()-2; ++i){
            //cout << "Check (" << path[i]->m->getName() << ", " << path[i+1]->m->getName() << ").\n";
            EdgeXDP* e;
            for(size_t j=0; j<path[i]->outList.size(); ++j){
                if(path[i]->outList[j]->v == path[i+1]){
                    e = path[i]->outList[j];
                    break;
                }
            }
            _delEdge(e);

            VertexXDP* u = path[i]->m->_vertexXDPVertical;
            VertexXDP* v = path[i+1]->m->_vertexXDPVertical;
            EdgeXDP* e_vertical;
            if(u->m->getY() <= v->m->getY()){
                e_vertical = _addEdge(u, v, u->m->getHeight()+_minChannelSpacing);
            }
            else{
                e_vertical = _addEdge(v, u, v->m->getHeight()+_minChannelSpacing);
            }
            
            int vertical_length;
            vector<VertexXDP*> vertical_path;
            _computeLongestPath(VERTICAL, vertical_length, vertical_path);
            if(vertical_length <= (length_limit)){
                //cout << "Move (" << path[i]->m->getName() << ", " << path[i+1]->m->getName() << ") to vertical is available.\n";
                if(e_vertical->v->m->getY()-(e_vertical->u->m->getY()+e_vertical->u->m->getHeight()) > max_slack){
                    max_pair = pair<VertexXDP*, VertexXDP*>(path[i], path[i+1]);
                    max_slack = e_vertical->v->m->getY()-(e_vertical->u->m->getY()+e_vertical->u->m->getHeight());
                }
            }
            _delEdge(e_vertical);
            _addEdge(path[i], path[i+1], path[i]->m->getWidth()+_powerplanWidth);
        }
        
        if(max_slack == INT32_MIN){
            cout << "Unable to refine\n";
            return;
        }
        EdgeXDP* e;
        for(size_t j=0; j<max_pair.first->outList.size(); ++j){
            if(max_pair.first->outList[j]->v == max_pair.second){
                e = max_pair.first->outList[j];
                break;
            }
        }
        _delEdge(e);
        VertexXDP* u = max_pair.first->m->_vertexXDPVertical;
        VertexXDP* v = max_pair.second->m->_vertexXDPVertical;
        //EdgeXDP* e_vertical;
        if(u->m->getY() <= v->m->getY()){
            _addEdge(u, v, u->m->getHeight()+_minChannelSpacing);
        }
        else{
            _addEdge(v, u, v->m->getHeight()+_minChannelSpacing);
        }
        cout << "Move (" << max_pair.first->m->getName() << ", " << max_pair.second->m->getName() << ") to vertical.\n";
    }
    else{
        pair<VertexXDP*, VertexXDP*> max_pair;
        int max_slack = INT32_MIN;
        for(size_t i=1; i<path.size()-2; ++i){
            //cout << "Check (" << path[i]->m->getName() << ", " << path[i+1]->m->getName() << ").\n";
            EdgeXDP* e;
            for(size_t j=0; j<path[i]->outList.size(); ++j){
                if(path[i]->outList[j]->v == path[i+1]){
                    e = path[i]->outList[j];
                    break;
                }
            }
            _delEdge(e);

            VertexXDP* u = path[i]->m->_vertexXDPHorizontal;
            VertexXDP* v = path[i+1]->m->_vertexXDPHorizontal;
            EdgeXDP* e_horizontal;
            if(u->m->getX() <= v->m->getX()){
                e_horizontal = _addEdge(u, v, u->m->getWidth()+_powerplanWidth);
            }
            else{
                e_horizontal = _addEdge(v, u, v->m->getWidth()+_powerplanWidth);
            }
            
            int horizontal_length;
            vector<VertexXDP*> horizontal_path;
            _computeLongestPath(HORIZONTAL, horizontal_length, horizontal_path);
            if(horizontal_length <= (length_limit)){
                cout << "Move (" << path[i]->m->getName() << ", " << path[i+1]->m->getName() << ") to horizontal is available.\n";
                if(e_horizontal->v->m->getX()-(e_horizontal->u->m->getX()+e_horizontal->u->m->getWidth()) > max_slack){
                    max_pair = pair<VertexXDP*, VertexXDP*>(path[i], path[i+1]);
                    max_slack = e_horizontal->v->m->getX()-(e_horizontal->u->m->getX()+e_horizontal->u->m->getWidth());
                }
            }
            else{
                cout << "Move (" << path[i]->m->getName() << ", " << path[i+1]->m->getName() << ") to horizontal is not available.\n";
            }
            _delEdge(e_horizontal);
            _addEdge(path[i], path[i+1], path[i]->m->getHeight()+_minChannelSpacing);
        }
        
        if(max_slack == INT32_MIN){
            cout << "Unable to refine\n";
            return;
        }
        EdgeXDP* e;
        for(size_t j=0; j<max_pair.first->outList.size(); ++j){
            if(max_pair.first->outList[j]->v == max_pair.second){
                e = max_pair.first->outList[j];
                break;
            }
        }
        _delEdge(e);
        VertexXDP* u = max_pair.first->m->_vertexXDPHorizontal;
        VertexXDP* v = max_pair.second->m->_vertexXDPHorizontal;
        //EdgeXDP* e_horizontal;
        if(u->m->getX() <= v->m->getX()){
            _addEdge(u, v, u->m->getWidth()+_powerplanWidth);
        }
        else{
            _addEdge(v, u, v->m->getWidth()+_powerplanWidth);
        }
        cout << "Move (" << max_pair.first->m->getName() << ", " << max_pair.second->m->getName() << ") to horizontal.\n";
    }
}

void Legalizer::_resetBoundaryEdgeWeight(){
    for(size_t i=0; i<_numMacros; ++i){
        if(!_macroList[i]->isFixed()){
            VertexXDP* v_h = _macroList[i]->_vertexXDPHorizontal;
            VertexXDP* v_v = _macroList[i]->_vertexXDPVertical;
            for(size_t j=0; j<v_h->outList.size(); ++j){
                if(v_h->outList[j]->v->m->getName() == "right"){
                    v_h->outList[j]->weight = _macroList[i]->getWidth();
                    break;
                }
            }
            for(size_t j=0; j<v_h->inList.size(); ++j){
                if(v_h->inList[j]->u->m->getName() == "left"){
                    v_h->inList[j]->weight = 0;
                    break;
                }
            }
            for(size_t j=0; j<v_v->outList.size(); ++j){
                if(v_v->outList[j]->v->m->getName() == "top"){
                    v_v->outList[j]->weight = _macroList[i]->getHeight();
                    break;
                }
            }
            for(size_t j=0; j<v_v->inList.size(); ++j){
                if(v_v->inList[j]->u->m->getName() == "bottom"){
                    v_v->inList[j]->weight = 0;
                    break;
                }
            }
        }
    }
}

bool Legalizer::_checkBufferArea(){
    vector<Macro*> _violatingList;
    for(size_t i=0; i<_numMacros; ++i){
        
    }
    return true;
}

void Legalizer::_buildCellArea(){
    for(size_t i=0; i<_macroList.size()-4; ++i){
        _buildCellArea(_macroList[i]);
    }
}

void Legalizer::_buildCellArea(Macro* m){
    for(int x=m->getBestX(); x<m->getBestX()+m->getWidth(); x+=500){
        for(int y=m->getBestY(); y<m->getBestX()+m->getHeight(); y+=500){
            _noCellGrids.insert(pair<int, int>(x, y));
        }
    }
}

double Legalizer::_computeDisplacementCurrent() const{
    double dx = 0, dy = 0;
    for(size_t i=0; i<_numMacros; ++i){
        dx += (double)(abs(_macroList[i]->getX() - _macroList[i]->getInitX())) / (double)_unit;
        dy += (double)(abs(_macroList[i]->getY() - _macroList[i]->getInitY())) / (double)_unit;
    }
    //cout << "Displacement = " << dx + dy << "\n";
    return dx + dy;
}

double Legalizer::_computeDisplacementBest() const{
    double dx = 0, dy = 0;
    for(size_t i=0; i<_numMacros; ++i){
        dx += (double)(abs(_macroList[i]->getBestX() - _macroList[i]->getInitX())) / (double)_unit;
        dy += (double)(abs(_macroList[i]->getBestY() - _macroList[i]->getInitY())) / (double)_unit;
    }
    cout << "Displacement = " << dx + dy << "\n";
    return dx + dy;
}

double Legalizer::_computeNoCellAreaCurrent() const{
    double noPlaceAreaHorizontal = 0, noPlaceAreaVertical = 0;
    for(size_t i=0; i<_macroList.size()-2; ++i){
        VertexXDP* u = _macroList[i]->_vertexXDPHorizontal;
        for(size_t j=0; j<u->outList.size(); ++j){
            VertexXDP* v = u->outList[j]->v;
            if(_deTransitiveGraphHorizontal[u->idx][v->idx]){
                int x1 = u->m->getX() + u->m->getWidth();
                int x2 = v->m->getX();
                int y11 = u->m->getY();
                int y12 = u->m->getY() + u->m->getHeight();
                int y21 = v->m->getY();
                int y22 = v->m->getY() + v->m->getHeight();
                if(y12 <= y21 || y22 <= y11){
                    continue;
                }
                if(x2 - x1 < _powerplanWidth){
                    if(y11 <= y21){
                        if(y12 <= y22){
                            noPlaceAreaHorizontal += ((double)(y12 - y21)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit));
                        }
                        else{
                            noPlaceAreaHorizontal += ((double)(y22 - y21)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit));
                        }
                    }
                    else{
                        if(y12 <= y22){
                            noPlaceAreaHorizontal += ((double)(y12 - y11)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit));
                        }
                        else{
                            noPlaceAreaHorizontal += ((double)(y22 - y11)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit));
                        }
                    }
                }
            }
        }
    }
    for(size_t i=0; i<_macroList.size()-4; ++i){
        VertexXDP* u = _macroList[i]->_vertexXDPVertical;
        for(size_t j=0; j<u->outList.size(); ++j){
            VertexXDP* v = u->outList[j]->v;
            if(_deTransitiveGraphVertical[u->idx][v->idx]){
                int y1 = u->m->getY() + u->m->getHeight();
                int y2 = v->m->getY();
                int x11 = u->m->getX();
                int x12 = u->m->getX() + u->m->getWidth();
                int x21 = v->m->getX();
                int x22 = v->m->getX() + v->m->getWidth();
                if(x12 <= x21 || x22 <= x11){
                    continue;
                }
                if(y2 - y1 < _powerplanWidth){
                    if(x11 <= x21){
                        if(x12 <= x22){
                            noPlaceAreaVertical += ((double)(x12 - x21)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                        else{
                            noPlaceAreaVertical += ((double)(x22 - x21)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                    }
                    else{
                        if(x12 <= x22){
                            noPlaceAreaVertical += ((double)(x12 - x11)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                        else{
                            noPlaceAreaVertical += ((double)(x22 - x11)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                    }
                }
            }
        }
    }
    for(size_t i=_macroList.size()-2; i<_macroList.size(); ++i){
        VertexXDP* u = _macroList[i]->_vertexXDPVertical;
        for(size_t j=0; j<u->outList.size(); ++j){
            VertexXDP* v = u->outList[j]->v;
            if(_deTransitiveGraphVertical[u->idx][v->idx]){
                int y1 = u->m->getY() + u->m->getHeight();
                int y2 = v->m->getY();
                int x11 = u->m->getX();
                int x12 = u->m->getX() + u->m->getWidth();
                int x21 = v->m->getX();
                int x22 = v->m->getX() + v->m->getWidth();
                if(x12 <= x21 || x22 <= x11){
                    continue;
                }
                if(y2 - y1 < _powerplanWidth){
                    if(x11 <= x21){
                        if(x12 <= x22){
                            noPlaceAreaVertical += ((double)(x12 - x21)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                        else{
                            noPlaceAreaVertical += ((double)(x22 - x21)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                    }
                    else{
                        if(x12 <= x22){
                            noPlaceAreaVertical += ((double)(x12 - x11)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                        else{
                            noPlaceAreaVertical += ((double)(x22 - x11)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                    }
                }
            }
        }
    }
    //cout << "Horizontal no place area = " << noPlaceAreaHorizontal << "\n";
    //cout << "Vertical no place area = " << noPlaceAreaVertical << "\n";
    return noPlaceAreaHorizontal + noPlaceAreaVertical;
}

double Legalizer::_computeNoCellAreaBest() const{
    double noPlaceAreaHorizontal = 0, noPlaceAreaVertical = 0;
    for(size_t i=0; i<_macroList.size()-2; ++i){
        VertexXDP* u = _macroList[i]->_vertexXDPHorizontal;
        for(size_t j=0; j<u->outList.size(); ++j){
            VertexXDP* v = u->outList[j]->v;
            if(_deTransitiveGraphHorizontal[u->idx][v->idx]){
                int x1 = u->m->getBestX() + u->m->getWidth();
                int x2 = v->m->getBestX();
                int y11 = u->m->getBestY();
                int y12 = u->m->getBestY() + u->m->getHeight();
                int y21 = v->m->getBestY();
                int y22 = v->m->getBestY() + v->m->getHeight();
                if(y12 <= y21 || y22 <= y11){
                    continue;
                }
                if(x2 - x1 < _powerplanWidth){
                    if(y11 <= y21){
                        if(y12 <= y22){
                            noPlaceAreaHorizontal += ((double)(y12 - y21)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit));
                            if(((double)(y12 - y21)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit)) > 0){
                                //cout << u->m->getName() << " and " << v->m->getName() << " area " << ((double)(y12 - y21)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit)) << "\n";
                            }
                        }
                        else{
                            noPlaceAreaHorizontal += ((double)(y22 - y21)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit));
                            if(((double)(y22 - y21)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit)) > 0){
                                //cout << u->m->getName() << " and " << v->m->getName() << " area " << ((double)(y22 - y21)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit)) << "\n";
                            }
                        }
                    }
                    else{
                        if(y12 <= y22){
                            noPlaceAreaHorizontal += ((double)(y12 - y11)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit));
                            if(((double)(y12 - y11)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit)) > 0){
                                //cout << u->m->getName() << " and " << v->m->getName() << " area " << ((double)(y12 - y11)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit)) << "\n";
                            }
                        }
                        else{
                            noPlaceAreaHorizontal += ((double)(y22 - y11)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit));
                            if(((double)(y22 - y11)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit)) > 0){
                                //cout << u->m->getName() << " and " << v->m->getName() << " area " << ((double)(y22 - y11)/(double)(_unit)) * ((double)(x2 - x1)/(double)(_unit)) << "\n";
                            }
                        }
                    }
                }
            }
        }
    }
    for(size_t i=0; i<_macroList.size()-4; ++i){
        VertexXDP* u = _macroList[i]->_vertexXDPVertical;
        for(size_t j=0; j<u->outList.size(); ++j){
            VertexXDP* v = u->outList[j]->v;
            if(_deTransitiveGraphVertical[u->idx][v->idx]){
                int y1 = u->m->getBestY() + u->m->getHeight();
                int y2 = v->m->getBestY();
                int x11 = u->m->getBestX();
                int x12 = u->m->getBestX() + u->m->getWidth();
                int x21 = v->m->getBestX();
                int x22 = v->m->getBestX() + v->m->getWidth();
                if(x12 <= x21 || x22 <= x11){
                    continue;
                }
                if(y2 - y1 < _powerplanWidth){
                    if(x11 <= x21){
                        if(x12 <= x22){
                            noPlaceAreaVertical += ((double)(x12 - x21)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                        else{
                            noPlaceAreaVertical += ((double)(x22 - x21)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                    }
                    else{
                        if(x12 <= x22){
                            noPlaceAreaVertical += ((double)(x12 - x11)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                        else{
                            noPlaceAreaVertical += ((double)(x22 - x11)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                    }
                }
            }
        }
    }
    for(size_t i=_macroList.size()-2; i<_macroList.size(); ++i){
        VertexXDP* u = _macroList[i]->_vertexXDPVertical;
        for(size_t j=0; j<u->outList.size(); ++j){
            VertexXDP* v = u->outList[j]->v;
            if(_deTransitiveGraphVertical[u->idx][v->idx]){
                int y1 = u->m->getBestY() + u->m->getHeight();
                int y2 = v->m->getBestY();
                int x11 = u->m->getBestX();
                int x12 = u->m->getBestX() + u->m->getWidth();
                int x21 = v->m->getBestX();
                int x22 = v->m->getBestX() + v->m->getWidth();
                if(x12 <= x21 || x22 <= x11){
                    continue;
                }
                if(y2 - y1 < _powerplanWidth){
                    if(x11 <= x21){
                        if(x12 <= x22){
                            noPlaceAreaVertical += ((double)(x12 - x21)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                        else{
                            noPlaceAreaVertical += ((double)(x22 - x21)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                    }
                    else{
                        if(x12 <= x22){
                            noPlaceAreaVertical += ((double)(x12 - x11)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                        else{
                            noPlaceAreaVertical += ((double)(x22 - x11)/(double)(_unit)) * ((double)(y2 - y1)/(double)(_unit));
                        }
                    }
                }
            }
        }
    }
    cout << "Horizontal no place area = " << noPlaceAreaHorizontal << "\n";
    cout << "Vertical no place area = " << noPlaceAreaVertical << "\n";
    return noPlaceAreaHorizontal + noPlaceAreaVertical;
}

double Legalizer::_computeCostCurrent() const{
    double d_cost, a_cost;
    d_cost = _computeDisplacementCurrent();
    a_cost = _computeNoCellAreaCurrent();
    cout << "Total cost = " << _alpha * d_cost + _beta * sqrt(a_cost) << "\n";
    return _alpha * d_cost + _beta * sqrt(a_cost);
}

double Legalizer::_computeCostBest() const{
    double d_cost, a_cost;
    d_cost = _computeDisplacementBest();
    a_cost = _computeNoCellAreaBest();
    cout << "Total cost = " << _alpha * d_cost + _beta * sqrt(a_cost) << "\n";
    return _alpha * d_cost + _beta * sqrt(a_cost);
}

void Legalizer::_simulatedAnnealing(){
    double orig_cost;
    double new_cost;
    int consecutive_uphill_count = 0;
    int length_h, length_v;
    vector<VertexXDP*> path_h, path_v;
    //double temperature = _computeCostCurrent();
    double prob, random_prob;
    double random;
    double ratio_swap = 0.0;
    double ratio_reverse = 0.0;
    int count = 0;
    while(consecutive_uphill_count < 3){
        cout << "\n";
        ++count;
        cout << "Iteration " << count << "\n";
        Macro* m1, *m2;
        EdgeXDP* e, *new_edge;
        random = (double)rand()/(RAND_MAX+1.0);
        cout << "Best cost = " << _bestCost << "\n";
        cout << "Original ";
        orig_cost = _computeCostCurrent();
        //perturb
        if(random < ratio_swap){
            m1 = _randomMacro();
            m2 = _randomMacro();
            while(m1 == m2){
                m2 = _randomMacro();
            }
            _swap(m1, m2);
        }
        else if(random < ratio_swap + ratio_reverse){
            e = _randomEdge();
            new_edge = _reverse(e);
        }
        else{
            e = _randomEdge();
            new_edge = _move(e);
        }
        //reset
        _resetBoundaryEdgeWeight();
        _buildDeTransitiveGraph();
        _topologicalSortHorizontal(_topologicalSeqH);
        _topologicalSortVertical(_topologicalSeqV);
        //evaluation
        _computeLongestPath(HORIZONTAL, length_h, path_h);
        _computeLongestPath(VERTICAL, length_v, path_v);
        //_reportLongestPath(length_h, path_h);
        //_reportLongestPath(length_v, path_v);
        if(!(length_h > (_chipRight - _chipLeft)) && !(length_v > (_chipTop - _chipBottom))){
            _updateLB(HORIZONTAL);
            _updateLB(VERTICAL);
            _updateRT(HORIZONTAL);
            _updateRT(VERTICAL);
            _determineMacroLocation();
            cout << "New ";
            new_cost = _computeCostCurrent();
            if(new_cost <= orig_cost){
                cout << "Accept directly\n";
                if(new_cost < _bestCost){
                    _storeBestLocation();
                    _bestCost = new_cost;
                }
                consecutive_uphill_count = 0;
                continue;
            }
            else{
                ++consecutive_uphill_count;
                //cout << "delta cost = " << new_cost - orig_cost << " temperature = " << temperature << "\n";
                prob = 0.12;//exp(-100*(new_cost - orig_cost)/temperature);
                random_prob = (double)rand()/(RAND_MAX+1.0);
                if(random_prob <= prob){
                    cout << "Accept with probability \n";
                    continue;
                }
                cout << "Reject with probability \n";
            }
        }
        else{
            cout << "Reject directly\n";
        }
        //recover
        //cout << "Recover...\n";
        if(random < ratio_swap){
            _swap(m1, m2);
        }
        else if(random < ratio_swap + ratio_reverse){
            _reverse(new_edge);
        }
        else{
            _move(new_edge);
        }
        _resetBoundaryEdgeWeight();
        _buildDeTransitiveGraph();
        _topologicalSortHorizontal(_topologicalSeqH);
        _topologicalSortVertical(_topologicalSeqV);
        _updateLB(HORIZONTAL);
        _updateLB(VERTICAL);
        _updateRT(HORIZONTAL);
        _updateRT(VERTICAL);
        _determineMacroLocation();
    }
    _resetBoundaryEdgeWeight();
    _buildDeTransitiveGraph();
    _topologicalSortHorizontal(_topologicalSeqH);
    _topologicalSortVertical(_topologicalSeqV);
    _updateLB(HORIZONTAL);
    _updateLB(VERTICAL);
    _updateRT(HORIZONTAL);
    _updateRT(VERTICAL);
    _determineMacroLocation();
}

void Legalizer::_perturb(){
    double random = (double)rand()/(RAND_MAX+1.0);
    double ratio_swap = 0.1;
    double ratio_reverse = 0.2;
    if(random < ratio_swap){
        Macro* m1 = _randomMacro();
        Macro* m2 = _randomMacro();
        while(m1 == m2){
            m2 = _randomMacro();
        }
        _swap(m1, m2);
    }
    else if(random < ratio_swap + ratio_reverse){
        EdgeXDP* e = _randomEdge();
        _reverse(e);
    }
    else{
        EdgeXDP* e = _randomEdge();
        _move(e);
    }
}

void Legalizer::_swap(Macro* a, Macro* b){
    VertexXDP* a_h = a->_vertexXDPHorizontal;
    VertexXDP* a_v = a->_vertexXDPVertical;
    int a_idx = a->_vertexXDPHorizontal->idx;
    int b_idx = b->_vertexXDPHorizontal->idx;
    a->_vertexXDPHorizontal = b->_vertexXDPHorizontal;
    a->_vertexXDPVertical = b->_vertexXDPVertical;
    b->_vertexXDPHorizontal = a_h;
    b->_vertexXDPVertical = a_v;
    
    a->_vertexXDPHorizontal->idx = a_idx;
    a->_vertexXDPVertical->idx = a_idx;
    a->_vertexXDPHorizontal->m = a;
    a->_vertexXDPVertical->m = a;
    b->_vertexXDPHorizontal->idx = b_idx;
    b->_vertexXDPVertical->idx = b_idx;
    b->_vertexXDPHorizontal->m = b;
    b->_vertexXDPVertical->m = b;
    cout << "Swap " << a->getName() << " and " << b->getName() << "\n";
}

EdgeXDP* Legalizer::_reverse(EdgeXDP* e){
    VertexXDP* u = e->u;
    VertexXDP* v = e->v;
    EdgeXDP* new_edge;
    if(u == u->m->_vertexXDPHorizontal){
        new_edge = _addEdge(v, u, v->m->getWidth() + _powerplanWidth);
        cout << "Reverse horizontal edge (" << e->u->m->getName() << ", " << e->v->m->getName() << ")\n";
    }
    else{
        assert(u == u->m->_vertexXDPVertical);
        new_edge = _addEdge(v, u, v->m->getHeight() + _minChannelSpacing);
        cout << "Reverse vertical edge (" << e->u->m->getName() << ", " << e->v->m->getName() << ")\n";
    }
    _delEdge(e);
    return new_edge;
}

EdgeXDP* Legalizer::_move(EdgeXDP* e){
    VertexXDP* u = (e->u == e->u->m->_vertexXDPHorizontal)? e->u->m->_vertexXDPVertical : e->u->m->_vertexXDPHorizontal;
    VertexXDP* v = (e->u == e->u->m->_vertexXDPHorizontal)? e->v->m->_vertexXDPVertical : e->v->m->_vertexXDPHorizontal;
    EdgeXDP* new_edge;
    if(u == u->m->_vertexXDPHorizontal){
        if(u->m->getX() > v->m->getX()){
            new_edge = _addEdge(v, u, v->m->getWidth() + _powerplanWidth);
            cout << "Move vertical edge (" << e->u->m->getName() << ", " << e->v->m->getName() << ") to horizontal edge (" << e->v->m->getName() << ", " << e->u->m->getName() << ")\n";
        }
        else{
            new_edge = _addEdge(u, v, u->m->getWidth() + _powerplanWidth);
            cout << "Move vertical edge (" << e->u->m->getName() << ", " << e->v->m->getName() << ") to horizontal edge (" << e->u->m->getName() << ", " << e->v->m->getName() << ")\n";
        }
    }
    else{
        assert(u == u->m->_vertexXDPVertical);
        if(u->m->getY() > v->m->getY()){
            new_edge = _addEdge(v, u, v->m->getHeight() + _minChannelSpacing);
            cout << "Move horizontal edge (" << e->u->m->getName() << ", " << e->v->m->getName() << ") to vertical edge (" << e->v->m->getName() << ", " << e->u->m->getName() << ")\n";
        }
        else{
            new_edge = _addEdge(u, v, u->m->getHeight() + _minChannelSpacing);
            cout << "Move horizontal edge (" << e->u->m->getName() << ", " << e->v->m->getName() << ") to vertical edge (" << e->u->m->getName() << ", " << e->v->m->getName() << ")\n";
        }
    }
    _delEdge(e);
    return new_edge;
}

Macro* Legalizer::_randomMacro() const{
    int idx;
    do{
        idx = rand() % _numMacros;
    }while(_macroList[idx]->isFixed());
    return _macroList[idx];
}

EdgeXDP* Legalizer::_randomEdge() const{
    int idx_1, idx_2;
    bool horizontal = rand() % 2;
    EdgeXDP* e;
    VertexXDP* u;
    VertexXDP* v;
    do{
        idx_1 = rand() % _numMacros;
    }while(_macroList[idx_1]->isFixed());
    size_t count = 0;
    if(horizontal){
        u = _macroList[idx_1]->_vertexXDPHorizontal;
        do{
            idx_2 = rand() % u->outList.size();
            v = u->outList[idx_2]->v;
            ++count;
            if(count == _numMacros){
                u = _macroList[idx_1]->_vertexXDPVertical;
            }
            if(count == _numMacros * 2){
                e = _randomEdge();
                return e;
            }
        }while(v->m->isFixed());
    }
    else{
        u = _macroList[idx_1]->_vertexXDPVertical;
        do{
            idx_2 = rand() % u->outList.size();
            v = u->outList[idx_2]->v;
            ++count;
            if(count == _numMacros){
                u = _macroList[idx_1]->_vertexXDPHorizontal;
            }
            if(count == _numMacros * 2){
                e = _randomEdge();
                return e;
            }
        }while(v->m->isFixed());
    }
    e = u->outList[idx_2];
    return e;
}

void Legalizer::_sortMacros(){
    _sortMacroList.resize(0);
    for(size_t i=0; i<_numMacros; ++i){
        _sortMacroList.push_back(_macroList[i]);
    }
    std::sort(_sortMacroList.begin(), _sortMacroList.end(), compare_macro());
}

void Legalizer::writeOutput(fstream& out){
    out << "VERSION " << _version << " ;\n";
    out << "DESIGN " << _designName << " ;\n";
    out << "UNITS DISTANCE MICRONS " << _unit << " ;\n";
    out << "\n";
    out << "DIEAREA ";
    for(unsigned i=0; i<_chipBoundary.size(); ++i){
        out << "( " << _chipBoundary[i].first << " " << _chipBoundary[i].second << " ) ";
    }
    out << ";\n";
    out << "\n";
    out << "COMPONENTS " << _numMacros << " ;\n";
    for(unsigned i=0; i<_numMacros; ++i){
        out << "   - " << _macroList[i]->getName() << " " << _macroList[i]->getType() << " \n";
        out << "      + ";
        if(_macroList[i]->isFixed()) out << "FIXED ";
        else out << "PLACED ";
        out << "( " << _macroList[i]->getBestPos().first << " " << _macroList[i]->getBestPos().second << " ) N ;\n";
    }
    out << "END COMPONENTS\n";
    out << "\n\nEND DESIGN\n\n\n";
}

void Legalizer::plotResult(const string outfilename, bool isPrompt){
    ofstream outfile( outfilename.c_str() , ios::out );
    outfile << " " << endl;
    outfile << "set title \"" << _designName << "\"" << endl;
    outfile << "set size ratio 1" << endl;
    outfile << "set nokey" << endl << endl;
    outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << endl << endl;
    outfile << "# bounding box" << endl;
    if(_chipBoundary.size() == 2){
        _plotBoxPLT( outfile, _chipBoundary[0].first, _chipBoundary[0].second, _chipBoundary[1].first, _chipBoundary[1].second );
    }
    else{
        for(size_t i=0; i<_chipBoundary.size(); ++i){
            outfile << _chipBoundary[i].first << ", " << _chipBoundary[i].second << "\n";
        }
    }
    outfile << _chipBoundary[0].first << ", " << _chipBoundary[0].second << "\n\n";
    outfile << "EOF" << endl;
    outfile << "# modules" << endl << "0.00, 0.00" << endl << endl;
    for( size_t i = 0; i < _numMacros; ++i ){
        Macro* macro = _macroList[i];
        _plotBoxPLT(outfile, macro->getBestPos().first, macro->getBestPos().second, macro->getBestPos().first + macro->getWidth(), macro->getBestPos().second + macro->getHeight());
    }
    outfile << "EOF" << endl;
    outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    if( isPrompt ){
        char cmd[ 200 ];
        sprintf( cmd, "gnuplot %s", outfilename.c_str() );
        if( !system( cmd ) ) { cout << "Fail to execute: \"" << cmd << "\"." << endl; }
    }
}

void Legalizer::_plotBoxPLT( ofstream& stream, int x1, int y1, int x2, int y2 ){
    stream << x1 << ", " << y1 << endl << x2 << ", " << y1 << endl
           << x2 << ", " << y2 << endl << x1 << ", " << y2 << endl
           << x1 << ", " << y1 << endl << endl;
}
