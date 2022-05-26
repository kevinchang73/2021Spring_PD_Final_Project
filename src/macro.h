#ifndef MACRO_H
#define MACRO_H

#include <vector>
#include <string>
using namespace std;

class VertexXDP;

class Macro
{
public:
    Macro() {}
    Macro(string name, string type, int width, int height, int x, int y, bool fixed, bool psuedo = false):
        _name(name), _type(type), _psuedo(psuedo), _fixed(fixed), _w(width), _h(height), 
        _currentX(x), _currentY(y), _initX(x), _initY(y), _bestX(x), _bestY(y)
    {}
    ~Macro() {}

    void setCurrentPos(int x, int y) { 
        _currentX = x; 
        _currentY = y; 
    }
    void storeBestPos() { 
        _bestX = _currentX; 
        _bestY = _currentY; 
    }

    string getName(){
        return _name;
    }
    string getType(){
        return _type;
    }
    bool isFixed(){
        return _fixed;
    }
    pair<int, int> getBestPos(){
        return pair<int, int>(_bestX, _bestY);
    }
    int getBestX() { return _bestX; }
    int getBestY() { return _bestY; }
    int getX() { return _currentX; }
    int getY() { return _currentY; }
    int getInitX() { return _initX; }
    int getInitY() { return _initY; }
    int getWidth() { return _w; }
    int getHeight() { return _h; }
    int getSize() { return _w * _h; }

    void reportInfo() {
        cout << _name << " " << _w << " " << _h << " " << _currentX << " " << _currentY << "\n";
    }

private:
    string _name;
    string _type;
    bool _psuedo;
    bool _fixed;
    int _w;
    int _h;
    //current position (left-bottom)
    int _currentX;
    int _currentY;
    //initial position (left-bottom)
    int _initX;
    int _initY;
    //best position (left-bottom)
    int _bestX;
    int _bestY;
public:
    //constraint graph
    VertexXDP* _vertexXDPVertical;
    VertexXDP* _vertexXDPHorizontal;
};

class EdgeXDP
{
public:
    EdgeXDP() {}
    EdgeXDP(VertexXDP* u, VertexXDP* v, int w):
        u(u), v(v), weight(w)
    {}

    VertexXDP* u;
    VertexXDP* v;
    int weight;
    int slack;
};

class VertexXDP
{
public:
    VertexXDP() {}
    VertexXDP(Macro* m, bool v, int i):
        m(m), idx(i), isVertical(v)
    {}
    VertexXDP(VertexXDP* v){
        m = v->m;
        isVertical = v->isVertical;
        for(unsigned i=0; i<v->outList.size(); ++i) outList.push_back(v->outList[i]);
        for(unsigned i=0; i<v->inList.size(); ++i) inList.push_back(v->inList[i]);
    }
    int slack() { return (RT - LB); }

    Macro* m;
    int idx;
    bool isVertical; //1:Gv; 0:Gh
    int LB; //left / bottom
    int RT; //right / top
    vector<EdgeXDP*> outList;
    vector<EdgeXDP*> inList;
};

#endif  // MACRO_H
