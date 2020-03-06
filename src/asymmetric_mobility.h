#include <iostream>
#include <fstream>
#include <cstdlib>
#include <list>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <stdlib.h>
#include <time.h>

using namespace std;

class AMLocation;
class AMGraph;

class AMLocation{
private:

public:
  string name;
  int S;
  int E;
  int I;
  int Ia;
  int R;
  int N;

  vector<int> visitorsS;
  vector<int> visitorsE;
  vector<int> visitorsI;
  vector<int> visitorsIa;
  vector<int> visitorsR;

  AMLocation(string name_, int Shome);

  void print();
  void seir_step(float beta, float a, float gamma, float asymptomaticProb, float asymptomaticRelativeInfectiousness, int &de2);
};



class AMGraph{
public:
  vector<AMLocation> locations;

  AMGraph();
  void add_node(string name, int Shome);
  void copy_graph(AMGraph G);
  void print();
  void count_everyone(string msg);
};
