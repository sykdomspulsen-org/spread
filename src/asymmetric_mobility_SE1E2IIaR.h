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

class AMNLocation;
class AMNGraph;

class AMNLocation{
private:

public:
  string name;
  int S;
  int E1;
  int E2;
  int I;
  int Ia;
  int R;
  int N;

  vector<int> visitorsS;
  vector<int> visitorsE1;
  vector<int> visitorsE2;
  vector<int> visitorsI;
  vector<int> visitorsIa;
  vector<int> visitorsR;

  AMNLocation(string name_, int Shome);

  void print();
  void seir_step(float beta, float a1, float a2, float gamma, float presymptomaticRelativeInfectiousness, float asymptomaticProb, float asymptomaticRelativeInfectiousness, int &de2, int &dea);
};



class AMNGraph{
public:
  vector<AMNLocation> locations;

  AMNGraph();
  void add_node(string name, int Shome);
  void copy_graph(AMNGraph G);
  void print();
  void count_everyone(string msg);
};
