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

class AMNLocation2strains;
class AMNGraph2strains;

class AMNLocation2strains{
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
  int E1_b;
  int E2_b;
  int I_b;
  int Ia_b;

  vector<int> visitorsS;
  vector<int> visitorsE1;
  vector<int> visitorsE2;
  vector<int> visitorsI;
  vector<int> visitorsIa;
  vector<int> visitorsR;
  vector<int> visitorsE1_b;
  vector<int> visitorsE2_b;
  vector<int> visitorsI_b;
  vector<int> visitorsIa_b;

  AMNLocation2strains(string name_, int Shome);

  void print();
  void seir_step(float beta, float a1, float a2, float gamma, float relativeInfectiousnessB, float presymptomaticRelativeInfectiousness, float asymptomaticProb, float asymptomaticRelativeInfectiousness, int &de2, int &dea, int &de2_b, int &dea_b);
};



class AMNGraph2strains{
public:
  vector<AMNLocation2strains> locations;

  AMNGraph2strains();
  void add_node(string name, int Shome);
  void copy_graph(AMNGraph2strains G);
  void print();
  void count_everyone(string msg);
};
