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

class Location;
class Link;
class Graph;

class Location{
	private:

	public:
		string name;
		int S;
		int E;
		int I;
		int Ia;
		int R;
		int N;
		vector<Link*> out_links;
		vector<Link*> in_links;


	Location(string name_, int Shome);
	void add_inlink(Link *link);
	void add_outlink(Link *link);
	void print();
	void seir_step_day(float beta, float a, float gamma, float asymptomaticProb, float asymptomaticRelativeInfectiousness, int &de2);
	void seir_step_night(float beta, float a, float gamma, float asymptomaticProb, float asymptomaticRelativeInfectiousness, int &de2);
};

class Link{

	public:
		int S;
		int E;
		int I;
		int Ia;
		int R;
		Location *from;
		int from_index;
		Location *to;
		int to_index;
	Link(Location *from_, Location *to_, int S_, int E_, int I_, int Ia_, int R_);
	void print();
};


class Graph{
	public:
		vector<Location> locations;
		vector<Link> edges;
	Graph();
	void add_node(string name, int Shome);
	void add_edge(string name1, string name2, int S, int E, int I, int Ia, int R);
	void add_edge_index(int i1, int i2, int S, int E, int I, int Ia, int R);
	void inform_locations_of_edges();
	void copy_graph(Graph G);
	void print();
};
