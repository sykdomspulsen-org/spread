#include <Rcpp.h>
#include "solveig_shared.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include "commuter.h"

using namespace Rcpp;
using namespace std;




Location::Location(string name_, int Shome){
  name = name_;
  S = Shome;
  E = 0;
  I = 0;
  Ia = 0;
  R = 0;
  N = 0;
}

void Location::add_inlink(Link *link){
  // This function adds a pointer to a newly created link that points to this location
  in_links.push_back(link);
}

void Location::add_outlink(Link  *link){
  // This function adds a pointer to a newly created link that points from this location
  out_links.push_back(link);

}

void Location::print(){
  // This is a debug function
  Rcout << "Location " << name << " with S= "<< S << ", E=" << E <<", I=" << I << ", Ia=" << Ia << ", R= " << R << ". " << endl;
  Rcout << "out_links: "<< out_links.size() <<endl;
  for(unsigned int i = 0; i < out_links.size(); ++i){
    Rcout << "Trying to print from out_link " << i << endl;
    out_links[i]->print();
    Rcout.flush();
  }
  Rcout << "in_links: "<< in_links.size()<< endl;
  for(unsigned int i = 0; i < in_links.size(); ++i){
    Rcout << "Trying to print from in_link " << i << endl;
    in_links[i]->print();
    Rcout.flush();
  }

}

void Location::seir_step_day(
    // Function to run a day time step of the model.
    // The commuters are all sent to their work locations.
    float beta,
    float a,
    float gamma,
    float asymptomaticProb,
    float asymptomaticRelativeInfectiousness,
    int &de2){ // Return symptomatic incidence
  int S_tmp = S;
  int E_tmp = E;
  int Ia_tmp = Ia;
  int I_tmp = I;
  int R_tmp = R;
  double pop_tmp = S + E + Ia + I + R;
  int num = in_links.size();

  //Vectors with the number of people in the respective compartment,
  //on each link, and in the home population
  //used to distribute the transitions between compartments
  //between the commuters on the different links and the home population
  //The in-edges are on the first num elements, the home population last
  int *S_probs = new int[num+1];
  int *E_probs = new int[num+1];
  int *I_probs = new int[num+1];
  int *Ia_probs = new int[num+1];
  int *R_probs = new int[num+1];
  for(unsigned int i = 0; i < in_links.size(); ++i){
    S_tmp += in_links[i]->S;
    E_tmp += in_links[i]->E;
    Ia_tmp += in_links[i]->Ia;
    I_tmp += in_links[i]->I;
    R_tmp += in_links[i]->R;
    pop_tmp += in_links[i]->S + in_links[i]->E + in_links[i]->Ia + in_links[i]->I + in_links[i]->R;

    S_probs[i] = in_links[i]->S;
    I_probs[i] = in_links[i]->I;
    Ia_probs[i] = in_links[i]->Ia;
    E_probs[i] = in_links[i]->E;
    R_probs[i] = in_links[i]->R;
  }

  S_probs[num] = S;
  I_probs[num] = I;
  Ia_probs[num] = Ia;
  E_probs[num] = E;
  R_probs[num] = R;
  int ds; int de1; int dia; int di;

  // Run the SEIR step
  seir_sim(ds, de1, de2, dia, di, S_tmp, E_tmp, Ia_tmp, I_tmp, beta, a, gamma, asymptomaticProb, asymptomaticRelativeInfectiousness, pop_tmp, 12.0/24.0);

  //Distribute the transitions
  double *probs = new double[num+1];
  double *probs_cum = new double[num+1];
  double randomnumber;
  int index = -1;
  for(int i = 0; i < dia; ++i){
    for(int k = 0; k < num+1; ++k) probs[k] =Ia_probs[k]*1.0/ Ia_tmp; // Oneline for-loop
    randomnumber = R::runif(0.0,1.0);
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = R::runif(0.0,1.0);
    }
    cumulative_sum(&probs_cum, &probs, num+1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      Ia -= 1;
      R += 1;
    }
    else{
      in_links[index]->Ia -= 1;
      in_links[index]->R += 1;
    }

    Ia_probs[index] -= 1;
    Ia_tmp -= 1;
  }

  for(int i = 0; i < di; ++i){
    for(int k = 0; k < num+1; ++k) probs[k] =I_probs[k]*1.0/ I_tmp; // Oneline for-loop
    randomnumber = R::runif(0.0,1.0);
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = R::runif(0.0,1.0);
    }
    cumulative_sum(&probs_cum, &probs, num+1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      I -= 1;
      R += 1;
    }
    else {
      in_links[index]->I -= 1;
      in_links[index]->R += 1;
    }
    I_probs[index] -= 1;
    I_tmp -= 1;
  }


  for(int i = 0; i < de1; ++i){
    for(int k = 0; k < num+1; ++k) probs[k] =E_probs[k]*1.0/ E_tmp; // Oneline for-loop
    randomnumber = R::runif(0.0,1.0);
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = R::runif(0.0,1.0);
    }
    cumulative_sum(&probs_cum, &probs, num+1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      E -= 1;
      Ia += 1;
    }
    else{
      in_links[index]->E -= 1;
      in_links[index]->Ia += 1;
    }

    E_probs[index] -= 1;
    E_tmp -= 1;

  }

  for(int i = 0; i < de2; ++i){
    for(int k = 0; k < num+1; ++k) probs[k] =E_probs[k]*1.0/ E_tmp; // Oneline for-loop
    randomnumber = R::runif(0.0,1.0);
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = R::runif(0.0,1.0);
    }
    cumulative_sum(&probs_cum, &probs, num+1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      E -= 1;
      I += 1;
    }
    else{
      in_links[index]->E -= 1;
      in_links[index]->I += 1;
    }
    E_probs[index] -= 1;
    E_tmp -= 1;

  }


  for(int i = 0; i < ds; ++i){
    for(int k = 0; k < num+1; ++k) probs[k] =S_probs[k]*1.0/S_tmp; // Oneline for-loop
    randomnumber = R::runif(0.0,1.0);
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = R::runif(0.0,1.0);
    }
    cumulative_sum(&probs_cum, &probs, num+1);
    for(int h = 0; h < num+1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      S -= 1;
      E += 1;
    }
    else{
      in_links[index]->S -= 1;
      in_links[index]->E += 1;
    }
    S_probs[index] -= 1;
    S_tmp -= 1;

  }


  delete[] probs;
  delete[] probs_cum;

  delete[] S_probs;
  delete[] E_probs;
  delete[] I_probs;
  delete[] Ia_probs;
  delete[] R_probs;
}

void Location::seir_step_night(
    // Function to run a night time step of the model.
    // The commuters are all sent to their home locations.
    float beta,
    float a,
    float gamma,
    float asymptomaticProb,
    float asymptomaticRelativeInfectiousness,
    int &de2){ //Return symptomatic incidence
  int S_tmp = S;
  int E_tmp = E;
  int Ia_tmp = Ia;
  int I_tmp = I;
  int R_tmp = R;
  double pop_tmp = S + E + Ia + I + R;
  unsigned int num = out_links.size();
  //Vectors with the number of people in the respective compartment,
  //on each link, and in the home population
  //used to distribute the transitions between compartments
  //between the commuters on the different links and the home population
  //The out-edges are on the first num elements, the home population last
  int *S_probs = new int[num+1];
  int *E_probs = new int[num+1];
  int *I_probs = new int[num+1];
  int *Ia_probs = new int[num+1];
  int *R_probs = new int[num+1];
  for(unsigned int i = 0; i < out_links.size(); ++i){
    S_tmp += out_links[i]->S;
    E_tmp += out_links[i]->E;
    Ia_tmp += out_links[i]->Ia;
    I_tmp += out_links[i]->I;
    R_tmp += out_links[i]->R;
    pop_tmp += out_links[i]->S + out_links[i]->E + out_links[i]->Ia + out_links[i]->I + out_links[i]->R;

    S_probs[i] = out_links[i]->S;
    I_probs[i] = out_links[i]->I;
    Ia_probs[i] = out_links[i]->Ia;
    E_probs[i] = out_links[i]->E;
    R_probs[i] = out_links[i]->R;
  }
  S_probs[num] = S;
  I_probs[num] = I;
  Ia_probs[num] = Ia;
  E_probs[num] = E;
  R_probs[num] = R;
  int ds; int de1; int dia; int di;

  // Run the SEIR step
  seir_sim(ds, de1, de2, dia, di, S_tmp, E_tmp, Ia_tmp, I_tmp, beta, a, gamma, asymptomaticProb, asymptomaticRelativeInfectiousness, pop_tmp, 12.0/24.0);

  //Distribute the transitions
  double *probs = new double[num+1];
  double *probs_cum = new double[num+1];
  double randomnumber;
  unsigned int index = 0;
  for(int i = 0; i < dia; ++i){
    for(unsigned int k = 0; k < num+1; ++k) probs[k] =Ia_probs[k]*1.0/ Ia_tmp; // Oneline for-loop
    randomnumber = R::runif(0.0,1.0);
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = R::runif(0.0,1.0);
    }
    cumulative_sum(&probs_cum, &probs, num+1);
    for(unsigned int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      Ia -= 1;
      R += 1;
    }
    else if(index < out_links.size()){
      out_links[index]->Ia -= 1;
      out_links[index]->R += 1;
    }

    Ia_probs[index] -= 1;
    Ia_tmp -= 1;

  }

  for(int i = 0; i < di; ++i){
    for(unsigned int k = 0; k < num+1; ++k) probs[k] =I_probs[k]*1.0/ I_tmp; // Oneline for-loop
    randomnumber = R::runif(0.0,1.0);
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = R::runif(0.0,1.0);
    }
    cumulative_sum(&probs_cum, &probs, num+1);
    for(unsigned int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      I -= 1;
      R += 1;
    }
    else if(index < out_links.size()){
      out_links[index]->I -= 1;
      out_links[index]->R += 1;
    }


    I_probs[index] -= 1;
    I_tmp -= 1;

  }

  for(int i = 0; i < de1; ++i){
    for(unsigned int k = 0; k < num+1; ++k) probs[k] =E_probs[k]*1.0/ E_tmp; // Oneline for-loop
    randomnumber = R::runif(0.0,1.0);
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = R::runif(0.0,1.0);
    }
    cumulative_sum(&probs_cum, &probs, num+1);
    for(unsigned int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      E -= 1;
      Ia += 1;
    }
    else if(index < out_links.size()){
      out_links[index]->E -= 1;
      out_links[index]->Ia += 1;
    }

    E_probs[index] -= 1;
    E_tmp -= 1;

  }

  for(int i = 0; i < de2; ++i){
    for(unsigned int k = 0; k < num+1; ++k) probs[k] =E_probs[k]*1.0/ E_tmp; // Oneline for-loop
    randomnumber = R::runif(0.0,1.0);
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = R::runif(0.0,1.0);
    }
    cumulative_sum(&probs_cum, &probs, num+1);
    for(unsigned int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      E -= 1;
      I += 1;
    }
    else if(index < out_links.size()){
      out_links[index]->E -= 1;
      out_links[index]->I += 1;
    }
    E_probs[index] -= 1;
    E_tmp -= 1;

  }

  for(int i = 0; i < ds; ++i){
    for(unsigned int k = 0; k < num+1; ++k) probs[k] =S_probs[k]*1.0/ S_tmp; // Oneline for-loop
    randomnumber = R::runif(0.0,1.0);
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = R::runif(0.0,1.0);
    }
    cumulative_sum(&probs_cum, &probs, num+1);
    for(unsigned int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }

    if(index == num){
      // Add to home pop
      S -= 1;
      E += 1;
    }
    else if(index < out_links.size()){
      out_links[index]->S -= 1;
      out_links[index]->E += 1;
    }

    S_probs[index] -= 1;
    S_tmp -= 1;

  }

  delete[] probs;
  delete[] probs_cum;

  delete[] S_probs;
  delete[] E_probs;
  delete[] I_probs;
  delete[] Ia_probs;
  delete[] R_probs;
}


Link::Link(Location *from_, Location *to_, int S_, int E_, int I_, int Ia_, int R_){
  from = from_;
  to = to_;
  S  = S_;
  E  = E_;
  I  = I_;
  Ia = Ia_;
  R  = R_;
}

void Link::print(){
  Rcout << "Link from " << from->name << ", to " << to->name << ", with S= "<< S << ", E=" << E <<", I=" << I << ", Ia=" << Ia << ", R= " << R << endl;
}


Graph::Graph(){
}

void Graph::add_node(string name, int Shome){
  Location newlocation(name, Shome);
  locations.push_back(newlocation);
}

void Graph::add_edge(string name1, string name2, int S, int E, int I, int Ia, int R){
  /// Find locations with names name1 and name2
  int name1_index = -1;
  int name2_index = -1;
  for (unsigned int i = 0; i < locations.size(); ++i){
    if(locations[i].name.compare(name1) == 0){
      name1_index = i;
    }
    if(locations[i].name.compare(name2) == 0){
      name2_index = i;
    }
  }
  if( ((name1_index == -1) || (name2_index == -1)) || (name1_index == name2_index)){
    Rcout << "Error in add edge with input: " << name1 << ", " << name2 << ", " << S << ", " << E << ", " << I << ", " << Ia << ", " << R << endl;
    Rcout << "name1_index = " << name1_index << ", name2_index = " << name2_index << endl;
    Rcpp::stop("Error!");
  }

  Link newlink(&(locations[name1_index]), &(locations[name2_index]), S, E, I, Ia, R);
  newlink.from_index = name1_index;
  newlink.to_index = name2_index;
  edges.push_back(newlink);
}


void Graph::add_edge_index(int i1, int i2, int S, int E, int I, int Ia, int R){
  Link newlink(&(locations[i1]), &(locations[i2]), S, E, I, Ia, R);
  edges.push_back(newlink);
}

void Graph::inform_locations_of_edges(){
  for (vector<Link>::iterator it = edges.begin() ; it != edges.end(); ++it){
    (*it).from->add_outlink(&(*it));
    (*it).to->add_inlink(&(*it));
  }
}


void Graph::print(){
  Rcout << endl << "Printing graph: " << endl;
  Rcout << "Links: " << endl;
  for(unsigned int i = 0; i < edges.size(); ++i){
    edges[i].print();
  }
  Rcout << "Locations: " << endl;
  for(unsigned int i = 0; i < locations.size(); ++i){
    locations[i].print();
    Rcout.flush();
  }

  Rcout << endl << endl;;
}

void Graph::copy_graph(Graph G){
  for(vector<Location>::iterator it = G.locations.begin() ; it != G.locations.end(); ++it){
    add_node(it->name, it->S);
  }
  for(vector<Link>::iterator it = G.edges.begin() ; it != G.edges.end(); ++it){
    add_edge_index(it->from_index, it->to_index, it->S, it->E, it->I, it->Ia, it->R);
  }
  inform_locations_of_edges();
}

// float beta;  // infection parameter, 0.6
// float a;  // 1/latent period, 1/1.9
// float gamma; // 1/infectious period, 1/3
// float asymptomaticProb; // Proportion/probability of asymptomatic given infectious
// float asymptomaticRelativeInfectiousness; // Relative infectiousness of asymptomatic infectious
// int N = 1; // Number of repetitions
// int M; // Number of days

//' commuter
//' @param seiiar_home Data frame
//' @param seiiar_commuters Data frame
//' @param beta Float, infection parameter, 0.6
//' @param a Float, 1/latent period, 1/1.9
//' @param gamma Float, 1/infectious period, 1/3
//' @param asymptomaticProb Float, Proportion/probability of asymptomatic given infectious
//' @param asymptomaticRelativeInfectiousness Float, Relative infectiousness of asymptomatic infectious
//' @param N Int = 1 int, Number of repetitions
//' @param M Int, Number of days
//' @param verbose Bool
//' @export
// [[Rcpp::export]]
DataFrame commuter_cpp(
    DataFrame seiiar_home,
    DataFrame seiiar_commuters,
    float beta,
    float a,
    float gamma,
    float asymptomaticProb,
    float asymptomaticRelativeInfectiousness,
    int N=1,
    int M=56,
    bool verbose=1) {

  int n=0; //Number of locations

  Graph G;
  // Read node file with the names and population sizes without commuters
  StringVector names = seiiar_home[0] ;
  IntegerVector home_S = seiiar_home[1] ;
  IntegerVector home_E = seiiar_home[2] ;
  IntegerVector home_I = seiiar_home[3] ;
  IntegerVector home_Ia = seiiar_home[4] ;
  IntegerVector home_R = seiiar_home[5] ;

  for (int i = 0; i < seiiar_home.rows(); i++) {
    string name = std::string(names[i]);
    G.add_node(name, 0);
    n+= 1;
  }

  // Read edges file with the number of commuters
  // Should not include edges with zeros for computational efficiency
  int safecount = 0;
  int n_edges = 0;

  StringVector names_from = seiiar_commuters[0] ;
  StringVector names_to = seiiar_commuters[1] ;
  IntegerVector commuters_S = seiiar_commuters[2] ;
  IntegerVector commuters_E = seiiar_commuters[3] ;
  IntegerVector commuters_I = seiiar_commuters[4] ;
  IntegerVector commuters_Ia = seiiar_commuters[5] ;
  IntegerVector commuters_R = seiiar_commuters[6] ;

  if(verbose) Rcout << "Starting to add edges, printing every 1000 edge" << endl;
  for (int i = 0; i < seiiar_commuters.rows(); i++) {
    string name_from = std::string(names_from[i]);
    string name_to = std::string(names_to[i]);
    int c_S = commuters_S[i];
    int c_E = commuters_E[i];
    int c_I = commuters_I[i];
    int c_Ia = commuters_Ia[i];
    int c_R = commuters_R[i];

    if (c_S != 0 || c_E !=0 || c_I != 0 || c_Ia != 0 || c_R != 0){
      safecount += 1;
      // name_from, name_to, S, E, I, Ia, R
      G.add_edge(name_from, name_to, c_S, c_E, c_I, c_Ia, c_R);

      if(verbose & (safecount % 1000 == 0)){
        Rcout << safecount << " ";
        Rcout.flush();
      }
      n_edges += 1;
    }
  }
  if(verbose) Rcout << "Found " << n_edges << " edges" << endl;


  G.inform_locations_of_edges();
  // Storage for statistics

  //The state in each location at each time point
  int ***values = new int**[n];

  // Peak dates in each location
  int **peak_date = new int*[n];

  // Peak number infected in each location
  int **peak_val = new int*[n];

  //Initial dates in each location
  //defined as first day where the symptomatic prevalence has been
  //more than 1% for 7 consecutive days
  int **start_date = new int*[n];

  //Dummy vector, to find peak and initial dates
  int **I_this = new int*[n];

  //Final number infected in each location
  int **final_size = new int*[n];

  //Save the prevalence curves (infectious + infectious asymptomatic) for each location
  // For each run, in order to make confidence curves over the N simulations
  int **bonds = new int*[N];

  /// Values has first index for position, second for time and third for value.
  /// Third index 0=S, 1=E, 2=I, 3=Ia, 4=R; 5= symptomatic incidence;
  for(int i = 0; i < n; ++i){
    values[i] = new int*[M*2]; // 2*M, stored for both day and night time.
    peak_date[i] = new int[N];
    peak_val[i] = new int[N];
    start_date[i] = new int[N];
    I_this[i] = new int[M*2];
    final_size[i] = new int[N];
    for(int k = 0; k < M*2; ++k){
      values[i][k] = new int[6];
      values[i][k][0] = 0;
      values[i][k][1] = 0;
      values[i][k][2] = 0;
      values[i][k][3] = 0;
      values[i][k][4] = 0;
      values[i][k][5] = 0;
    }
  }

  for (int i = 0; i < N; ++i){
    bonds[i] = new int[2*M]; //2*M, stored for both day and night time
    for (int j = 0; j < 2*M; ++j){
      bonds[i][j] = 0;
    }
  }

  if(verbose) Rcout << "Running " << N << " simulations of " << M << " days" << endl << endl;
  Progress p(N*M, verbose);
  for(int i_sim = 0; i_sim < N; ++i_sim){

    for (int i = 0; i < n; ++i){
      final_size[i][i_sim] = 0;
      for (int k = 0; k < M*2; ++k){
        I_this[i][k] = 0;
      }
    }
    Graph G_current;
    G_current.copy_graph(G);

    // Seed the epidemic
    for(int i = 0; i < n; ++i){
      G_current.locations[i].S = home_S[i];
      G_current.locations[i].E = home_E[i];
      G_current.locations[i].I = home_I[i];
      G_current.locations[i].Ia = home_Ia[i];
      G_current.locations[i].R = home_R[i];
    }
    //Rcout << "Starting simulation " << i_sim+1 << "/" << N << endl;

    for(int i_day = 0; i_day < M; ++i_day){
      //Rcout << "Day " << i_day+1 << "/" << M << endl;
      p.increment(); // update progress
      if (Progress::check_abort() )
        return -1.0;

      for (int i = 0; i < n; ++i){
        int de2 = 0;
        // Let commuters mix at their work location during day time
        G_current.locations[i].seir_step_day(beta, a, gamma, asymptomaticProb, asymptomaticRelativeInfectiousness, de2);
        values[i][2*i_day][0] += G_current.locations[i].S;
        values[i][2*i_day][1] += G_current.locations[i].E;
        values[i][2*i_day][2] += G_current.locations[i].I;
        values[i][2*i_day][3] += G_current.locations[i].Ia;
        values[i][2*i_day][4] += G_current.locations[i].R;
        I_this[i][2*i_day] += G_current.locations[i].I;
        int num = G_current.locations[i].out_links.size();
        values[i][2*i_day][5] += de2;
        bonds[i_sim][2*i_day] += G_current.locations[i].I + G_current.locations[i].Ia;
        for(int j = 0; j < num; ++j){
          values[i][2*i_day][0] += G_current.locations[i].out_links[j]->S;
          values[i][2*i_day][1] += G_current.locations[i].out_links[j]->E;
          values[i][2*i_day][2] += G_current.locations[i].out_links[j]->I;
          values[i][2*i_day][3] += G_current.locations[i].out_links[j]->Ia;
          values[i][2*i_day][4] += G_current.locations[i].out_links[j]->R;
          I_this[i][2*i_day] += G_current.locations[i].out_links[j] -> I;
          bonds[i_sim][2*i_day] += G_current.locations[i].out_links[j]->I + G_current.locations[i].out_links[j]->Ia;
        }
      }
      for (int i = 0; i < n; ++i){
        int de2 = 0;
        // Let commuters mix in their home location during night time
        G_current.locations[i].seir_step_night(beta, a, gamma, asymptomaticProb, asymptomaticRelativeInfectiousness, de2);
        if(i_day == (M-1)){
          final_size[i][i_sim] += G_current.locations[i].R;
        }
        values[i][2*i_day+1][0] += G_current.locations[i].S;
        values[i][2*i_day+1][1] += G_current.locations[i].E;
        values[i][2*i_day+1][2] += G_current.locations[i].I;
        values[i][2*i_day+1][3] += G_current.locations[i].Ia;
        values[i][2*i_day+1][4] += G_current.locations[i].R;
        I_this[i][2*i_day+1] += G_current.locations[i].I;
        values[i][2*i_day+1][5] += de2;
        bonds[i_sim][2*i_day+1] += G_current.locations[i].I + G_current.locations[i].Ia;
        int num = G_current.locations[i].out_links.size();
        for(int j = 0; j < num; ++j){
          if (i_day == (M-1)){
            final_size[i][i_sim] += G_current.locations[i].out_links[j]->R;
          }
          values[i][2*i_day+1][0] += G_current.locations[i].out_links[j]->S;
          values[i][2*i_day+1][1] += G_current.locations[i].out_links[j]->E;
          values[i][2*i_day+1][2] += G_current.locations[i].out_links[j]->I;
          values[i][2*i_day+1][3] += G_current.locations[i].out_links[j]->Ia;
          values[i][2*i_day+1][4] += G_current.locations[i].out_links[j]->R;
          I_this[i][2*i_day+1] += G_current.locations[i].out_links[j]->I;
          bonds[i_sim][2*i_day+1] += G_current.locations[i].out_links[j]->I + G_current.locations[i].out_links[j]->Ia;
        }
      }
    }
    // Initial dates and peak dates.
    float baseline;
    int pop = 0;
    int pday = 0;
    int pmax = 0;
    int count;
    int startday = 0;
    for(int i = 0; i < n; i++){
      pday = 0;
      pmax = 0;
      startday = 0;
      count = 0;
      pop = G_current.locations[i].S +  G_current.locations[i].E +  G_current.locations[i].I +  G_current.locations[i].Ia +  G_current.locations[i].R;
      int num = G.locations[i].out_links.size();
      for (int j = 0; j<num; ++j){
        pop += G.locations[i].out_links[j]->S + G.locations[i].out_links[j]->E + G.locations[i].out_links[j]->I + G.locations[i].out_links[j]->Ia + G.locations[i].out_links[j]->R;
      }
      if(pop < 100){
        baseline = 1;
      }
      else{
        baseline = 0.01*pop;
      }
      for(int i_day = 0; i_day < M; ++i_day){
        if(I_this[i][2*i_day] > pmax){
          pday = i_day;
          pmax = I_this[i][2*i_day];
        }
        if(I_this[i][2*i_day] > baseline){
          count += 1;
          if (count > 6){
            startday = i_day;
            count = -20000;
          }
        }
      }
      peak_date[i][i_sim] = pday;
      peak_val[i][i_sim] = pmax;
      start_date[i][i_sim] = startday;
    }

    //Rcout << "Finished simulation " << i_sim+1 << "/" << N << endl;
  }

  if(verbose) Rcout << endl << "Finished all simulations" << endl;

  StringVector res_names(n*2*M);
  IntegerVector res_week(n*2*M);
  IntegerVector res_day(n*2*M);
  LogicalVector res_6pm(n*2*M);
  IntegerVector res_S(n*2*M);
  IntegerVector res_E(n*2*M);
  IntegerVector res_I(n*2*M);
  IntegerVector res_Ia(n*2*M);
  IntegerVector res_R(n*2*M);
  IntegerVector res_INCIDENCE(n*2*M);

  int index=0;
  for(int i=0; i < n; ++i){
    for (int k = 0; k<2*M; ++k){

      res_names[index] = names[i];
      res_week[index] = k/14+1;
      res_day[index] = k/2+1;
      res_6pm[index] = k%2;
      res_S[index] = values[i][k][0]*1.0/N;
      res_E[index] = values[i][k][1]*1.0/N;
      res_I[index] = values[i][k][2]*1.0/N;
      res_Ia[index] = values[i][k][3]*1.0/N;
      res_R[index] = values[i][k][4]*1.0/N;
      res_INCIDENCE[index] = values[i][k][5]*1.0/N;

      index++;
    }
  }

  // return a new data frame
  DataFrame df = DataFrame::create(
    _["location_code"]= res_names,
    _["week"]=res_week,
    _["day"]=res_day,
    _["is_6pm"]=res_6pm,
    _["S"]= res_S,
    _["E"]= res_E,
    _["I"]= res_I,
    _["Ia"]= res_Ia,
    _["R"]= res_R,
    _["incidence"]= res_INCIDENCE
  );

  df.attr("class") = Rcpp::CharacterVector::create("data.table", "data.frame");

  return(df);
}

/*** R
x <- spread:::commuter_convert_seiiar(
  seiiar=spread::norway_seiiar_oslo_2017_b2020,
  commuters=spread::norway_commuters_2017_b2020
)

d <- commuter_cpp(
  seiiar_home=x[["seiiar_home"]],
  seiiar_commuters=x[["seiiar_commuters"]],
  beta=1,
  a=1,
  gamma=1,
  asymptomaticProb = 1,
  asymptomaticRelativeInfectiousness = 1,
  N=98,
  M=99)
d
*/

