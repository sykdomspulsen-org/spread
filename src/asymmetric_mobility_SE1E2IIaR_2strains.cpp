#include <Rcpp.h>

//~ * Compile with:
//~ g++ -O3 -Wall -I/usr/local/Cellar/gsl/1.16/include -c commuter_mobphone.cpp
//~ * g++ -O3 commuter_mobphone.o -lgsl
//~ ./a.out
//~ * */
#include "asymmetric_mobility_SE1E2IIaR_2strains.h"
using namespace Rcpp;
using namespace std;


void cumulative_sum2strains(double **outarray, double **array, int n){
  /// REMEMBER to delete outarray after use
  (*outarray)[0] = (*array)[0];
  for(int i=1; i < n; ++i){
    (*outarray)[i] = (*outarray)[i-1] + (*array)[i];
  }
}

void se1e2iiar_sim(int &ds, int &de1e2, int &de1ia, int &de2i, int &dia, int &di,
                   int &ds_b, int & de1e2_b, int &de1ia_b, int &de2i_b, int &dia_b, int &di_b,     // Outputs
                   int S, int E1, int E2, int Ia, int I, int E1_b, int E2_b, int Ia_b, int I_b,
                   float beta, float a1, float a2, float gamma,
                   float relativeInfectiousnessB,
                   float presymptomaticRelativeInfectiousness,
                   float asymptomaticProb, float asymptomaticRelativeInfectiousness,
                   int pop, double delta_t){ //Inputs
  ds = 0; de1e2 = 0; de1ia = 0; de2i = 0; dia = 0; di = 0;
  ds_b = 0; de1e2_b = 0; de1ia_b = 0; de2i_b = 0; dia_b = 0; di_b = 0;
  // Function to run one time step of the seir model
  // ds are susceptible going to exposed 1
  // de1e2 are the exposed 1 going to exposed with symptoms, exposed 2
  // de1ia are the exposed 1 going to asymptomatic infectious
  // de2i are the exposed 2 going to symptomatic
  // dia are the infectious asymptomatic going to recovered
  // di are the infectious symptomatic going to recovered
  // with _b are the ones who are infected with strain b.
  int de = 0;
  int de_b = 0;
  if(I != 0 || Ia != 0 || E1 != 0 || E2 != 0 || I_b != 0 || Ia_b != 0 || E1_b != 0 || E2_b != 0){
    if(E1 == 0){
      de = 0;
    }
    else{
      de = R::rbinom(E1, a1 * delta_t);
      if(de != 0){
        de1ia = R::rbinom(de, asymptomaticProb);
        de1e2 = de-de1ia;
      }
    }
    if(E1_b == 0){
      de_b = 0;
    }
    else{
      de_b = R::rbinom(E1_b, a1 * delta_t);
      if(de_b != 0){
        de1ia_b = R::rbinom(de_b, asymptomaticProb);
        de1e2_b = de_b-de1ia_b;
      }
    }
    if(E2 == 0){
      de2i = 0;
    }
    else{
      de2i = R::rbinom(E2, a2 * delta_t);
    }
    if(E2_b == 0){
      de2i_b = 0;
    }
    else{
      de2i_b = R::rbinom(E2_b, a2 * delta_t);
    }

    if(I == 0){
      di = 0;
    }
    else{
      di = R::rbinom(I, gamma*delta_t);
    }
    if(I_b == 0){
      di_b = 0;
    }
    else{
      di_b = R::rbinom(I_b, gamma*delta_t);
    }

    if(Ia == 0){
      dia = 0;
    }
    else{
      dia = R::rbinom(Ia, gamma*delta_t);
    }
    if(Ia_b == 0){
      dia_b = 0;
    }
    else{
      dia_b = R::rbinom(Ia_b, gamma*delta_t);
    }

    double inf_a = beta*delta_t*I/pop + asymptomaticRelativeInfectiousness*beta*delta_t*Ia/pop + presymptomaticRelativeInfectiousness * beta * delta_t * E2 / pop;
    double inf_b = relativeInfectiousnessB * (beta*delta_t*I_b/pop + asymptomaticRelativeInfectiousness*beta*delta_t*Ia_b/pop + presymptomaticRelativeInfectiousness * beta * delta_t * E2_b / pop);
    double p_b = (inf_b)/(inf_a + inf_b);
    ds = R::rbinom(S, inf_a + inf_b);
    ds_b = R::rbinom(ds, p_b);
    ds = ds - ds_b;
  }
}


void rng_mvhyper2strains(const int n[], int sum, int k, int **x){
  int m = 10;
  int n_otr;
  unsigned int n1;
  unsigned int n2;
  unsigned int t;
  unsigned int res;
  n_otr = sum - n[0];
  n1 = n[0];
  n2 = n_otr;
  t = k;
  res = R::rhyper(n1, n2, t);
  //res = gsl_ran_hypergeometric(r, n1, n2, t);
  (*x)[0] = res;
  for (int i = 1; i < m-1; i++) {
    n_otr = n_otr - n[i];
    n2 = n_otr;
    n1 = n[i];
    k = k - (*x)[i-1];
    t = k;
    res = R::rhyper(n1, n2, t);
    //res = gsl_ran_hypergeometric(r, n1, n2, t);
    (*x)[i] = res;
  }

  (*x)[m-1]=k-(*x)[m-2];
}

AMNLocation2strains::AMNLocation2strains(string name_, int Shome){
  name = name_;
  S = Shome;
  E1 = 0;
  E2 = 0;
  I = 0;
  Ia = 0;
  R = 0;
  N = 0;
  E1_b = 0;
  E2_b = 0;
  I_b = 0;
  Ia_b = 0;
}


void AMNLocation2strains::print(){
  // This is a debug function
  int num;
  Rcout << "AMLocation " << name << " with S= "<< S << ", E1=" << E1 << ", E2=" << E2 << ", I=" << I << ", Ia=" << Ia << ", R= " << R << ", E1_b= " << E1_b << ", E2_b= " << E2_b << ", I_b= " << I_b << ", Ia_b = "  << Ia_b << ". " << endl;
  num = visitorsS.size();
  for (int j = 0; j < num; ++j){
    Rcout << "S " << visitorsS[j] << endl;
    Rcout << "E1 " << visitorsE1[j] << endl;
    Rcout << "E2 " << visitorsE2[j] << endl;
    Rcout << "I " << visitorsI[j] << endl;
    Rcout << "Ia " << visitorsIa[j] << endl;
    Rcout << "R " << visitorsR[j] << endl;
    Rcout << "E1_b " << visitorsE1_b[j] << endl;
    Rcout << "E2_b " << visitorsE2_b[j] << endl;
    Rcout << "I_b " << visitorsI_b[j] << endl;
    Rcout << "Ia_b " << visitorsIa_b[j] << endl;
  }
}

void AMNLocation2strains::seir_step(
    // Function to run a day time step of the model.
    // The commuters are all sent to their work locations.
    float beta,
    float a1,
    float a2,
    float gamma,
    float relativeInfectiousnessB,
    float presymptomaticRelativeInfectiousness,
    float asymptomaticProb,
    float asymptomaticRelativeInfectiousness,
    int &de2,
    int &dea,
    int &de2_b,
    int &dea_b){ // Return symptomatic and asymptomatic incidence
  int S_tmp = S;
  int E1_tmp = E1;
  int E2_tmp = E2;
  int Ia_tmp = Ia;
  int I_tmp = I;
  int R_tmp = R;
  int E1_b_tmp = E1_b;
  int E2_b_tmp = E2_b;
  int I_b_tmp = I_b;
  int Ia_b_tmp = Ia_b;
  double pop_tmp = S + E1 + E2 + Ia + I + R + E2_b + E1_b + I_b + Ia_b;
  int num = visitorsS.size();

  //Vectors with the number of people in the respective compartment,
  //first the the visitors, and then the home population
  //used to distribute the transitions between compartments
  //between the commuters on the different links and the home population
  //the visitors are on the first num elements,
  //the home population last
  int *S_probs = new int[num + 1];
  int *E1_probs = new int[num + 1];
  int *E2_probs = new int[num + 1];
  int *I_probs = new int[num + 1];
  int *Ia_probs = new int[num + 1];
  int *R_probs = new int[num + 1];
  int *E1_b_probs = new int[num + 1];
  int *E2_b_probs = new int[num + 1];
  int *I_b_probs = new int[num + 1];
  int *Ia_b_probs = new int[num + 1];

  for (int i = 0; i < num; ++i){
    S_tmp += visitorsS[i];
    E1_tmp += visitorsE1[i];
    E2_tmp += visitorsE2[i];
    Ia_tmp += visitorsIa[i];
    I_tmp += visitorsI[i];
    R_tmp += visitorsR[i];
    E1_b_tmp += visitorsE1_b[i];
    E2_b_tmp += visitorsE2_b[i];
    Ia_b_tmp += visitorsIa_b[i];
    I_b_tmp += visitorsI_b[i];
    pop_tmp += visitorsS[i] + visitorsE1[i] + visitorsE2[i] + visitorsIa[i] + visitorsI[i] + visitorsR[i] + visitorsE1_b[i] + visitorsE2_b[i] + visitorsIa_b[i] + visitorsI_b[i];

    S_probs[i] = visitorsS[i];
    I_probs[i] = visitorsI[i];
    Ia_probs[i] = visitorsIa[i];
    E1_probs[i] = visitorsE1[i];
    E2_probs[i] = visitorsE2[i];
    I_b_probs[i] = visitorsI_b[i];
    Ia_b_probs[i] = visitorsIa_b[i];
    E1_b_probs[i] = visitorsE1_b[i];
    E2_b_probs[i] = visitorsE2_b[i];
    R_probs[i] = visitorsR[i];
  }

  S_probs[num] = S;
  I_probs[num] = I;
  Ia_probs[num] = Ia;
  E1_probs[num] = E1;
  E2_probs[num] = E2;
  R_probs[num] = R;
  I_b_probs[num] = I_b;
  Ia_b_probs[num] = Ia_b;
  E1_b_probs[num] = E1_b;
  E2_b_probs[num] = E2_b;
  int ds; int de1e2; int de1ia; int de2i; int dia; int di; int ds_b; int de1e2_b; int de1ia_b; int de2i_b; int dia_b; int di_b;

  // Run the SEIR step
  se1e2iiar_sim(ds, de1e2, de1ia, de2i, dia, di, ds_b, de1e2_b, de1ia_b, de2i_b, dia_b, di_b, S_tmp, E1_tmp, E2_tmp, Ia_tmp, I_tmp, E1_b_tmp, E2_b_tmp, Ia_b_tmp, I_b_tmp, beta, a1, a2, gamma, relativeInfectiousnessB, presymptomaticRelativeInfectiousness, asymptomaticProb, asymptomaticRelativeInfectiousness, pop_tmp, 6.0/24.0);
  de2 = de2i;
  dea = de1ia;
  de2_b = de2i_b;
  dea_b = de1ia_b;
  //Distribute the transitions
  double *probs = new double[num + 1];
  double *probs_cum = new double[num + 1];
  double randomnumber;
  int index = -1;
  for(int i = 0; i < dia; ++i){
    for(int k = 0; k < num + 1; ++k) probs[k] = Ia_probs[k] * 1.0 / Ia_tmp; // Oneline for-loop
    randomnumber = rand() * 1.0 / RAND_MAX;
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = rand() * 1.0 / RAND_MAX;
    }
    cumulative_sum2strains(&probs_cum, &probs, num + 1);
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
      // Add to those who are visiting long term
      visitorsIa[index] -= 1;
      visitorsR[index] += 1;
    }

    Ia_probs[index] -= 1;
    Ia_tmp -= 1;
  }

  for(int i = 0; i < dia_b; ++i){
    for(int k = 0; k < num + 1; ++k) probs[k] = Ia_b_probs[k] * 1.0 / Ia_b_tmp; // Oneline for-loop
    randomnumber = rand() * 1.0 / RAND_MAX;
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = rand() * 1.0 / RAND_MAX;
    }
    cumulative_sum2strains(&probs_cum, &probs, num + 1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      Ia_b -= 1;
      R += 1;
    }
    else{
      // Add to those who are visiting long term
      visitorsIa_b[index] -= 1;
      visitorsR[index] += 1;
    }

    Ia_b_probs[index] -= 1;
    Ia_b_tmp -= 1;
  }

  for(int i = 0; i < di; ++i){
    for(int k = 0; k < num + 1; ++k) probs[k] = I_probs[k] * 1.0 / I_tmp; // Oneline for-loop
    randomnumber = rand() * 1.0 / RAND_MAX;
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = rand() * 1.0 / RAND_MAX;
    }
    cumulative_sum2strains(&probs_cum, &probs, num + 1);
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
    else{
      visitorsI[index] -= 1;
      visitorsR[index] += 1;
    }
    I_probs[index] -= 1;
    I_tmp -= 1;
  }

  for(int i = 0; i < di_b; ++i){
    for(int k = 0; k < num + 1; ++k) probs[k] = I_b_probs[k] * 1.0 / I_b_tmp; // Oneline for-loop
    randomnumber = rand() * 1.0 / RAND_MAX;
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = rand() * 1.0 / RAND_MAX;
    }
    cumulative_sum2strains(&probs_cum, &probs, num + 1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      I_b -= 1;
      R += 1;
    }
    else{
      visitorsI_b[index] -= 1;
      visitorsR[index] += 1;
    }
    I_b_probs[index] -= 1;
    I_b_tmp -= 1;
  }

  for(int i = 0; i < de2i; ++i){
    for(int k = 0; k < num + 1; ++k) probs[k] = E2_probs[k] * 1.0 / E2_tmp; // Oneline for-loop
    randomnumber = rand() * 1.0 / RAND_MAX;
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = rand() * 1.0 / RAND_MAX;
    }
    cumulative_sum2strains(&probs_cum, &probs, num + 1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      E2 -= 1;
      I += 1;
    }
    else{
      visitorsE2[index] -= 1;
      visitorsI[index] += 1;
    }

    E2_probs[index] -= 1;
    E2_tmp -= 1;

  }

  for(int i = 0; i < de2i_b; ++i){
    for(int k = 0; k < num + 1; ++k) probs[k] = E2_b_probs[k] * 1.0 / E2_b_tmp; // Oneline for-loop
    randomnumber = rand() * 1.0 / RAND_MAX;
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = rand() * 1.0 / RAND_MAX;
    }
    cumulative_sum2strains(&probs_cum, &probs, num + 1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      E2_b -= 1;
      I_b += 1;
    }
    else{
      visitorsE2_b[index] -= 1;
      visitorsI_b[index] += 1;
    }

    E2_b_probs[index] -= 1;
    E2_b_tmp -= 1;

  }

  for(int i = 0; i < de1e2; ++i){
    for(int k = 0; k < num + 1; ++k) probs[k] = E1_probs[k] * 1.0 / E1_tmp; // Oneline for-loop
    randomnumber = rand() * 1.0 / RAND_MAX;
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = rand() * 1.0 / RAND_MAX;
    }
    cumulative_sum2strains(&probs_cum, &probs, num+1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      E1 -= 1;
      E2 += 1;
    }

    else{
      visitorsE1[index] -= 1;
      visitorsE2[index] += 1;
    }
    E1_probs[index] -= 1;
    E1_tmp -= 1;
  }


  for(int i = 0; i < de1e2_b; ++i){
    for(int k = 0; k < num + 1; ++k) probs[k] = E1_b_probs[k] * 1.0 / E1_b_tmp; // Oneline for-loop
    randomnumber = rand() * 1.0 / RAND_MAX;
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = rand() * 1.0 / RAND_MAX;
    }
    cumulative_sum2strains(&probs_cum, &probs, num+1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      E1_b -= 1;
      E2_b += 1;
    }

    else{
      visitorsE1_b[index] -= 1;
      visitorsE2_b[index] += 1;
    }
    E1_b_probs[index] -= 1;
    E1_b_tmp -= 1;
  }

  for(int i = 0; i < de1ia; ++i){
    for(int k = 0; k < num + 1; ++k) probs[k] = E1_probs[k] * 1.0 / E1_tmp; // Oneline for-loop
    randomnumber = rand() * 1.0 / RAND_MAX;
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = rand() * 1.0 / RAND_MAX;
    }
    cumulative_sum2strains(&probs_cum, &probs, num+1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      E1 -= 1;
      Ia += 1;
    }

    else{
      visitorsE1[index] -= 1;
      visitorsIa[index] += 1;
    }
    E1_probs[index] -= 1;
    E1_tmp -= 1;
  }

  for(int i = 0; i < de1ia_b; ++i){
    for(int k = 0; k < num + 1; ++k) probs[k] = E1_b_probs[k] * 1.0 / E1_b_tmp; // Oneline for-loop
    randomnumber = rand() * 1.0 / RAND_MAX;
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = rand() * 1.0 / RAND_MAX;
    }
    cumulative_sum2strains(&probs_cum, &probs, num+1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      E1_b -= 1;
      Ia_b += 1;
    }

    else{
      visitorsE1_b[index] -= 1;
      visitorsIa_b[index] += 1;
    }
    E1_b_probs[index] -= 1;
    E1_b_tmp -= 1;
  }

  for(int i = 0; i < ds; ++i){
    for(int k = 0; k < num + 1; ++k) probs[k] = S_probs[k] * 1.0 / S_tmp; // Oneline for-loop
    randomnumber = rand() * 1.0 / RAND_MAX;
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = rand() * 1.0/RAND_MAX;
    }
    cumulative_sum2strains(&probs_cum, &probs, num + 1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      S -= 1;
      E1 += 1;
    }

    else{
      visitorsS[index] -= 1;
      visitorsE1[index] += 1;
    }
    S_probs[index] -= 1;
    S_tmp -= 1;
  }

  for(int i = 0; i < ds_b; ++i){
    for(int k = 0; k < num + 1; ++k) probs[k] = S_probs[k] * 1.0 / S_tmp; // Oneline for-loop
    randomnumber = rand() * 1.0 / RAND_MAX;
    while(randomnumber == 0 || randomnumber == 1){
      randomnumber = rand() * 1.0/RAND_MAX;
    }
    cumulative_sum2strains(&probs_cum, &probs, num + 1);
    for(int h = 0; h < num + 1; ++h){
      if(randomnumber < probs_cum[h]){
        index = h;
        break;
      }
    }
    if(index == num){
      // Add to home pop
      S -= 1;
      E1_b += 1;
    }

    else{
      visitorsS[index] -= 1;
      visitorsE1_b[index] += 1;
    }
    S_probs[index] -= 1;
    S_tmp -= 1;
  }


  delete[] probs;
  delete[] probs_cum;
  delete[] S_probs;
  delete[] E1_probs;
  delete[] E2_probs;
  delete[] I_probs;
  delete[] Ia_probs;
  delete[] R_probs;
  delete[] E1_b_probs;
  delete[] E2_b_probs;
  delete[] I_b_probs;
  delete[] Ia_b_probs;
}


AMNGraph2strains::AMNGraph2strains(){
}

void AMNGraph2strains::add_node(string name, int Shome){
  AMNLocation2strains newlocation(name, Shome);
  locations.push_back(newlocation);
}




void AMNGraph2strains::print(){
  Rcout << endl << "Printing graph: " << endl;
  Rcout << "AMLocations: " << endl;
  for(int i = 0; i < locations.size(); ++i){
    locations[i].print();
    Rcout.flush();
  }

  Rcout << endl << endl;;
}

void AMNGraph2strains::copy_graph(AMNGraph2strains G){
  int counter = 0;
  for(vector<AMNLocation2strains>::iterator it = G.locations.begin() ; it != G.locations.end(); ++it){
    add_node(it->name, it->S);
    locations[counter].visitorsS = it->visitorsS;
    locations[counter].visitorsE1 = it->visitorsE1;
    locations[counter].visitorsE2 = it->visitorsE2;
    locations[counter].visitorsI = it->visitorsI;
    locations[counter].visitorsIa = it->visitorsIa;
    locations[counter].visitorsR = it->visitorsR;
    locations[counter].visitorsE1_b = it->visitorsE1_b;
    locations[counter].visitorsE2_b = it->visitorsE2_b;
    locations[counter].visitorsI_b = it->visitorsI_b;
    locations[counter].visitorsIa_b = it->visitorsIa_b;

    locations[counter].S = it->S;
    locations[counter].E1 = it->E1;
    locations[counter].E2 = it->E2;
    locations[counter].I = it->I;
    locations[counter].Ia = it->Ia;
    locations[counter].R = it->R;
    locations[counter].E1_b = it->E1_b;
    locations[counter].E2_b = it->E2_b;
    locations[counter].I_b = it->I_b;
    locations[counter].Ia_b = it->Ia_b;

    counter += 1;
  }
}

void AMNGraph2strains::count_everyone(string msg){
  int sum = 0;
  for(int i = 0; i < locations.size(); ++i){
    sum += locations[i].S + locations[i].E1 + locations[i].E2 + locations[i].I + locations[i].Ia + locations[i].R + locations[i].E1_b + locations[i].E2_b + locations[i].I_b + locations[i].Ia_b;
    for (int j = 0; j < locations.size(); ++j){
      sum+= locations[i].visitorsS[j] + locations[i].visitorsE1[j] + locations[i].visitorsE2[j] + locations[i].visitorsI[j] + locations[i].visitorsIa[j] + locations[i].visitorsR[j] + locations[i].visitorsE1_b[j] + locations[i].visitorsE2_b[j] + locations[i].visitorsI_b[j] + locations[i].visitorsIa_b[j];
    }
  }

  Rcout << msg << ", SUM IS " << sum << endl;
}

//' asymmetric_mobility_se1e2iiar_2strains_cpp
//'
//' Raw CPP function. Should not be called directly.
//'
//' @param se1e2iiar_2strains_pop Data.frame
//' @param mobility_matrix List of data.frames
//' @param seed_matrix matrix of seeding cases of strain a per date per geographical location
//' @param seed_matrix_b matrix of seeding cases of strain b per date per geographical location
//' @param betas matrix of floats, number of time intervals times number of locations, infection parameter, 0.6
//' @param inputSeed Integer, input seed
//' @param a1 Float, 1/latent period, 1/2.0
//' @param a2 Float, 1/presymptomatic period, 1/3.0
//' @param gamma Float, 1/infectious period, 1/5.0
//' @param relativeInfectiousnessB Float, Relative infectiousness of strain B
//' @param presymptomaticRelativeInfectiousness Float, Relative infectiousness of presymptomatic infectious
//' @param asymptomaticProb Float, Proportion/probability of asymptomatic given infectious
//' @param asymptomaticRelativeInfectiousness Float, Relative infectiousness of asymptomatic infectious
//' @param N Int = 1 int, Number of repetitions
//' @param M Int, Number of days
//' @export
// [[Rcpp::export]]
List asymmetric_mobility_se1e2iiar_2strains_cpp(
    DataFrame se1e2iiar_2strains_pop,
    List mobility_matrix,
    NumericMatrix seed_matrix,
    NumericMatrix seed_matrix_b,
    NumericMatrix betas,
    int inputSeed,
    float a1,
    float a2,
    float gamma,
    float relativeInfectiousnessB,
    float presymptomaticRelativeInfectiousness,
    float asymptomaticProb,
    float asymptomaticRelativeInfectiousness,
    int N=1,
    int M=56
){
  int n=0; //Number of locations
  // return a new data frame
  DataFrame empty = DataFrame::create(
    _["empty"]= 1
  );

  StringVector names = se1e2iiar_2strains_pop[0] ;
  IntegerVector pop_S = se1e2iiar_2strains_pop[1] ;
  IntegerVector pop_E1 = se1e2iiar_2strains_pop[2] ;
  IntegerVector pop_E2 = se1e2iiar_2strains_pop[3] ;
  IntegerVector pop_I = se1e2iiar_2strains_pop[4] ;
  IntegerVector pop_Ia = se1e2iiar_2strains_pop[5] ;
  IntegerVector pop_R = se1e2iiar_2strains_pop[6] ;
  IntegerVector pop_E1_b = se1e2iiar_2strains_pop[7] ;
  IntegerVector pop_E2_b = se1e2iiar_2strains_pop[8] ;
  IntegerVector pop_I_b = se1e2iiar_2strains_pop[9] ;
  IntegerVector pop_Ia_b = se1e2iiar_2strains_pop[10] ;

  AMNGraph2strains G;
  string tmpstr;
  string tmpstr2;
  int pop;

  for (int i = 0; i < se1e2iiar_2strains_pop.rows(); i++) {
    string name = std::string(names[i]);
    G.add_node(name, pop_S[i]+pop_E1[i] + pop_E2[i] +pop_I[i]+pop_Ia[i]+pop_R[i]+pop_E1_b[i] + pop_E2_b[i] +pop_I_b[i]+pop_Ia_b[i]);
    n+= 1;
  }

  // Add the visitors
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j){
      G.locations[i].visitorsS.push_back(0);
      G.locations[i].visitorsE1.push_back(0);
      G.locations[i].visitorsE2.push_back(0);
      G.locations[i].visitorsI.push_back(0);
      G.locations[i].visitorsIa.push_back(0);
      G.locations[i].visitorsR.push_back(0);
      G.locations[i].visitorsE1_b.push_back(0);
      G.locations[i].visitorsE2_b.push_back(0);
      G.locations[i].visitorsI_b.push_back(0);
      G.locations[i].visitorsIa_b.push_back(0);
    }
  }

  int S; int E1 = 0; int E2 = 0; int I = 0; int Ia = 0; int R = 0; int E1_b = 0; int E2_b = 0; int I_b = 0; int Ia_b = 0;

  // Add the travellers for the first time step,
  // they are assumed to be susceptible
  int name1_index;
  int name2_index;
  int safecount = 0;

  DataFrame mtrx = DataFrame(mobility_matrix[0]);
  StringVector mtrx_from = mtrx["from"];
  StringVector mtrx_to = mtrx["to"];
  NumericVector mtrx_n = mtrx["n"];

  for (int i_m = 0; i_m < mtrx.rows(); i_m++) {
    tmpstr = mtrx_from[i_m];
    tmpstr2 = mtrx_to[i_m];
    S = mtrx_n[i_m];

    safecount += 1;
    name1_index = -1;
    name2_index = -1;
    for (int i = 0; i < G.locations.size(); ++i){
      if(G.locations[i].name.compare(tmpstr) == 0){
        name1_index = i;
      }
      if(G.locations[i].name.compare(tmpstr2) == 0){
        name2_index = i;
      }
    }
    if (name1_index == -1 || name2_index == -1){
      Rcout << "Error in add edge with input: " << tmpstr << ", " << tmpstr2 << ", " << S << ", " << E1 << ", " << E2 <<  ", " << I << ", " << Ia << ", " << R << E1_b << E2_b << I_b << Ia_b << endl;
      Rcout << "name1_index = " << name1_index << ", name_index2 = " << name2_index << endl;
      return(empty);
    }
    else{
      G.locations[name2_index].visitorsS[name1_index] = S;
      G.locations[name1_index].S -= S;
    }
  }

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
  /// Third index 0=S, 1=E1, 2 = E2, 3=I, 4=Ia, 5=R; 6=E1_b; 7 = E2_b; 8 = I_b; 9 = Ia_b for belonging to kommune
  /// 10 = S, 11 = E1, 12 = E2,  13 = I, 14 = Ia, 15 = R, 16 = E1_b, 17 = E2_b, 18 = I_b, 19 = Ia_b, for currently in kommune.
  /// 20 = symptomatic incidence occurring in a kommune of strain A;
  /// 21 = asymptomatic incidence occuring in a kommune of strain A;
  /// 22 = symptomatic incidence occurring in a kommune of strain B;
  /// 23 = asymptomatic incidence occuring in a kommune of strain B;
  /// 24 = symptomatic incidence occurring in a kommune total;
  /// 25 = asymptomatic incidence occuring in a kommune total;
  for(int i = 0; i < n; ++i){
    values[i] = new int*[M * 4]; // 4 * M, stored for all 6-hour intervals.
    peak_date[i] = new int[N];
    peak_val[i] = new int[N];
    start_date[i] = new int[N];
    I_this[i] = new int[4 * M];
    final_size[i] = new int[N];
    for(int k = 0; k < 4 * M; ++k){
      values[i][k] = new int[26];
      values[i][k][0] = 0;
      values[i][k][1] = 0;
      values[i][k][2] = 0;
      values[i][k][3] = 0;
      values[i][k][4] = 0;
      values[i][k][5] = 0;
      values[i][k][6] = 0;
      values[i][k][7] = 0;
      values[i][k][8] = 0;
      values[i][k][9] = 0;
      values[i][k][10] = 0;
      values[i][k][11] = 0;
      values[i][k][12] = 0;
      values[i][k][13] = 0;
      values[i][k][14] = 0;
      values[i][k][15] = 0;
      values[i][k][16] = 0;
      values[i][k][17] = 0;
      values[i][k][18] = 0;
      values[i][k][19] = 0;
      values[i][k][20] = 0;
      values[i][k][21] = 0;
      values[i][k][22] = 0;
      values[i][k][23] = 0;
      values[i][k][24] = 0;
      values[i][k][25] = 0;
    }
  }

  for (int i = 0; i < N; ++i){
    bonds[i] = new int[4 * M]; //4 * M, stored for all 6-hour intervals.
    for (int j = 0; j < 4 * M; ++j){
      bonds[i][j] = 0;
    }
  }

  unsigned int Seed2 = inputSeed; // for random seed.
  srand(Seed2);
  /*
   //extern const gsl_rng_type *gsl_rng_default;
   default_random_engine generator(time(0));
   gsl_rng_env_setup();
   const gsl_rng_type * T2;
   gsl_rng * r2;
   //T2 = gsl_rng_default;
   r2 = gsl_rng_alloc( gsl_rng_borosh13);
   gsl_rng_set (r2, Seed2);
   */
  int Nk; // Number of travellers on link on current day.
  float sum;
  int Nk_prev; // Number of travellers and long term stayers on link on previous day.
  int leftover; // Number of people we have to take from the other kommuner

  int p[10]; //To distribute travellers according to proportion in the six compartments
  int *x = new int [10];

  for(int i_sim = 0; i_sim < N; ++i_sim){
    for (int i = 0; i < n; ++i){
      final_size[i][i_sim] = 0;
      for (int k = 0; k < M * 4; ++k){
        I_this[i][k] = 0;
      }
    }
    AMNGraph2strains G_current;

    G_current.copy_graph(G);

    // Seed the epidemic
    for(int i = 0; i < n; ++i){
      G_current.locations[i].S -= pop_E1[i] + pop_E2[i] + pop_I[i] + pop_Ia[i] + pop_R[i] + pop_E1_b[i] + pop_E2_b[i] + pop_I_b[i] + pop_Ia_b[i];
      G_current.locations[i].E1  = pop_E1[i];
      G_current.locations[i].E2  = pop_E2[i];
      G_current.locations[i].I  = pop_I[i];
      G_current.locations[i].Ia = pop_Ia[i];
      G_current.locations[i].R  = pop_R[i];
      G_current.locations[i].E1_b  = pop_E1_b[i];
      G_current.locations[i].E2_b  = pop_E2_b[i];
      G_current.locations[i].I_b  = pop_I_b[i];
      G_current.locations[i].Ia_b = pop_Ia_b[i];
    }

    for(int i_t = 0; i_t < 4 * M; ++i_t){
      if (i_t%4 == 0){
        for (int i = 0; i < n; ++i){
          if (G_current.locations[i].S != 0){
            G_current.locations[i].S -= seed_matrix(i_t/4, i);
            G_current.locations[i].I += seed_matrix(i_t/4, i);
            G_current.locations[i].I_b += seed_matrix_b(i_t/4, i);
          }
        }
      }
      for (int i = 0; i < n; ++i){
        int de2 = 0;
        int dea = 0;
        int de2_b = 0;
        int dea_b = 0;
        G_current.locations[i].seir_step(betas(i_t, i), a1, a2, gamma, relativeInfectiousnessB, presymptomaticRelativeInfectiousness, asymptomaticProb, asymptomaticRelativeInfectiousness, de2, dea, de2_b, dea_b);
        values[i][i_t][0] += G_current.locations[i].S;
        values[i][i_t][1] += G_current.locations[i].E1;
        values[i][i_t][2] += G_current.locations[i].E2;
        values[i][i_t][3] += G_current.locations[i].I;
        values[i][i_t][4] += G_current.locations[i].Ia;
        values[i][i_t][5] += G_current.locations[i].R;
        values[i][i_t][6] += G_current.locations[i].E1_b;
        values[i][i_t][7] += G_current.locations[i].E2_b;
        values[i][i_t][8] += G_current.locations[i].I_b;
        values[i][i_t][9] += G_current.locations[i].Ia_b;
        values[i][i_t][10] += G_current.locations[i].S;
        values[i][i_t][11] += G_current.locations[i].E1;
        values[i][i_t][12] += G_current.locations[i].E2;
        values[i][i_t][13] += G_current.locations[i].I;
        values[i][i_t][14] += G_current.locations[i].Ia;
        values[i][i_t][15] += G_current.locations[i].R;
        values[i][i_t][16] += G_current.locations[i].E1_b;
        values[i][i_t][17] += G_current.locations[i].E2_b;
        values[i][i_t][18] += G_current.locations[i].I_b;
        values[i][i_t][19] += G_current.locations[i].Ia_b;
        I_this[i][i_t] += G_current.locations[i].I + G_current.locations[i].I_b;

        if(i_t == (4 * M-1)){
          final_size[i][i_sim] += G_current.locations[i].R;
        }

        values[i][i_t][20] += de2;
        values[i][i_t][21] += dea;
        values[i][i_t][22] += de2_b;
        values[i][i_t][23] += dea_b;
        values[i][i_t][24] += de2 + de2_b;
        values[i][i_t][25] += dea + dea_b;
        bonds[i_sim][i_t] += G_current.locations[i].I + G_current.locations[i].Ia + G_current.locations[i].I_b + G_current.locations[i].Ia_b;

        int num = G_current.locations[i].visitorsS.size();
        for (int j = 0; j < num; ++j){
          values[j][i_t][0] += G_current.locations[i].visitorsS[j];
          values[j][i_t][1] += G_current.locations[i].visitorsE1[j];
          values[j][i_t][2] += G_current.locations[i].visitorsE2[j];
          values[j][i_t][3] += G_current.locations[i].visitorsI[j];
          values[j][i_t][4] += G_current.locations[i].visitorsIa[j];
          values[j][i_t][5] += G_current.locations[i].visitorsR[j];
          values[j][i_t][6] += G_current.locations[i].visitorsE1_b[j];
          values[j][i_t][7] += G_current.locations[i].visitorsE2_b[j];
          values[j][i_t][8] += G_current.locations[i].visitorsI_b[j];
          values[j][i_t][9] += G_current.locations[i].visitorsIa_b[j];
          I_this[j][i_t] += G_current.locations[i].visitorsI[j] + G_current.locations[i].visitorsI_b[j];
          bonds[i_sim][i_t] += G_current.locations[i].visitorsI[j] + G_current.locations[i].visitorsIa[j] + G_current.locations[i].visitorsI_b[j] + G_current.locations[i].visitorsIa_b[j];
          values[i][i_t][10] += G_current.locations[i].visitorsS[j];
          values[i][i_t][11] += G_current.locations[i].visitorsE1[j];
          values[i][i_t][12] += G_current.locations[i].visitorsE2[j];
          values[i][i_t][13] += G_current.locations[i].visitorsI[j];
          values[i][i_t][14] += G_current.locations[i].visitorsIa[j];
          values[i][i_t][15] += G_current.locations[i].visitorsR[j];
          values[i][i_t][16] += G_current.locations[i].visitorsE1_b[j];
          values[i][i_t][17] += G_current.locations[i].visitorsE2_b[j];
          values[i][i_t][18] += G_current.locations[i].visitorsI_b[j];
          values[i][i_t][19] += G_current.locations[i].visitorsIa_b[j];
        }
      }

      AMNGraph2strains G_prev; // AMGraph at previous time point
      G_prev.copy_graph(G_current);

      int safecount = 0;

      if(i_t < (4 * M - 1)){
        mtrx = mobility_matrix[i_t + 1];
        mtrx_from = mtrx["from"];
        mtrx_to = mtrx["to"];
        mtrx_n = mtrx["n"];

        for (int i_m = 0; i_m < mtrx.rows(); i_m++) {

          tmpstr = mtrx_from[i_m];
          tmpstr2 = mtrx_to[i_m];
          S = mtrx_n[i_m];

          int name1_index = -1;
          int name2_index = -1;
          safecount += 1;
          for (int i = 0; i < G_current.locations.size(); ++i){
            if(G_current.locations[i].name.compare(tmpstr) == 0){
              name1_index = i;
            }
            if(G_current.locations[i].name.compare(tmpstr2) == 0){
              name2_index = i;
            }
          }
          if (name1_index == -1 || name2_index == -1){
            Rcout << "Error in add edge with input: " << tmpstr << ", " << tmpstr2 << ", " << S << ", " << E1 << ", " << E2 << ", " << I << ", " << Ia << ", " << R << ", " << E1_b << ", " << E2_b << ", " << I_b << ", " << Ia_b << endl;
            Rcout << "name1_index = " << name1_index << ", name2_index = " << name2_index << endl;
            Rcout << "Problem with filecounter " << i_t << endl;
            return(empty);
          }
          else{
            // Fill commuter links
            //Make new outlinks
            // Count number of people from name2 that are currently in name 1, those we wish to send home
            Nk_prev = G_prev.locations[name1_index].visitorsS[name2_index] + G_prev.locations[name1_index].visitorsE1[name2_index] + G_prev.locations[name1_index].visitorsE2[name2_index] + G_prev.locations[name1_index].visitorsI[name2_index] + G_prev.locations[name1_index].visitorsIa[name2_index] + G_prev.locations[name1_index].visitorsR[name2_index]  + G_prev.locations[name1_index].visitorsE1_b[name2_index] + G_prev.locations[name1_index].visitorsE2_b[name2_index] + G_prev.locations[name1_index].visitorsI_b[name2_index] + G_prev.locations[name1_index].visitorsIa_b[name2_index];
            Nk = S;
            if(Nk_prev == Nk){
              G_current.locations[name2_index].S +=  G_prev.locations[name1_index].visitorsS[name2_index];
              G_current.locations[name2_index].E1 +=  G_prev.locations[name1_index].visitorsE1[name2_index];
              G_current.locations[name2_index].E2 +=  G_prev.locations[name1_index].visitorsE2[name2_index];
              G_current.locations[name2_index].I +=  G_prev.locations[name1_index].visitorsI[name2_index];
              G_current.locations[name2_index].Ia +=  G_prev.locations[name1_index].visitorsIa[name2_index];
              G_current.locations[name2_index].R +=  G_prev.locations[name1_index].visitorsR[name2_index];
              G_current.locations[name2_index].E1_b +=  G_prev.locations[name1_index].visitorsE1_b[name2_index];
              G_current.locations[name2_index].E2_b +=  G_prev.locations[name1_index].visitorsE2_b[name2_index];
              G_current.locations[name2_index].I_b +=  G_prev.locations[name1_index].visitorsI_b[name2_index];
              G_current.locations[name2_index].Ia_b +=  G_prev.locations[name1_index].visitorsIa_b[name2_index];

              G_current.locations[name1_index].visitorsS[name2_index] -= G_prev.locations[name1_index].visitorsS[name2_index];
              G_current.locations[name1_index].visitorsE1[name2_index] -= G_prev.locations[name1_index].visitorsE1[name2_index];
              G_current.locations[name1_index].visitorsE2[name2_index] -= G_prev.locations[name1_index].visitorsE2[name2_index];
              G_current.locations[name1_index].visitorsI[name2_index] -= G_prev.locations[name1_index].visitorsI[name2_index];
              G_current.locations[name1_index].visitorsIa[name2_index] -= G_prev.locations[name1_index].visitorsIa[name2_index];
              G_current.locations[name1_index].visitorsR[name2_index] -= G_prev.locations[name1_index].visitorsR[name2_index];
              G_current.locations[name1_index].visitorsE1_b[name2_index] -= G_prev.locations[name1_index].visitorsE1_b[name2_index];
              G_current.locations[name1_index].visitorsE2_b[name2_index] -= G_prev.locations[name1_index].visitorsE2_b[name2_index];
              G_current.locations[name1_index].visitorsI_b[name2_index] -= G_prev.locations[name1_index].visitorsI_b[name2_index];
              G_current.locations[name1_index].visitorsIa_b[name2_index] -= G_prev.locations[name1_index].visitorsIa_b[name2_index];

              G_prev.locations[name1_index].visitorsS[name2_index] = 0;
              G_prev.locations[name1_index].visitorsE1[name2_index] = 0;
              G_prev.locations[name1_index].visitorsE2[name2_index] = 0;
              G_prev.locations[name1_index].visitorsI[name2_index] = 0;
              G_prev.locations[name1_index].visitorsIa[name2_index] = 0;
              G_prev.locations[name1_index].visitorsR[name2_index] = 0;
              G_prev.locations[name1_index].visitorsE1_b[name2_index] = 0;
              G_prev.locations[name1_index].visitorsE2_b[name2_index] = 0;
              G_prev.locations[name1_index].visitorsI_b[name2_index] = 0;
              G_prev.locations[name1_index].visitorsIa_b[name2_index] = 0;
            }

            else if(Nk_prev < Nk){
              // We send everyone home, plus the additional people, selected at random from host population
              G_current.locations[name2_index].S +=  G_prev.locations[name1_index].visitorsS[name2_index];
              G_current.locations[name2_index].E1 +=  G_prev.locations[name1_index].visitorsE1[name2_index];
              G_current.locations[name2_index].E2 +=  G_prev.locations[name1_index].visitorsE2[name2_index];
              G_current.locations[name2_index].I +=  G_prev.locations[name1_index].visitorsI[name2_index];
              G_current.locations[name2_index].Ia +=  G_prev.locations[name1_index].visitorsIa[name2_index];
              G_current.locations[name2_index].R +=  G_prev.locations[name1_index].visitorsR[name2_index];
              G_current.locations[name2_index].E1_b +=  G_prev.locations[name1_index].visitorsE1_b[name2_index];
              G_current.locations[name2_index].E2_b +=  G_prev.locations[name1_index].visitorsE2_b[name2_index];
              G_current.locations[name2_index].I_b +=  G_prev.locations[name1_index].visitorsI_b[name2_index];
              G_current.locations[name2_index].Ia_b +=  G_prev.locations[name1_index].visitorsIa_b[name2_index];

              G_current.locations[name1_index].visitorsS[name2_index] -= G_prev.locations[name1_index].visitorsS[name2_index];
              G_current.locations[name1_index].visitorsE1[name2_index] -= G_prev.locations[name1_index].visitorsE1[name2_index];
              G_current.locations[name1_index].visitorsE2[name2_index] -= G_prev.locations[name1_index].visitorsE2[name2_index];
              G_current.locations[name1_index].visitorsI[name2_index] -= G_prev.locations[name1_index].visitorsI[name2_index];
              G_current.locations[name1_index].visitorsIa[name2_index] -= G_prev.locations[name1_index].visitorsIa[name2_index];
              G_current.locations[name1_index].visitorsR[name2_index] -= G_prev.locations[name1_index].visitorsR[name2_index];
              G_current.locations[name1_index].visitorsE1_b[name2_index] -= G_prev.locations[name1_index].visitorsE1_b[name2_index];
              G_current.locations[name1_index].visitorsE2_b[name2_index] -= G_prev.locations[name1_index].visitorsE2_b[name2_index];
              G_current.locations[name1_index].visitorsI_b[name2_index] -= G_prev.locations[name1_index].visitorsI_b[name2_index];
              G_current.locations[name1_index].visitorsIa_b[name2_index] -= G_prev.locations[name1_index].visitorsIa_b[name2_index];

              G_prev.locations[name1_index].visitorsS[name2_index] = 0;
              G_prev.locations[name1_index].visitorsE1[name2_index] = 0;
              G_prev.locations[name1_index].visitorsE2[name2_index] = 0;
              G_prev.locations[name1_index].visitorsI[name2_index] = 0;
              G_prev.locations[name1_index].visitorsIa[name2_index] = 0;
              G_prev.locations[name1_index].visitorsR[name2_index] = 0;
              G_prev.locations[name1_index].visitorsE1_b[name2_index] = 0;
              G_prev.locations[name1_index].visitorsE2_b[name2_index] = 0;
              G_prev.locations[name1_index].visitorsI_b[name2_index] = 0;
              G_prev.locations[name1_index].visitorsIa_b[name2_index] = 0;

              sum = G_prev.locations[name1_index].S + G_prev.locations[name1_index].E1 + G_prev.locations[name1_index].E2 + G_prev.locations[name1_index].I + G_prev.locations[name1_index].Ia + G_prev.locations[name1_index].R + G_prev.locations[name1_index].E1_b + G_prev.locations[name1_index].E2_b + G_prev.locations[name1_index].I_b + G_prev.locations[name1_index].Ia_b;

              if(sum >= Nk - Nk_prev){
                // Enough people present from home location in kommune
                if (sum == G_prev.locations[name1_index].S){
                  G_current.locations[name2_index].visitorsS[name1_index] += Nk - Nk_prev;
                  G_current.locations[name1_index].S -= Nk - Nk_prev;
                  G_prev.locations[name1_index].S -= Nk - Nk_prev;
                }
                else{
                  p[0] = G_prev.locations[name1_index].S;
                  p[1] = G_prev.locations[name1_index].E1;
                  p[2] = G_prev.locations[name1_index].E2;
                  p[3] = G_prev.locations[name1_index].I;
                  p[4] = G_prev.locations[name1_index].Ia;
                  p[5] = G_prev.locations[name1_index].R;
                  p[6] = G_prev.locations[name1_index].E1_b;
                  p[7] = G_prev.locations[name1_index].E2_b;
                  p[8] = G_prev.locations[name1_index].I_b;
                  p[9] = G_prev.locations[name1_index].Ia_b;
                  // Draw the number of travellers from each compartment
                  rng_mvhyper2strains(p, sum, Nk - Nk_prev, &x);
                  G_current.locations[name2_index].visitorsS[name1_index] += x[0];
                  G_current.locations[name1_index].S -= x[0];
                  G_prev.locations[name1_index].S -= x[0];
                  G_current.locations[name2_index].visitorsE1[name1_index] += x[1];
                  G_current.locations[name1_index].E1 -= x[1];
                  G_prev.locations[name1_index].E1 -= x[1];
                  G_current.locations[name2_index].visitorsE2[name1_index] += x[2];
                  G_current.locations[name1_index].E2 -= x[2];
                  G_prev.locations[name1_index].E2 -= x[2];
                  G_current.locations[name2_index].visitorsI[name1_index] += x[3];
                  G_current.locations[name1_index].I -= x[3];
                  G_prev.locations[name1_index].I -= x[3];
                  G_current.locations[name2_index].visitorsIa[name1_index] += x[4];
                  G_current.locations[name1_index].Ia -= x[4];
                  G_prev.locations[name1_index].Ia -= x[4];
                  G_current.locations[name2_index].visitorsR[name1_index] += x[5];
                  G_current.locations[name1_index].R -= x[5];
                  G_prev.locations[name1_index].R -= x[5];
                  G_current.locations[name2_index].visitorsE1_b[name1_index] += x[6];
                  G_current.locations[name1_index].E1_b -= x[6];
                  G_prev.locations[name1_index].E1_b -= x[6];
                  G_current.locations[name2_index].visitorsE2_b[name1_index] += x[7];
                  G_current.locations[name1_index].E2_b -= x[7];
                  G_prev.locations[name1_index].E2_b -= x[7];
                  G_current.locations[name2_index].visitorsI_b[name1_index] += x[8];
                  G_current.locations[name1_index].I_b -= x[8];
                  G_prev.locations[name1_index].I_b -= x[8];
                  G_current.locations[name2_index].visitorsIa_b[name1_index] += x[9];
                  G_current.locations[name1_index].Ia_b -= x[9];
                  G_prev.locations[name1_index].Ia_b -= x[9];
                }
              }
              else{
                // We have to take from the rest of the people present
                leftover = Nk - Nk_prev - sum;
                G_current.locations[name2_index].visitorsS[name1_index] += G_prev.locations[name1_index].S;
                G_current.locations[name1_index].S -= G_prev.locations[name1_index].S;
                G_prev.locations[name1_index].S = 0;
                G_current.locations[name2_index].visitorsE1[name1_index] += G_prev.locations[name1_index].E1;
                G_current.locations[name1_index].E1 -= G_prev.locations[name1_index].E1;
                G_prev.locations[name1_index].E1 = 0;
                G_current.locations[name2_index].visitorsE2[name1_index] += G_prev.locations[name1_index].E2;
                G_current.locations[name1_index].E2 -= G_prev.locations[name1_index].E2;
                G_prev.locations[name1_index].E2 = 0;
                G_current.locations[name2_index].visitorsI[name1_index] += G_prev.locations[name1_index].I;
                G_current.locations[name1_index].I -= G_prev.locations[name1_index].I;
                G_prev.locations[name1_index].I = 0;
                G_current.locations[name2_index].visitorsIa[name1_index] += G_prev.locations[name1_index].Ia;
                G_current.locations[name1_index].Ia -= G_prev.locations[name1_index].Ia;
                G_prev.locations[name1_index].Ia = 0;
                G_current.locations[name2_index].visitorsR[name1_index] += G_prev.locations[name1_index].R;
                G_current.locations[name1_index].R -= G_prev.locations[name1_index].R;
                G_prev.locations[name1_index].R = 0;
                G_current.locations[name2_index].visitorsE1_b[name1_index] += G_prev.locations[name1_index].E1_b;
                G_current.locations[name1_index].E1_b -= G_prev.locations[name1_index].E1_b;
                G_prev.locations[name1_index].E1_b = 0;
                G_current.locations[name2_index].visitorsE2_b[name1_index] += G_prev.locations[name1_index].E2_b;
                G_current.locations[name1_index].E2_b -= G_prev.locations[name1_index].E2_b;
                G_prev.locations[name1_index].E2_b = 0;
                G_current.locations[name2_index].visitorsI_b[name1_index] += G_prev.locations[name1_index].I_b;
                G_current.locations[name1_index].I_b -= G_prev.locations[name1_index].I_b;
                G_prev.locations[name1_index].I_b = 0;
                G_current.locations[name2_index].visitorsIa_b[name1_index] += G_prev.locations[name1_index].Ia_b;
                G_current.locations[name1_index].Ia_b -= G_prev.locations[name1_index].Ia_b;
                G_prev.locations[name1_index].Ia_b = 0;

                int num = G_prev.locations[name1_index].visitorsS.size();

                //Vectors with the number of people in the respective compartment,
                //as visitors from each location
                //used to distribute the extra travellers
                // As there are not enough people from the home locations currently present in name1_index
                int S_tmp = 0; int E1_tmp = 0; int E2_tmp = 0; int Ia_tmp = 0; int I_tmp = 0; int R_tmp = 0; int E1_b_tmp = 0; int E2_b_tmp = 0; int Ia_b_tmp = 0; int I_b_tmp = 0;

                int *S_probs = new int[num];
                int *E1_probs = new int[num];
                int *E2_probs = new int[num];
                int *I_probs = new int[num];
                int *Ia_probs = new int[num];
                int *R_probs = new int[num];
                int *E1_b_probs = new int[num];
                int *E2_b_probs = new int[num];
                int *I_b_probs = new int[num];
                int *Ia_b_probs = new int[num];

                for(int j = 0; j < num; ++j){
                  S_tmp += G_prev.locations[name1_index].visitorsS[j];
                  E1_tmp += G_prev.locations[name1_index].visitorsE1[j];
                  E2_tmp += G_prev.locations[name1_index].visitorsE2[j];
                  Ia_tmp += G_prev.locations[name1_index].visitorsIa[j];
                  I_tmp += G_prev.locations[name1_index].visitorsI[j];
                  R_tmp += G_prev.locations[name1_index].visitorsR[j];
                  E1_b_tmp += G_prev.locations[name1_index].visitorsE1_b[j];
                  E2_b_tmp += G_prev.locations[name1_index].visitorsE2_b[j];
                  Ia_b_tmp += G_prev.locations[name1_index].visitorsIa_b[j];
                  I_b_tmp += G_prev.locations[name1_index].visitorsI_b[j];

                  S_probs[j] = G_prev.locations[name1_index].visitorsS[j];
                  I_probs[j] = G_prev.locations[name1_index].visitorsI[j];
                  Ia_probs[j] = G_prev.locations[name1_index].visitorsIa[j];
                  E1_probs[j] = G_prev.locations[name1_index].visitorsE1[j];
                  E2_probs[j] = G_prev.locations[name1_index].visitorsE2[j];
                  I_b_probs[j] = G_prev.locations[name1_index].visitorsI_b[j];
                  Ia_b_probs[j] = G_prev.locations[name1_index].visitorsIa_b[j];
                  E1_b_probs[j] = G_prev.locations[name1_index].visitorsE1_b[j];
                  E2_b_probs[j] = G_prev.locations[name1_index].visitorsE2_b[j];
                  R_probs[j] = G_prev.locations[name1_index].visitorsR[j];
                }

                p[0] = S_tmp;
                p[1] = E1_tmp;
                p[2] = E2_tmp;
                p[3] = I_tmp;
                p[4] = Ia_tmp;
                p[5] = R_tmp;
                p[6] = E1_b_tmp;
                p[7] = E2_b_tmp;
                p[8] = I_b_tmp;
                p[9] = Ia_b_tmp;

                // Draw the number of travellers from each compartment
                if(S_tmp + E1_tmp + E2_tmp + I_tmp + Ia_tmp + R_tmp + E1_b_tmp + E2_b_tmp + I_b_tmp + Ia_b_tmp > leftover){
                  rng_mvhyper2strains(p, S_tmp + E1_tmp + E2_tmp + I_tmp + Ia_tmp + R_tmp + E1_b_tmp + E2_b_tmp + I_b_tmp + Ia_b_tmp, leftover, &x);
                }
                else{
                  x = p;
                }

                //Distribute the travellers
                double *probs = new double[num];
                double *probs_cum = new double[num];
                double randomnumber;
                int index = -1;
                // Distribute susceptible
                for(int i = 0; i < x[0]; ++i){
                  for(int k = 0; k < num; ++k) probs[k] = S_probs[k]*1.0/ S_tmp; // Oneline for-loop
                  randomnumber = rand()*1.0/RAND_MAX;
                  while(randomnumber == 0 || randomnumber == 1){
                    randomnumber = rand()*1.0/RAND_MAX;
                  }
                  cumulative_sum2strains(&probs_cum, &probs, num);
                  for(int h = 0; h < num; ++h){
                    if(randomnumber < probs_cum[h]){
                      index = h;
                      break;
                    }
                  }

                  G_current.locations[name1_index].visitorsS[index] -= 1;
                  G_prev.locations[name1_index].visitorsS[index] -= 1;
                  G_current.locations[name2_index].visitorsS[index] += 1;

                  S_probs[index] -= 1;
                  S_tmp -= 1;
                }

                //Distribute exposed 1
                for(int i = 0; i < x[1]; ++i){
                  for(int k = 0; k < num; ++k) probs[k] = E1_probs[k]*1.0/ E1_tmp; // Oneline for-loop
                  randomnumber = rand()*1.0/RAND_MAX;
                  while(randomnumber == 0 || randomnumber == 1){
                    randomnumber = rand()*1.0/RAND_MAX;
                  }
                  cumulative_sum2strains(&probs_cum, &probs, num);
                  for(int h = 0; h < num; ++h){
                    if(randomnumber < probs_cum[h]){
                      index = h;
                      break;
                    }
                  }

                  G_current.locations[name1_index].visitorsE1[index] -= 1;
                  G_prev.locations[name1_index].visitorsE1[index] -= 1;
                  G_current.locations[name2_index].visitorsE1[index] += 1;

                  E1_probs[index] -= 1;
                  E1_tmp -= 1;
                }

                //Distribute exposed 2
                for(int i = 0; i < x[2]; ++i){
                  for(int k = 0; k < num; ++k) probs[k] = E2_probs[k]*1.0/ E2_tmp; // Oneline for-loop
                  randomnumber = rand()*1.0/RAND_MAX;
                  while(randomnumber == 0 || randomnumber == 1){
                    randomnumber = rand()*1.0/RAND_MAX;
                  }
                  cumulative_sum2strains(&probs_cum, &probs, num);
                  for(int h = 0; h < num; ++h){
                    if(randomnumber < probs_cum[h]){
                      index = h;
                      break;
                    }
                  }

                  G_current.locations[name1_index].visitorsE2[index] -= 1;
                  G_prev.locations[name1_index].visitorsE2[index] -= 1;
                  G_current.locations[name2_index].visitorsE2[index] += 1;

                  E2_probs[index] -= 1;
                  E2_tmp -= 1;
                }

                // Distribute infectious
                for(int i = 0; i < x[3]; ++i){
                  for(int k = 0; k < num; ++k) probs[k] = I_probs[k] * 1.0/ I_tmp; // Oneline for-loop
                  randomnumber = rand()*1.0/RAND_MAX;
                  while(randomnumber == 0 || randomnumber == 1){
                    randomnumber = rand()*1.0/RAND_MAX;
                  }
                  cumulative_sum2strains(&probs_cum, &probs, num);
                  for(int h = 0; h < num; ++h){
                    if(randomnumber < probs_cum[h]){
                      index = h;
                      break;
                    }
                  }

                  G_current.locations[name1_index].visitorsI[index] -= 1;
                  G_prev.locations[name1_index].visitorsI[index] -= 1;
                  G_current.locations[name2_index].visitorsI[index] += 1;

                  I_probs[index] -= 1;
                  I_tmp -= 1;
                }

                // Distribute asymptomatic infectious
                for(int i = 0; i < x[4]; ++i){
                  for(int k = 0; k < num; ++k) probs[k] = Ia_probs[k] * 1.0/ Ia_tmp; // Oneline for-loop
                  randomnumber = rand()*1.0/RAND_MAX;
                  while(randomnumber == 0 || randomnumber == 1){
                    randomnumber = rand()*1.0/RAND_MAX;
                  }
                  cumulative_sum2strains(&probs_cum, &probs, num);
                  for(int h = 0; h < num; ++h){
                    if(randomnumber < probs_cum[h]){
                      index = h;
                      break;
                    }
                  }

                  G_current.locations[name1_index].visitorsIa[index] -= 1;
                  G_prev.locations[name1_index].visitorsIa[index] -= 1;
                  G_current.locations[name2_index].visitorsIa[index] += 1;

                  Ia_probs[index] -= 1;
                  Ia_tmp -= 1;
                }

                // Distribute recovered
                for(int i = 0; i < x[5]; ++i){
                  for(int k = 0; k < num; ++k) probs[k] = R_probs[k] * 1.0/ R_tmp; // Oneline for-loop
                  randomnumber = rand()*1.0/RAND_MAX;
                  while(randomnumber == 0 || randomnumber == 1){
                    randomnumber = rand()*1.0/RAND_MAX;
                  }
                  cumulative_sum2strains(&probs_cum, &probs, num);
                  for(int h = 0; h < num; ++h){
                    if(randomnumber < probs_cum[h]){
                      index = h;
                      break;
                    }
                  }

                  G_current.locations[name1_index].visitorsR[index] -= 1;
                  G_prev.locations[name1_index].visitorsR[index] -= 1;
                  G_current.locations[name2_index].visitorsR[index] += 1;

                  R_probs[index] -= 1;
                  R_tmp -= 1;
                }

                //Distribute exposed 1 b
                for(int i = 0; i < x[6]; ++i){
                  for(int k = 0; k < num; ++k) probs[k] = E1_b_probs[k]*1.0/ E1_b_tmp; // Oneline for-loop
                  randomnumber = rand()*1.0/RAND_MAX;
                  while(randomnumber == 0 || randomnumber == 1){
                    randomnumber = rand()*1.0/RAND_MAX;
                  }
                  cumulative_sum2strains(&probs_cum, &probs, num);
                  for(int h = 0; h < num; ++h){
                    if(randomnumber < probs_cum[h]){
                      index = h;
                      break;
                    }
                  }

                  G_current.locations[name1_index].visitorsE1_b[index] -= 1;
                  G_prev.locations[name1_index].visitorsE1_b[index] -= 1;
                  G_current.locations[name2_index].visitorsE1_b[index] += 1;

                  E1_b_probs[index] -= 1;
                  E1_b_tmp -= 1;
                }

                //Distribute exposed 2 b
                for(int i = 0; i < x[7]; ++i){
                  for(int k = 0; k < num; ++k) probs[k] = E2_b_probs[k]*1.0/ E2_b_tmp; // Oneline for-loop
                  randomnumber = rand()*1.0/RAND_MAX;
                  while(randomnumber == 0 || randomnumber == 1){
                    randomnumber = rand()*1.0/RAND_MAX;
                  }
                  cumulative_sum2strains(&probs_cum, &probs, num);
                  for(int h = 0; h < num; ++h){
                    if(randomnumber < probs_cum[h]){
                      index = h;
                      break;
                    }
                  }

                  G_current.locations[name1_index].visitorsE2_b[index] -= 1;
                  G_prev.locations[name1_index].visitorsE2_b[index] -= 1;
                  G_current.locations[name2_index].visitorsE2_b[index] += 1;

                  E2_b_probs[index] -= 1;
                  E2_b_tmp -= 1;
                }

                // Distribute infectious b
                for(int i = 0; i < x[8]; ++i){
                  for(int k = 0; k < num; ++k) probs[k] = I_b_probs[k] * 1.0/ I_b_tmp; // Oneline for-loop
                  randomnumber = rand()*1.0/RAND_MAX;
                  while(randomnumber == 0 || randomnumber == 1){
                    randomnumber = rand()*1.0/RAND_MAX;
                  }
                  cumulative_sum2strains(&probs_cum, &probs, num);
                  for(int h = 0; h < num; ++h){
                    if(randomnumber < probs_cum[h]){
                      index = h;
                      break;
                    }
                  }

                  G_current.locations[name1_index].visitorsI_b[index] -= 1;
                  G_prev.locations[name1_index].visitorsI_b[index] -= 1;
                  G_current.locations[name2_index].visitorsI_b[index] += 1;

                  I_b_probs[index] -= 1;
                  I_b_tmp -= 1;
                }

                // Distribute asymptomatic infectious b
                for(int i = 0; i < x[9]; ++i){
                  for(int k = 0; k < num; ++k) probs[k] = Ia_b_probs[k] * 1.0/ Ia_b_tmp; // Oneline for-loop
                  randomnumber = rand()*1.0/RAND_MAX;
                  while(randomnumber == 0 || randomnumber == 1){
                    randomnumber = rand()*1.0/RAND_MAX;
                  }
                  cumulative_sum2strains(&probs_cum, &probs, num);
                  for(int h = 0; h < num; ++h){
                    if(randomnumber < probs_cum[h]){
                      index = h;
                      break;
                    }
                  }

                  G_current.locations[name1_index].visitorsIa_b[index] -= 1;
                  G_prev.locations[name1_index].visitorsIa_b[index] -= 1;
                  G_current.locations[name2_index].visitorsIa_b[index] += 1;

                  Ia_b_probs[index] -= 1;
                  Ia_b_tmp -= 1;
                }

                delete[] probs;
                delete[] probs_cum;

                delete[] S_probs;
                delete[] E1_probs;
                delete[] E2_probs;
                delete[] I_probs;
                delete[] Ia_probs;
                delete[] R_probs;
                delete[] E1_b_probs;
                delete[] E2_b_probs;
                delete[] I_b_probs;
                delete[] Ia_b_probs;


                if(S_tmp + E1_tmp + E2_tmp + I_tmp + Ia_tmp + R_tmp + E1_b_tmp + E2_b_tmp + I_b_tmp + Ia_b_tmp < leftover){
                  if( G_current.locations[name2_index].visitorsS[name1_index]  +
                      G_current.locations[name2_index].visitorsE1[name1_index]  +
                      G_current.locations[name2_index].visitorsE2[name1_index]  +
                      G_current.locations[name2_index].visitorsI[name1_index]  +
                      G_current.locations[name2_index].visitorsIa[name1_index] +
                      G_current.locations[name2_index].visitorsR[name1_index] +
                      G_current.locations[name2_index].visitorsE1_b[name1_index]  +
                      G_current.locations[name2_index].visitorsE2_b[name1_index]  +
                      G_current.locations[name2_index].visitorsI_b[name1_index]  +
                      G_current.locations[name2_index].visitorsIa_b[name1_index] < Nk)
                  {
                    G_current.locations[name2_index].visitorsS[name1_index] += leftover - S_tmp - E1_tmp - E2_tmp - I_tmp - Ia_tmp - R_tmp - E1_b_tmp - E2_b_tmp - I_b_tmp - Ia_b_tmp;
                    Rcout << "---------------" << endl;
                    Rcout << " Added extra people in " << tmpstr << " index " << name1_index << endl;
                    Rcout << " Name1 index " << name1_index << " Name2 index " << name2_index << endl;
                    Rcout << " Number added " << leftover - S_tmp - E1_tmp - E2_tmp - I_tmp - Ia_tmp - R_tmp - E1_b_tmp - E2_b_tmp - I_b_tmp - Ia_b_tmp << endl;
                  }
                }
              }
            }
            else if(Nk_prev > Nk){
              // Fewer people are going back, hence we draw them randomly
              p[0] = G_prev.locations[name1_index].visitorsS[name2_index];
              p[1] = G_prev.locations[name1_index].visitorsE1[name2_index];
              p[2] = G_prev.locations[name1_index].visitorsE2[name2_index];
              p[3] = G_prev.locations[name1_index].visitorsI[name2_index];
              p[4] = G_prev.locations[name1_index].visitorsIa[name2_index];
              p[5] = G_prev.locations[name1_index].visitorsR[name2_index];
              p[6] = G_prev.locations[name1_index].visitorsE1_b[name2_index];
              p[7] = G_prev.locations[name1_index].visitorsE2_b[name2_index];
              p[8] = G_prev.locations[name1_index].visitorsI_b[name2_index];
              p[9] = G_prev.locations[name1_index].visitorsIa_b[name2_index];
              sum = p[0] + p[1] + p[2] + p[3] + p[4] + p[5] + p[6] + p[7] + p[8] + p[9];
              // Draw the number of travellers from each compartment
              rng_mvhyper2strains(p, sum, Nk, &x);
              G_current.locations[name2_index].S += x[0];
              G_current.locations[name2_index].E1 += x[1];
              G_current.locations[name2_index].E2 += x[2];
              G_current.locations[name2_index].I += x[3];
              G_current.locations[name2_index].Ia += x[4];
              G_current.locations[name2_index].R += x[5];
              G_current.locations[name2_index].E1_b += x[6];
              G_current.locations[name2_index].E2_b += x[7];
              G_current.locations[name2_index].I_b += x[8];
              G_current.locations[name2_index].Ia_b += x[9];

              G_current.locations[name1_index].visitorsS[name2_index]  -= x[0];
              G_current.locations[name1_index].visitorsE1[name2_index]  -= x[1];
              G_current.locations[name1_index].visitorsE2[name2_index]  -= x[2];
              G_current.locations[name1_index].visitorsI[name2_index]  -= x[3];
              G_current.locations[name1_index].visitorsIa[name2_index] -= x[4];
              G_current.locations[name1_index].visitorsR[name2_index]  -= x[5];
              G_current.locations[name1_index].visitorsE1_b[name2_index]  -= x[6];
              G_current.locations[name1_index].visitorsE2_b[name2_index]  -= x[7];
              G_current.locations[name1_index].visitorsI_b[name2_index]  -= x[8];
              G_current.locations[name1_index].visitorsIa_b[name2_index] -= x[9];

              G_prev.locations[name1_index].visitorsS[name2_index]  -=  x[0];
              G_prev.locations[name1_index].visitorsE1[name2_index]  -=  x[1];
              G_prev.locations[name1_index].visitorsE2[name2_index]  -=  x[2];
              G_prev.locations[name1_index].visitorsI[name2_index]  -=  x[3];
              G_prev.locations[name1_index].visitorsIa[name2_index] -=  x[4];
              G_prev.locations[name1_index].visitorsR[name2_index]  -=  x[5];
              G_prev.locations[name1_index].visitorsE1_b[name2_index]  -=  x[6];
              G_prev.locations[name1_index].visitorsE2_b[name2_index]  -=  x[7];
              G_prev.locations[name1_index].visitorsI_b[name2_index]  -=  x[8];
              G_prev.locations[name1_index].visitorsIa_b[name2_index] -=  x[9];
            }
          }
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
      pop = G_current.locations[i].S +  G_current.locations[i].E1 + G_current.locations[i].E2 +  G_current.locations[i].I +  G_current.locations[i].Ia +  G_current.locations[i].R +  G_current.locations[i].E1_b + G_current.locations[i].E2_b +  G_current.locations[i].I_b +  G_current.locations[i].Ia_b;

      int num = G.locations[i].visitorsS.size();
      for (int j = 0; j < num; ++j){
        pop += G.locations[j].visitorsS[i];
        pop += G.locations[j].visitorsE1[i];
        pop += G.locations[j].visitorsE2[i];
        pop += G.locations[j].visitorsI[i];
        pop += G.locations[j].visitorsIa[i];
        pop += G.locations[j].visitorsR[i];
        pop += G.locations[j].visitorsE1_b[i];
        pop += G.locations[j].visitorsE2_b[i];
        pop += G.locations[j].visitorsI_b[i];
        pop += G.locations[j].visitorsIa_b[i];
      }
      if(pop < 100){
        baseline = 1;
      }
      else{
        baseline = 0.01*pop;
      }
      for(int i_t = 0; i_t < M; ++i_t){
        if(I_this[i][i_t] > pmax){
          pday = i_t;
          pmax = I_this[i][i_t];
        }
        if(I_this[i][i_t] > baseline){
          count += 1;
          if (count > 6 * 4){
            startday = i_t;
            count = -20000;
          }
        }
      }
      peak_date[i][i_sim] = pday;
      peak_val[i][i_sim] = pmax;
      start_date[i][i_sim] = startday;
    }
  }

  /// Values has first index for position, second for time and third for value.
  /// Third index 0=S, 1=E1, 2 = E2, 3=I, 4=Ia, 5=R; 6=E1_b; 7 = E2_b; 8 = I_b; 9 = Ia_b for belonging to kommune
  /// 10 = S, 11 = E1, 12 = E2,  13 = I, 14 = Ia, 15 = R, 16 = E1_b, 17 = E2_b, 18 = I_b, 19 = Ia_b, for currently in kommune.
  /// 20 = symptomatic incidence occurring in a kommune of strain A;
  /// 21 = asymptomatic incidence occuring in a kommune of strain A;
  /// 22 = symptomatic incidence occurring in a kommune of strain B;
  /// 23 = asymptomatic incidence occuring in a kommune of strain B;
  /// 24 = symptomatic incidence occurring in a kommune total;
  /// 25 = asymptomatic incidence occuring in a kommune total;

  StringVector res_names(n*4*M);
  IntegerVector res_week(n*4*M);
  IntegerVector res_day(n*4*M);
  IntegerVector res_time(n*4*M);
  IntegerVector res_B_S(n*4*M);
  IntegerVector res_B_E1(n*4*M);
  IntegerVector res_B_E2(n*4*M);
  IntegerVector res_B_I(n*4*M);
  IntegerVector res_B_Ia(n*4*M);
  IntegerVector res_B_R(n*4*M);
  IntegerVector res_B_E1_b(n*4*M);
  IntegerVector res_B_E2_b(n*4*M);
  IntegerVector res_B_I_b(n*4*M);
  IntegerVector res_B_Ia_b(n*4*M);
  IntegerVector res_C_S(n*4*M);
  IntegerVector res_C_E1(n*4*M);
  IntegerVector res_C_E2(n*4*M);
  IntegerVector res_C_I(n*4*M);
  IntegerVector res_C_Ia(n*4*M);
  IntegerVector res_C_R(n*4*M);
  IntegerVector res_C_E1_b(n*4*M);
  IntegerVector res_C_E2_b(n*4*M);
  IntegerVector res_C_I_b(n*4*M);
  IntegerVector res_C_Ia_b(n*4*M);
  IntegerVector res_INCIDENCE_a(n*4*M);
  IntegerVector res_as_incidence_a(n*4*M);
  IntegerVector res_INCIDENCE_b(n*4*M);
  IntegerVector res_as_incidence_b(n*4*M);
  IntegerVector res_INCIDENCE(n*4*M);
  IntegerVector res_as_incidence(n*4*M);

  int index=0;
  for(int i=0; i < n; ++i){
    for (int k = 0; k < 4 * M; ++k){
      res_names[index] = names[i];
      res_week[index] = k/(7*4)+1;
      res_day[index] = k/4+1;
      res_time[index] = (k%4)+1;

      res_B_S[index] = values[i][k][0]*1.0/N;
      res_B_E1[index] = values[i][k][1]*1.0/N;
      res_B_E2[index] = values[i][k][2]*1.0/N;
      res_B_I[index] = values[i][k][3]*1.0/N;
      res_B_Ia[index] = values[i][k][4]*1.0/N;
      res_B_R[index] = values[i][k][5]*1.0/N;
      res_B_E1_b[index] = values[i][k][6]*1.0/N;
      res_B_E2_b[index] = values[i][k][7]*1.0/N;
      res_B_I_b[index] = values[i][k][8]*1.0/N;
      res_B_Ia_b[index] = values[i][k][9]*1.0/N;

      res_C_S[index] = values[i][k][10]*1.0/N;
      res_C_E1[index] = values[i][k][11]*1.0/N;
      res_C_E2[index] = values[i][k][12]*1.0/N;
      res_C_I[index] = values[i][k][13]*1.0/N;
      res_C_Ia[index] = values[i][k][14]*1.0/N;
      res_C_R[index] = values[i][k][15]*1.0/N;
      res_C_E1_b[index] = values[i][k][16]*1.0/N;
      res_C_E2_b[index] = values[i][k][17]*1.0/N;
      res_C_I_b[index] = values[i][k][18]*1.0/N;
      res_C_Ia_b[index] = values[i][k][19]*1.0/N;

      res_INCIDENCE_a[index] = values[i][k][20]*1.0/N;
      res_as_incidence_a[index] = values[i][k][21]*1.0/N;

      res_INCIDENCE_b[index] = values[i][k][22]*1.0/N;
      res_as_incidence_b[index] = values[i][k][23]*1.0/N;

      res_INCIDENCE[index] = values[i][k][24]*1.0/N;
      res_as_incidence[index] = values[i][k][25]*1.0/N;

      index++;
    }
  }



  // return two data frames, as max 20 columns
  DataFrame df1 = DataFrame::create(
    _["location_code"]= res_names,
    _["week"]=res_week,
    _["day"]=res_day,
    _["time"]= res_time,
    _["b_S"]= res_B_S,
    _["b_E1"]= res_B_E1,
    _["b_E2"]= res_B_E2,
    _["b_I"]= res_B_I,
    _["b_Ia"]= res_B_Ia,
    _["b_R"]= res_B_R,
    _["b_E1_b"]= res_B_E1_b,
    _["b_E2_b"]= res_B_E2_b,
    _["b_I_b"]= res_B_I_b,
    _["b_Ia_b"]= res_B_Ia_b,
    _["c_S"]= res_C_S,
    _["c_E1"]= res_C_E1,
    _["c_E2"]= res_C_E2,
    _["c_I"]= res_C_I,
    _["c_Ia"]= res_C_Ia,
    _["c_R"]= res_C_R
  );


  // return a new data frame
  DataFrame df2 = DataFrame::create(
    _["c_E1_b"]= res_C_E1_b,
    _["c_E2_b"]= res_C_E2_b,
    _["c_I_b"]= res_C_I_b,
    _["c_Ia_b"]= res_C_Ia_b,
    _["c_symp_incidence_a"]= res_INCIDENCE_a,
    _["c_asymp_incidence_a"]= res_as_incidence_a,
    _["c_symp_incidence_b"]= res_INCIDENCE_b,
    _["c_asymp_incidence_b"]= res_as_incidence_b,
    _["c_symp_incidence"]= res_INCIDENCE,
    _["c_asymp_incidence"]= res_as_incidence
  );

  df1.attr("class") = Rcpp::CharacterVector::create("data.table", "data.frame");
  df2.attr("class") = Rcpp::CharacterVector::create("data.table", "data.frame");

  for(int i = 0; i < n; ++i){
    for(int k = 0; k < 4 * M; ++k){
      delete[] values[i][k];
    }
    delete[] values[i];
  }
  delete[] values;
  for(int i = 0; i < n; ++i){
    delete[] I_this[i];
    delete[] peak_date[i];
    delete[] peak_val[i];
    delete[] start_date[i];
    delete[] final_size[i];
  }
  for (int i = 0; i < N; ++i){
    delete[] bonds[i];
  }
  delete[] I_this;
  delete[] peak_date;
  delete[] start_date;
  delete[] peak_val;
  delete[] final_size;
  delete[] bonds;


  return(Rcpp::List::create(Rcpp::Named("df1") = df1, Rcpp::Named("df2") = df2));
}


/*** R
se1e2iiar_2strains_pop <- data.table::data.table(
  "location_code" = c("a","b","c"),
  "S" = c(1000,1000,2000),
  "E1" = c(0,0,0),
  "E2" = c(0,0,100),
  "I" = c(10,0,0),
  "Ia" = c(0,0,0),
  "R" = c(0,0,0),
  "E1_b" = c(0,0,0),
  "E2_b" = c(0,0,100),
  "I_b" = c(10,0,0),
  "Ia_b" = c(0,0,0)
)

temp <- data.table::data.table(
  from = c("a","a","b","b","c","c"),
  to = c("b","c","a","c","a","b"),
  n = c(10,10,10,10,10,3000)
)
mobility_matrix <- vector("list",length=10 * 4)
for(i in seq_along(mobility_matrix)){
  mobility_matrix[[i]] <- data.table::copy(temp)
  data.table::setnames(mobility_matrix[[i]],c("from","to","n"))
}

seed_matrix <- matrix(0, nrow = 10, ncol = 3)
seed_matrix_b <- matrix(0, nrow = 10, ncol = 3)

dayEach <- rep(1:10, each = 4)
timeEach <- rep(c(0, 6, 12, 18), 10)
location_codes <- c(rep("a", 10 * 4), rep("b", 10 * 4), rep("c", 10 * 4))
betas <- c(rep(0.6, 10 * 4), rep(0.4, 10 * 2), rep(0.2, 10 * 2), rep(0.7, 10 * 3), rep(0.4, 10))
days <- rep(dayEach, 3)
times <- rep(timeEach, 3)
betaDF <- data.table("location_code" = location_codes, "day" = days, "time" = times, "beta" = betas)
betaDF = betaDF[order(betaDF$day), ]
betas <- convert_beta_to_matrix(betaDF, location_codes = c("a", "b", "c"), days = 1:10, times = c(0, 6, 12, 18))

d <- asymmetric_mobility_se1e2iiar_2strains_cpp(
  se1e2iiar_2strains_pop = se1e2iiar_2strains_pop,
  mobility_matrix = mobility_matrix,
  seed_matrix = seed_matrix,
  seed_matrix_b = seed_matrix_b,
  betas=betas,
  inputSeed = 3,
  a1=1/2.0,
  a2 = 1/3.0,
  gamma= 1/5.0,
  relativeInfectiousnessB = 1.5,
  presymptomaticRelativeInfectiousness = 1.25,
  asymptomaticProb = 0.4,
  asymptomaticRelativeInfectiousness = 0.5,
  N=1,
  M=10
)
d
*/
