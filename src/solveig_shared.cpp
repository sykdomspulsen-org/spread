#include <Rcpp.h>

void cumulative_sum(double **outarray, double **array, int n){
  /// REMEMBER to delete outarray after use
  (*outarray)[0] = (*array)[0];
  for(int i=1; i < n; ++i){
    (*outarray)[i] = (*outarray)[i-1] + (*array)[i];
  }
}

void seir_sim(int &ds, int &de1, int &de2, int &dia, int &di,       // Outputs
              int S, int E, int Ia, int I, float beta, float a,  float gamma, float asymptomaticProb, float asymptomaticRelativeInfectiousness, int pop, double delta_t){ // Inputs
  ds = 0; de1 = 0; de2 = 0; dia = 0; di = 0;
  // Function to run one time step of the seir model
  // ds are susceptible going to exposed
  // de1 are the exposed going to asymptomatic
  // de2 are the exposed going to symptomatic
  // dia are the infectious asymptomatic going to recovered
  // di are the infectious symptomatic going to recovered
  int de = 0;
  if(I != 0 || Ia != 0 || E != 0){
    if( E == 0){
      de = 0;
    }
    else{
      de = R::rbinom(E, a*delta_t);
      if(de != 0){
        de1 = R::rbinom(de, asymptomaticProb);
        de2 = de-de1;
      }
    }
    if(I == 0){
      di = 0;
    }
    else{
      di = R::rbinom(I, gamma*delta_t);
    }
    if(Ia == 0){
      dia = 0;
    }
    else{
      dia = R::rbinom(Ia, gamma*delta_t);
    }
    ds = R::rbinom(S, beta*delta_t*I/pop + asymptomaticRelativeInfectiousness*beta*delta_t*Ia/pop);
  }
}
