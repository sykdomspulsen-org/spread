
#ifndef __SOLVEIG__
#define __SOLVEIG__

void cumulative_sum(double **outarray, double **array, int n);
void seir_sim(
    int &ds, int &de1, int &de2, int &dia, int &di,       // Outputs
                int S, int E, int Ia, int I, float beta, float a,  float gamma, float asymptomaticProb, float asymptomaticRelativeInfectiousness, int pop, double delta_t);

#endif
