// f-wave solver

#include <stdlib.h>
#include <math.h>

#define G 9.81

typedef struct {
    double h;
    double hu;
} State;


double* roe_eigenvals(State l, State r) {
    double* eigenvals = (double*) malloc(2 * sizeof(double));
    
    double u_l = l.hu / l.h;
    double u_r = r.hu / r.h;
    
    double sqrt_hl = sqrt(l.h);
    double sqrt_hr = sqrt(r.h);
    
    double h_roe = (l.h + r.h) / 2;
    double u_roe = (u_l * sqrt_hl + u_r * sqrt_hr) / (sqrt_hl + sqrt_hr);
    
    double c = sqrt( G * h_roe );
    
    eigenvals[0] = u_roe - c;
    eigenvals[1] = u_roe + c;
    
    return eigenvals;
}

