// f-wave solver

#include <math.h>

#define G 9.81


typedef struct {
    double x;
    double y;
} Vector2D;


typedef struct {
    double h;
    double hu;
} State;


Vector2D flux(State q) {
    Vector2D result;
    
    result.x = q.hu;
    result.y = (q.hu * q.hu) + G * (q.h * q.h) / 2;
    
    return result;
}


Vector2D roe_eigenvals(State l, State r) {
    Vector2D eigenvals;
    
    double u_l = l.hu / l.h;
    double u_r = r.hu / r.h;
    
    double sqrt_hl = sqrt(l.h);
    double sqrt_hr = sqrt(r.h);
    
    double h_roe = (l.h + r.h) / 2;
    double u_roe = (u_l * sqrt_hl + u_r * sqrt_hr) / (sqrt_hl + sqrt_hr);
    
    double c = sqrt( G * h_roe );
    
    eigenvals.x = u_roe - c;
    eigenvals.y = u_roe + c;
    
    return eigenvals;
}

