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


Vector2D roe_eigenvals(State ql, State qr) {
    Vector2D eigenvals;
    
    double u_l = ql.hu / ql.h;
    double u_r = qr.hu / qr.h;
    
    double sqrt_hl = sqrt(ql.h);
    double sqrt_hr = sqrt(qr.h);
    
    double h_roe = (ql.h + qr.h) / 2;
    double u_roe = (u_l * sqrt_hl + u_r * sqrt_hr) / (sqrt_hl + sqrt_hr);
    
    double c = sqrt( G * h_roe );
    
    eigenvals.x = u_roe - c;
    eigenvals.y = u_roe + c;
    
    return eigenvals;
}

Vector2D eigencoeffis(State ql, State qr, Vector2D roe_eigenvals) {
    Vector2D f1 = flux(ql);
    Vector2D f2 = flux(qr);
    
    Vector2D df;
    df.x = f1.x - f2.x;
    df.y = f1.y - f2.y;
    
    // matrix-vector multiplication:
    // /  lambda_2 -1 \-1
    // \ -lambda_1  1 /    * df
    
    assert(roe_eigenvals.y - roe_eigenvals.x != 0);
    double c = 1 / (roe_eigenvals.y - roe_eigenvals.x);
    
    Vector2D result;
    result.x = (roe_eigenvals.y * df.x - df.y) * c;
    result.y = (df.y - roe_eigenvals.x * df.x) * c;
    
    return result;
}
