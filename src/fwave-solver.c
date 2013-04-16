// f-wave solver

#include <stdlib.h>
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

Vector2D eigencoeffis(Vector2D roe_eigenvals, Vector2D df) {
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


Vector2D* calculate_updates(State ql, State qr) {
    Vector2D lambda_roe = roe_eigenvals(ql, qr);
    
    if(lambda_roe.x < 0.0 && lambda_roe.y < 0.0)
        lambda_roe.y = 0;
    else if(lambda_roe > 0.0 && lambda_roe > 0.0)
        lambda_roe.x = 0;
    
    Vector2D fl = flux(ql);
    Vector2D fr = flux(qr);
    
    Vector2D df;
    df.x = fr.x - fl.x;
    df.y = fr.y - fl.y;
    
    Vector2D alpha = eigencoeffis(lambda_roe, df);
    
    Vector2D[2] z;
    
    z[0].x = alpha.x;
    z[0].y = alpha.x * lambda_roe.x;
    
    z[1].x = alpha.y;
    z[1].y = alpha.y * lambda_roe.y;
    
    
    Vector2D* result = (Vector2D*) malloc(2 * sizeof(Vector2D));
    // A- dQ
    if(alpha.x < 0) {
        result[0].x = z[0].x;
        result[0].y = z[0].y;
    }
    else {
        result[0].x = 0;
        result[0].y = 0;
    }
    
    if(alpha.y < 0) {
        result[0].x += z[1].x;
        result[0].y += z[1].y;
    }
    
    
    // A+ dQ
    if(alpha.x > 0) {
        result[1].x = z[0].x;
        result[1].y = z[0].y;
    }
    else {
        result[1].x = 0;
        result[1].y = 0;
    }
    
    if(alpha.y > 0) {
        result[1].x += z[1].x;
        result[1].y += z[1].y;
    }
    
    return result;
}
