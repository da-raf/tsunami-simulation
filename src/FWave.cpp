#include "FWave.hpp"

#include <math.h>
#include <assert.h>

template<class T>
void FWave<T>::flux(T h, T hu, T& f1, T& f2)
{
    T u = hu / h;
    f1  = hu;
    f2  = h * (u * u) + G * (h * h) / 2;
    
    return;
}

template<class T>
void FWave<T>::roeEigenvals(
        T hl, T hr, T hul, T hur,
        T& roe1, T& roe2)
{
    assert(hl > 0.0);
    assert(hr > 0.0);

    T ul = hul / hl;
    T ur = hur / hr;
    
    T sqrt_hl = sqrt(hl);
    T sqrt_hr = sqrt(hr);
    
    T h_roe = (hl + hr) / 2;
    T u_roe = (ul * sqrt_hl + ur * sqrt_hr) / (sqrt_hl + sqrt_hr);
    
    T c = sqrt( G * h_roe );
    
    roe1 = u_roe - c;
    roe2 = u_roe + c;
    
    return;
}

template<class T>
void FWave<T>::eigencoeffis(
        T hl, T hr, T hul, T hur,
        T roe1, T roe2,
        T& alpha1, T& alpha2 )
{
    T f1l, f2l;
    T f1r, f2r;
    T df1, df2;
    
    flux(hr, hur, f1r, f2r);
    flux(hl, hul, f1l, f2l);
    
    df1 = f1r - f1l;
    df2 = f2r - f2l;
    
    // matrix-vector multiplication:
    // /        1 |        1 \-1
    // \ lambda_1 | lambda_2 /    * df
    
    assert(fabs(roe2 - roe1) > 0.01);
    
    T c = 1 / (roe2 - roe1);
    
    alpha1 = (roe2 * df1 - df2) * c;
    alpha2 = (df2 - roe1 * df1) * c;
    
    return;
}

template<class T>
void FWave<T>::computeNetUpdates(
		T hLeft,            T hRight,
		T huLeft,           T huRight,
		T bathLeft,         T bathRight,
		T& hNetUpdateLeft,  T& hNetUpdateRight,
		T& huNetUpdateLeft, T& huNetUpdateRight,
		T& maxEdgeSpeed )
{
    T roe[2];
    roeEigenvals(hLeft, hRight, huLeft, huRight, roe[0], roe[1]); // formel 3,4
    
    T alpha[2];
    eigencoeffis(hLeft, hRight, huLeft, huRight, roe[0], roe[1], alpha[0], alpha[1]); // Formel 8
    
    T h[2];
    T hu[2];
    
    h[0]  = alpha[0];
    hu[0] = alpha[0] * roe[0];
    
    h[1]  = alpha[1];
    hu[1] = alpha[1] * roe[1];
    
    
    hNetUpdateLeft   = 0.0;
    huNetUpdateLeft  = 0.0;
    hNetUpdateRight  = 0.0;
    huNetUpdateRight = 0.0;
    
    // A- dQ
    int i;
    for(i=0; i<2; i++) {
        if(roe[i] < 0.0) {
            hNetUpdateLeft  += h[i];
            huNetUpdateLeft += hu[i];
        }
    }
    
    // A+ dQ
    for(i=0; i<2; i++) {
        if(roe[i] > 0.0) {
            hNetUpdateRight  += h[i];
            huNetUpdateRight += hu[i];
        }
    }
    
    // calculate the maximum speed
    T lambda_l = roe[0];
    T lambda_r = roe[1];
    
    if(roe[0] < 0.0 && roe[1] < 0.0)
        lambda_r = 0.0;
    else if(roe[0] > 0.0 && roe[1] > 0.0)
        lambda_l = 0.0;
    
    lambda_l = fabs(lambda_l);
    lambda_r = fabs(lambda_r);

    if(lambda_l > lambda_r)
        maxEdgeSpeed = lambda_l;
    else
        maxEdgeSpeed = lambda_r;
    
    return;
}

