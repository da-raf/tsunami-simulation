#include "FWave.hpp"

#include <math.h>
#include <assert.h>

template<class T>
void FWave<T>::flux(const T &h, const T &hu, T *fl)
{
    fl[0] = hu;
    fl[1] = (hu * hu) / h + G * (h * h) / 2;
    
    return;
}

template<class T>
void FWave<T>::roeEigenvals(
        const T &hl,  const T &hr,
        const T &hul, const T &hur,
        T *lambda_roe)
{
    // assert that we got positiv values for the water columns
    assert(hl > 0.0);
    assert(hr > 0.0);

    // calculate the particle speeds
    T ul = hul / hl;
    T ur = hur / hr;
    
    T sqrt_hl = sqrt(hl);
    T sqrt_hr = sqrt(hr);
    
    T h_roe = (hl + hr) / 2;
    T u_roe = (ul * sqrt_hl + ur * sqrt_hr) / (sqrt_hl + sqrt_hr);
    
    T c = sqrt( G * h_roe );
    
    // calculate the roe eigenvalues
    lambda_roe[0] = u_roe - c;
    lambda_roe[1] = u_roe + c;
    
    return;
}

template<class T>
void FWave<T>::eigencoeffis(
        const T &hl,  const T &hr,
        const T &hul, const T &hur,
        const T *lambda_roe,
        T *alpha )
{
    // flux left and right
    T *fl = new T[2];
    T *fr = new T[2];
    
    // the difference of the two fluxes
    T *df = new T[2];
    
    // calculate the flux
    flux(hr, hur, fr);
    flux(hl, hul, fl);
    
    // calculate the difference
    df[0] = fr[0] - fl[0];
    df[1] = fr[1] - fl[1];
    
    
    // matrix-vector multiplication:
    // /        1 |        1 \-1
    // \ lambda_1 | lambda_2 /    * df
    
    // make sure we get a numerically stable result
    assert(fabs(lambda_roe[1] - lambda_roe[0]) > 0.01);
    
    // perform the matrix-vector multiplication
    T c = 1 / (lambda_roe[1] - lambda_roe[0]);
    
    alpha[0] = (lambda_roe[1] * df[0] - df[1]) * c;
    alpha[1] = (df[1] - lambda_roe[0] * df[0]) * c;
    
    return;
}

template<class T>
void FWave<T>::computeNetUpdates(
		const T &hLeft,     const T &hRight,
		const T &huLeft,    const T &huRight,
		const T &bathLeft,  const T &bathRight,
		T &hNetUpdateLeft,  T &hNetUpdateRight,
		T &huNetUpdateLeft, T &huNetUpdateRight,
		T &maxEdgeSpeed )
{
    // compute the eigenvalues
    T *lambda_roe = new T[2];
    roeEigenvals(hLeft, hRight, huLeft, huRight, lambda_roe); // formula 3 and 4
    
    // compute the eigencoefficients
    T *alpha = new T[2];
    eigencoeffis(hLeft, hRight, huLeft, huRight, lambda_roe, alpha); // formula 8
    
    // compute the waves
    T h[2];
    T hu[2];
    
    h[0]  = alpha[0];
    hu[0] = alpha[0] * lambda_roe[0];
    
    h[1]  = alpha[1];
    hu[1] = alpha[1] * lambda_roe[1];
    
    // compute the net-updates
    hNetUpdateLeft   = 0.0;
    huNetUpdateLeft  = 0.0;
    hNetUpdateRight  = 0.0;
    huNetUpdateRight = 0.0;
    
    // A- dQ
    int i;
    for(i=0; i<2; i++) {
        if(lambda_roe[i] < 0.0) {
            hNetUpdateLeft  += h[i];
            huNetUpdateLeft += hu[i];
        }
    }
    
    // A+ dQ
    for(i=0; i<2; i++) {
        if(lambda_roe[i] > 0.0) {
            hNetUpdateRight  += h[i];
            huNetUpdateRight += hu[i];
        }
    }
    
    // calculate the maximum speed
    T lambda_l = lambda_roe[0];
    T lambda_r = lambda_roe[1];
    
    if(lambda_roe[0] < 0.0 && lambda_roe[1] < 0.0)
        lambda_r = 0.0;
    else if(lambda_roe[0] > 0.0 && lambda_roe[1] > 0.0)
        lambda_l = 0.0;
    
    lambda_l = fabs(lambda_l);
    lambda_r = fabs(lambda_r);

    if(lambda_l > lambda_r)
        maxEdgeSpeed = lambda_l;
    else
        maxEdgeSpeed = lambda_r;
    
    return;
}

