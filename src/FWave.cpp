/**
 * @file FWave.cpp Implementation of an f-wave solver
 * 
 * @author Raphael Dümig
 */

#include "FWave.hpp"

#include <cmath>
//#include <math.h>
#include <cassert>

using namespace solver;

template<class T>
void FWave<T>::flux(const T &h, const T &hu,const T &u, T *fl)
{
    fl[0] = hu;
    fl[1] = (hu * u) + G * (h * h) / 2;
    
    return;
}

template<class T>
void FWave<T>::roeEigenvals(
        const T &hl,  const T &hr,
        const T &hul, const T &hur,
        const T &ul, const T &ur,
        T *lambda_roe)
{
    // assert that we got positiv values for the water columns
    assert(hl > 0.0);
    assert(hr > 0.0);

    // calculate the particle speeds
    
    T sqrt_hl = std::sqrt(hl);
    T sqrt_hr = std::sqrt(hr);
    
    T h_roe = (hl + hr) / 2;
    T u_roe = (ul * sqrt_hl + ur * sqrt_hr) / (sqrt_hl + sqrt_hr);
    
    T c = std::sqrt( G * h_roe );
    
    // calculate the roe eigenvalues
    lambda_roe[0] = u_roe - c;
    lambda_roe[1] = u_roe + c;
    
    return;
}

template<class T>
void FWave<T>::eigencoeffis(
        const T &hl,  const T &hr,
        const T &hul, const T &hur,
        const T &bl,  const T &br,
        const T &ul, const T &ur,
        const T *lambda_roe,
        T *alpha )
{
    // flux left and right
    T fl[2];
    T fr[2];
    
    // the difference of the two fluxes
    T df[2];
    
    // calculate the flux
    flux(hr, hur, ur, fr);
    flux(hl, hul, ul, fl);
    
    T bathymetry_influence = G * (bl - br) * (hl + hr) / 2;
    
    // calculate the difference
    df[0] = fr[0] - fl[0];
    df[1] = fr[1] - fl[1] - bathymetry_influence;
    
    
    // matrix-vector multiplication:
    // /        1 |        1 \-1
    // \ lambda_1 | lambda_2 /    * df
    
    // make sure we get a numerically stable result
    assert( std::abs(lambda_roe[1] - lambda_roe[0]) > 0.01);
    
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
        const T &uLeft, const T &uRight,
		T &hNetUpdateLeft,  T &hNetUpdateRight,
		T &huNetUpdateLeft, T &huNetUpdateRight,
		T &maxEdgeSpeed )
{
    
    T lambda_roe[2];
    T alpha[2];
    
    if(bathLeft < 0.0 && bathRight < 0.0) {
        // normal scenario: water on both sides
        
        // compute the eigenvalues
        roeEigenvals(hLeft, hRight, huLeft, huRight, uLeft, uRight, lambda_roe); // formula 3 and 4
        
        // compute the eigencoefficients
        eigencoeffis(hLeft, hRight, huLeft, huRight, bathLeft, bathRight, uLeft, uRight, lambda_roe, alpha); // formula 8
    }
    else if(bathRight < 0.0) {
        // the cell on the left is dry => reflection to the right
        roeEigenvals(hRight, hRight, (-1) * huRight, huRight, uLeft, uRight, lambda_roe);
        eigencoeffis(hRight, hRight, (-1) * huRight, huRight, bathRight, bathRight, uLeft, uRight, lambda_roe, alpha);
    }
    else if(bathLeft < 0.0) {
        // the cell on the right is dry => reflection to the left
        roeEigenvals(hLeft, hLeft, huLeft, (-1) * huLeft, uLeft, uRight, lambda_roe);
        eigencoeffis(hLeft, hLeft, huLeft, (-1) * huLeft, bathLeft, bathLeft, uLeft, uRight, lambda_roe, alpha);
    }
    
    
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
    
    if(bathLeft < 0.0)
    {
        // A- dQ
        int i;
        for(i=0; i<2; i++) {
            if(lambda_roe[i] < 0.0) {
                hNetUpdateLeft  += h[i];
                huNetUpdateLeft += hu[i];
            }
        }
    }
    
    if(bathRight < 0.0)
    {
        // A+ dQ
        int j;
        for(j=0; j<2; j++) {
            if(lambda_roe[j] > 0.0) {
                hNetUpdateRight  += h[j];
                huNetUpdateRight += hu[j];
            }
        }
    }
    
    // calculate the maximum speed
    T lambda_l = lambda_roe[0];
    T lambda_r = lambda_roe[1];
    
    if(lambda_roe[0] < 0.0 && lambda_roe[1] < 0.0)
        lambda_r = 0.0;
    else if(lambda_roe[0] > 0.0 && lambda_roe[1] > 0.0)
        lambda_l = 0.0;
    
    lambda_l = std::abs(lambda_l);
    lambda_r = std::abs(lambda_r);

    if(lambda_l > lambda_r)
        maxEdgeSpeed = lambda_l;
    else
        maxEdgeSpeed = lambda_r;
    
    return;
}

