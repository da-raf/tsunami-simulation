#include <assert.h>

#ifndef FWAVE_H
#define FWAVE_H

#define G 9.81

namespace solver {

template<class T>
class FWave
{
public:
	void computeNetUpdates(
		T hLeft,             T hRight,
		T huLeft,            T huRight,
		T bathLeft,          T bathRight,
		T& hNetUpdatesLeft,  T& hNetUpdatesRight,
		T& huNetUpdatesLeft, T& huNetUpdatesRight,
		T& maxEdgeSpeed );
private:
    /** calculate the result of the flux function for this quantity
     * @param q quantity
     */
    void flux(T h, T hu, T& f1, T& f2);
    
    /** calculate the Roe eigenvalues for the left and right quantities
     * @param ql left quantity
     * @param qr right quantity
     * @result Roe eigenvalues (x: \lambda_l, y: lambda_r)
     */
    void roeEigenvals(T hl, T hr, T hul, T hur, T& roe1, T& roe2);
    
    /** calculate the eigencoefficients
     */
    void eigencoeffis(
        T hl, T hr, T hul, T hur,
        T roe1, T roe2,
        T& alpha1, T& alpha2 );
};

#include "FWave.cpp"
// end of namespace solver
}

#endif

