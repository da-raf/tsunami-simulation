#include <assert.h>

#ifndef FWAVE_H
#define FWAVE_H

#define G 9.81

namespace solver {

template<class T>
class FWave
{
public:
    /**
     * calculate the net-updates for a simulation of the flow of water
     * 
     * This implementation will calculate the net-updates for a simulation of flow of water
     * using the height of the water column and its momentum as parameters and working with
     * the system of finite differences.
     * 
     * @param hLeft water column height on the left side
     * @param hRight water column height on the right side
     * @param huLeft momentum of the water on the left side
     * @param huRight momentum of the water on the right side
     * @param bathLeft elevation of the ocean floor on the left side
     * @param bathRight elevation of the ocean floor on the right side
     * @param hNetUpdateLeft reference to the variable that will store the update to the height of the water column on the left
     * @param hNetUpdateRight reference to the variable that will store the update to the height of the water column on the right
     * @param huNetUpdateLeft reference to the variable that will store the update to the momentum of the water on the left
     * @param huNetUpdateRight reference to the variable that will store the update to the momentum of the water on the right
     * @param maxEdgeSpeed reference to the variable that will store the maximum edge-speed
     * 
     */
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

