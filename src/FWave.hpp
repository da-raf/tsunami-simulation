#include <assert.h>

#ifndef FWAVE_H
#define FWAVE_H

#define G 9.81

namespace solver {

/**
 * This class will approximately solve the initial value problem for the <b>shallow water
 * equations</b> over time:
 * <pre>
 * / h  \   /    hu           \
 * |    | + |                 |  = S(x,t)
 * \ hu /   \ hu^2 + 1/2*gh^2 /
 * </pre>
 */
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
     * @param hNetUpdatesLeft reference to the variable that will store the update to the height of the water column on the left
     * @param hNetUpdatesRight reference to the variable that will store the update to the height of the water column on the right
     * @param huNetUpdatesLeft reference to the variable that will store the update to the momentum of the water on the left
     * @param huNetUpdatesRight reference to the variable that will store the update to the momentum of the water on the right
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
    /**
     * calculate the result of the flux function for this quantity
     * 
     * @param h height of the water column
     * @param hu momentum of the water current
     * @param f1 reference to the variable the result of the first component of the flux function will be written to
     * @param f2 reference to the variable the result of the second component of the flux function will be written to
     */
    void flux(T h, T hu, T& f1, T& f2);
    
    /** 
     * calculate the Roe eigenvalues for the left and right quantities
     * 
     * @param hl water column height on the left side
     * @param hr water column height on the right side
     * @param hul momentum of the water on the left side
     * @param hur momentum of the water on the right side
     * @param roe1 reference to the variable the first roe-eigenvalue will be written into
     * @param roe2 reference to the variable the second roe-eigenvalue will be written into
     * 
     */
    void roeEigenvals(T hl, T hr, T hul, T hur, T& roe1, T& roe2);
    
    /**
     * calculate the eigencoefficients for given states left and right and roe-eigenvalues
     * 
     * @param hl water column height on the left side
     * @param hr water column height on the right side
     * @param hul momentum of the water on the left side
     * @param hur momentum of the water on the right side
     * @param roe1 first roe-eigenvalue
     * @param roe2 second roe-eigenvalue
     * @param alpha1 reference to the variable that the first eigencoefficient will be written into
     * @param alpha2 reference to the variable that the second eigencoefficient will be written into
     * 
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

