/**
 * @file FWave.hpp f-wave solver
 * 
 * @author Raphael Dümig
 */

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
     * using the height of the water column and its momentum as parameters.
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
     */
    void computeNetUpdates(const T &hLeft, const T &hRight, const T &huLeft, const T &huRight,
                           const T &bathLeft, const T &bathRight,
                           T& hNetUpdatesLeft,  T& hNetUpdatesRight,
                           T& huNetUpdatesLeft, T& huNetUpdatesRight,
                           T& maxEdgeSpeed );
    
private:
    /**
     * calculate the result of the flux function for this quantity
     * 
     * @param h height of the water column
     * @param hu momentum of the water current
     * @param fl pointer to an array with at least two elements, where the resulting vector will be written into
     */
    void flux(const T &h, const T &hu, T *fl);
    
    /** 
     * calculate the Roe eigenvalues for the left and right quantities
     * 
     * @param hl water column height on the left side
     * @param hr water column height on the right side
     * @param hul momentum of the water on the left side
     * @param hur momentum of the water on the right side
     * @param lambda_roe an array with two elements, where the eigenvalues will be written into
     */
    void roeEigenvals(const T &hl, const T &hr, const T &hul, const T &hur, T *lambda_roe);
    
    /**
     * calculate the eigencoefficients for given states left and right and roe-eigenvalues
     * 
     * @param hl water column height on the left side
     * @param hr water column height on the right side
     * @param hul momentum of the water on the left side
     * @param hur momentum of the water on the right side
     * @param bl elevation of the ocean floor on the left side
     * @param br elevation of the ocean floor on the right side
     * @param lambda_roe an array containing the two roe eigenvalues
     * @param alpha an array with two elements, where the eigencoefficients will be written into
     */
    void eigencoeffis(const T &hl, const T &hr, const T &hul, const T &hur, const T &bl, const T &br,
                      const T *lambda_roe, T *alpha );
};

// end of namespace solver
};

#include "FWave.cpp"

#endif

