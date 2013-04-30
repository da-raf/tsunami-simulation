#include <cxxtest/TestSuite.h>
#define private public

#include "FWave.hpp"

#define TYPE double
#define MAXERROR 1.0E-9


using namespace solver;

/**
 * This is the cxxtest test-suite for the FWave template class.
 */
class FWaveSolverTest : public CxxTest::TestSuite
{
public:
    /**
     * this function will test the private function roeEigenvals
     * of the template FWave by performing two value checks
     */
    void testEigenvalues() {
        FWave<TYPE> solver;
        
        TYPE hl, hr;
        TYPE hul, hur;
        
        TYPE lambda[2];
        TYPE test_lambda[2];
        
        hl  = 1.0;
        hul = 1.0;
        hr  = 4.0;
        hur = 4.0;
        
        lambda[0] = -3.95227220577;
        lambda[1] =  5.95227220577;
        
        solver.roeEigenvals(hl, hr, hul, hur, test_lambda);
        
        TS_ASSERT_DELTA(lambda[0], test_lambda[0], MAXERROR);
        TS_ASSERT_DELTA(lambda[1], test_lambda[1], MAXERROR);
        
        
        hl  = 1.3;
        hul = 0.6;
        hr  = 2.9;
        hur = 4.1;
        
        lambda[0] = -3.50692249113;
        lambda[1] =  5.57074240685;
        
        solver.roeEigenvals(hl, hr, hul, hur, test_lambda);
        
        TS_ASSERT_DELTA(lambda[0], test_lambda[0], MAXERROR);
        TS_ASSERT_DELTA(lambda[1], test_lambda[1], MAXERROR);
    }
    
    /**
     * testFlux will perform a fast value check on the private function "flux"
     * of the template "FWave"
     */
    void testFlux() {
        FWave<TYPE> solver;
        
        TYPE h  = 4.3;
        TYPE hu = 2.5;
        
        TYPE fl[2];
        solver.flux(h, hu, fl);
        
        TS_ASSERT_DELTA(fl[0],           2.5, MAXERROR);
        TS_ASSERT_DELTA(fl[1], 92.1469383721, MAXERROR);
    }
    
    /**
     * testEigencoefficients will do a valueCheck on the function "eigencoeffis" of the
     * template "FWave"
     */
    void testEigencoefficients() {
        FWave<TYPE> solver;
        
        TYPE hul = 1.5;
        TYPE hur = 3.0;
        TYPE hl  = 1.0;
        TYPE hr  = 0.73572032979;
        
        TYPE eigenvals[2];
        TYPE ec[2];
        
        eigenvals[0] = -0.9;
        eigenvals[1] =  1.4;
        // df.x = 1.5;
        // df.y = 7.73290921233;
        
        solver.eigencoeffis(hl, hr, hul, hur, 0.0, 0.0, eigenvals, ec);
        
        TS_ASSERT_DELTA(ec[0], -2.44909096188, MAXERROR);
        TS_ASSERT_DELTA(ec[1],  3.94909096188, MAXERROR);
    }
    
    /**
     * testSteadyState will calculate the net updates for identical water columns and
     * momentum on both sides, which have to be equal to zero
     */
    void testSteadyState() {
        FWave<TYPE> solver;
        
        TYPE hl, hr, hul, hur;
        
        TYPE updateLeft[2], updateRight[2];
        TYPE maxSpeed;
        
        hl  = hr  = 2.4;
        hul = hur = 3.7;
        
        solver.computeNetUpdates(
                hl, hr, hul, hur, 0.f, 0.f,
                updateLeft[0], updateRight[0],
                updateLeft[1], updateRight[1],
                maxSpeed );
        
        TS_ASSERT_DELTA(updateLeft[0], 0.0, MAXERROR);
        TS_ASSERT_DELTA(updateLeft[1], 0.0, MAXERROR);
        TS_ASSERT_DELTA(updateRight[0], 0.0, MAXERROR);
        TS_ASSERT_DELTA(updateRight[1], 0.0, MAXERROR);
        
    }
    
    /**
     * testSupersonic will check the function computeNetUpdates for correct behavior
     * in the supersonic case (both eigenvalues greater/less than zero)
     * - greater than zero: the net-updates on the left have to be zero
     * - less than zero: the net-updates on the right have to be zero
     */
    void testSupersonic() {
        FWave<TYPE> solver;
        
        TYPE updateLeft[2], updateRight[2];
        TYPE maxSpeed;
        
        // supersonic wave from left to right
        TYPE hl = 1.0;
        TYPE hr = 1.0;
        TYPE hul = 4.0;
        TYPE hur = 6.0;
        
        solver.computeNetUpdates(hl, hr, hul, hur, 0.f, 0.f,
                                 updateLeft[0], updateRight[0],
                                 updateLeft[1], updateRight[1],
                                 maxSpeed );
        TS_ASSERT_DELTA(updateLeft[0], 0.0, MAXERROR);
        TS_ASSERT_DELTA(updateLeft[1], 0.0, MAXERROR);
        
        // supersonic wave from right to left
        hul = -4.0;
        hur = -6.0;
        
        solver.computeNetUpdates(hl, hr, hul, hur, 0.f, 0.f,
                                 updateLeft[0], updateRight[0],
                                 updateLeft[1], updateRight[1],
                                 maxSpeed );
        TS_ASSERT_DELTA(updateRight[0], 0.0, MAXERROR);
        TS_ASSERT_DELTA(updateRight[1], 0.0, MAXERROR);
        
    }
    
    /**
     * check the handling of bathymetric data
     * 
     * a steady state will be created, but with different water depths
     */
    void testBathymetry() {
        FWave<TYPE> solver;
        
        TYPE hl = 10.0;
        TYPE hr =  5.0;
        TYPE hul = 0.0;
        TYPE hur = 0.0;
        
        TYPE bl = -5.0;
        TYPE br =  0.0;
        // notice: hl + bl = 5 = hr + br
        
        TYPE updateLeft[2], updateRight[2];
        TYPE maxSpeed;
        
        solver.computeNetUpdates(hl, hr, hul, hur, bl, br,
                                 updateLeft[0], updateRight[0],
                                 updateLeft[1], updateRight[1],
                                 maxSpeed );
        TS_ASSERT_DELTA(updateLeft[0], 0.0, MAXERROR);
        TS_ASSERT_DELTA(updateLeft[1], 0.0, MAXERROR);
        TS_ASSERT_DELTA(updateRight[0], 0.0, MAXERROR);
        TS_ASSERT_DELTA(updateRight[1], 0.0, MAXERROR);
    }
    
    void testDryCells() {
        FWave<TYPE> solver;
        
        TYPE h  = 5.0;
        TYPE hu = 5.0;
        
        TYPE bl = -2.0;
        TYPE br =  1.0;
        
        TYPE updateLeft[2], updateRight[2];
        TYPE maxSpeed;
        
        solver.computeNetUpdates(h, 0.0, hu, 0.0, bl, br,
                                 updateLeft[0], updateRight[0],
                                 updateLeft[1], updateRight[1],
                                 maxSpeed );
        TS_ASSERT_DELTA(updateRight[0], 0.0, MAXERROR);
        TS_ASSERT_DELTA(updateRight[1], 0.0, MAXERROR);
    }
    
};
