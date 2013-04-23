#include <cxxtest/TestSuite.h>
#include "FWave.h"

#define PRIVATE_PUBLIC

using namespace std;

class FWaveSolverTest : public CxxTest::TestSuite
{
public:
    /* void xtestEigenvalues() {
        FWave<float> solver;
        
        float hl, hr;
        float hul, hur;
        
        float lambda[2];
        float test_lambda[2];
        
        hl  = 1.0;
        hul = 1.0;
        hr  = 4.0;
        hur = 4.0;
        
        lambda[0] = -3.95227220577;
        lambda[1] =  5.95227220577;
        
        solver.roeEigenvals(hl, hr, hul, hur, test_lambda[0], test_lambda[1]);
        
        TS_ASSERT_DELTA(lambda[0], test_lambda[0], 1.0E-9);
        TS_ASSERT_DELTA(lambda[1], test_lambda[1], 1.0E-9);
        
        
        hl  = 1.3;
        hul = 0.6;
        hr  = 2.9;
        hur = 4.1;
        
        lambda[0] = -3.50692249113;
        lambda[1] =  5.57074240685;
        
        solver.roeEigenvals(hl, hr, hul, hur, test_lambda[0], test_lambda[1]);
        
        TS_ASSERT_DELTA(lambda[0], test_lambda[0], 1.0E-9);
        TS_ASSERT_DELTA(lambda[1], test_lambda[1], 1.0E-9);
    }
    
    void xtestFlux() {
        FWave<float> solver;
        
        float h  = 4.3;
        float hu = 2.5;
        
        float fl[2];
        solver.flux(h, hu, fl[0], fl[1]);
        
        TS_ASSERT_DELTA(fl[0],      2.5, 1.0E-9);
        TS_ASSERT_DELTA(fl[1], 96.94345, 1.0E-9);
    }
    
    void xtestEigencoefficients() {
        FWave<float> solver;
        
        float hul = 1.5;
	    float hur = 3.0;
	    float hl  = 1.0;
	    float hr  = 0.73572032979;
        
	    float eigenvals[2];
        float ec[2];
        
        eigenvals[0] = -0.9;
        eigenvals[1] =  1.4;
        // df.x = 1.5;
        // df.y = 4.5;
        
        solver.eigencoeffis(hl, hr, hul, hur, eigenvals[0], eigenvals[1], ec[0], ec[1]);
        
        TS_ASSERT_DELTA(ec[0], -1.04347826087, 1.0E-9);
        TS_ASSERT_DELTA(ec[1],  2.54347826087, 1.0E-9);
    }
    */
    void testSteadyState() {
        FWave<float> solver;
        
        float hl, hr, hul, hur;
        
        float updateLeft[2], updateRight[2];
        float maxSpeed;
        
        hl  = hr  = 2.4;
        hul = hur = 3.7;
        
        solver.computeNetUpdates(
                hl, hr, hul, hur, 0.f, 0.f,
                updateLeft[0], updateRight[0],
                updateLeft[1], updateRight[1],
                maxSpeed );
        
        TS_ASSERT_DELTA(updateLeft[0], 0.0, 1.0E-9);
        TS_ASSERT_DELTA(updateLeft[1], 0.0, 1.0E-9);
        TS_ASSERT_DELTA(updateRight[0], 0.0, 1.0E-9);
        TS_ASSERT_DELTA(updateRight[1], 0.0, 1.0E-9);
        
    }
    
    // void testSupersonic();
};
