#include <cxxtest/TestSuite.h>
#include "fwave-solver.h"


class FWaveSolverTest : public CxxTest::TestSuite
{
public:
    void testEigenvalues() {
        State ql, qr;
        double lambda1, lambda2;
        
        ql.h  = 1.0;
        ql.hu = 1.0;
        qr.h  = 4.0;
        qr.hu = 4.0;
        
        lambda1 = 1 - 4.95227220577;
        lambda2 = 1 + 4.95227220577;
        
        Vector2D test_lambda = roe_eigenvals(ql, qr);
        
        TS_ASSERT_DELTA(lambda1, test_lambda.x, 1.0E-9);
        TS_ASSERT_DELTA(lambda2, test_lambda.y, 1.0E-9);
        
    }
    
    // void testSteadyState();
    // void testSupersonic();
};
