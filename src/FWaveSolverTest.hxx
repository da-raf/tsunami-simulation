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
        
        lambda1 = -3.95227220577;
        lambda2 = 5.95227220577;
        
        Vector2D test_lambda = roe_eigenvals(ql, qr);
        
        TS_ASSERT_DELTA(lambda1, test_lambda.x, 1.0E-9);
        TS_ASSERT_DELTA(lambda2, test_lambda.y, 1.0E-9);
        
        ql.h  = 1.3;
        ql.hu = 0.6;
        qr.h  = 2.9;
        qr.hu = 4.1;
        
        lambda1 = -3.50692249113;
        lambda2 = 5.57074240685;
        
        test_lambda = roe_eigenvals(ql, qr);
        
        TS_ASSERT_DELTA(lambda1, test_lambda.x, 1.0E-9);
        TS_ASSERT_DELTA(lambda2, test_lambda.y, 1.0E-9);
    }
    
    void testFlux() {
        State q;
        q.h = 4.3;
        q.hu = 2.5;
        
        Vector2D flux_test = flux(q);
        
        TS_ASSERT_DELTA(flux_test.x, 2.5, 1.0E-9);
        TS_ASSERT_DELTA(flux_test.y, 96.94345, 1.0E-9);
    }
    
    void testEigencoefficients() {
        Vector2D eigenvals, df;
        eigenvals.x = -0.9;
        eigenvals.y = 1.4;
        df.x = 1.5;
        df.y = 4.5;
        
        Vector2D ec = eigencoeffis(eigenvals, df);
        
        TS_ASSERT_DELTA(ec.x, -1.04347826087, 1.0E-9);
        TS_ASSERT_DELTA(ec.y, 2.54347826087, 1.0E-9);
    }
    
    void testSteadyState() {
        State ql, qr;
        
        ql.h = qr.h = 2.4;
        ql.hu = qr.hu = 3.7;
        
        Vector2D* z = calculate_updates(ql, qr);
        TS_ASSERT_DELTA(z[0].x, 0.0, 1.0E-9);
        TS_ASSERT_DELTA(z[0].y, 0.0, 1.0E-9);
        TS_ASSERT_DELTA(z[1].y, 0.0, 1.0E-9);
        TS_ASSERT_DELTA(z[1].y, 0.0, 1.0E-9);
        
    }
    
    // void testSupersonic();
};
