#define G 9.81


typedef struct {
    double x;
    double y;
} Vector2D;


typedef struct {
    double h;
    double hu;
} State;


/** calculate the result of the flux function for this quantity
 * @param q quantity
 */
Vector2D flux(State q);

/** calculate the Roe eigenvalues for the left and right quantities
 * @param ql left quantity
 * @param qr right quantity
 * @result Roe eigenvalues (x: \lambda_l, y: lambda_r)
 */
Vector2D roe_eigenvals(State ql, State qr);

/** calculate the eigencoefficients
 */
Vector2D eigencoeffis(State ql, State qr, Vector2D roe_eigenvals);

Vector2D* calculate_updates(State ql, State qr);
