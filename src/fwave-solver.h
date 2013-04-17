#define G 9.81


typedef struct {
    double x;
    double y;
} Vector2D;


typedef struct {
    double h;
    double hu;
} State;


Vector2D flux(State q);
Vector2D roe_eigenvals(State ql, State qr);
Vector2D eigencoeffis(Vector2D roe_eigenvals, Vector2D df);
Vector2D* calculate_updates(State ql, State qr);
