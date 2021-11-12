public interface Vector <V extends Vector<V>> {
    V copy();
    MatrixSimple toMatrix();

    double[] getComponents();

    boolean isZeroVector();
    boolean isUnitVector();
    boolean isNormalTo(V vector);

    double dot(V vector);
    double magnitude();
    double angleTo(V vector);

    int dim();

    V normalize();
    V scaleTo(double magnitude);
    V multiply(double scalar);
    V divide(double scalar);
    V add(V vector);
    V subtract(V vector);
    V project(V ontoVector);
    V rotate(math.Axis<V> u, math.Axis<V> v, double angleRadians);
    V set(V vector);
    V set(MatrixSimple vectorMatrix);
}
