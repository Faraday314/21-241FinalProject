public class Line<P extends Point<P,V>, V extends Vector<V>>{
    protected P startPoint, endPoint;
    protected V vector;
    public Line(P startPoint, P endPoint) {
        ExceptionChecker.assertEqual(startPoint.dim(), endPoint.dim(), new ArithmeticException("dimensionality of points do not match"));
        this.startPoint = startPoint;
        this.endPoint = endPoint;
        vector = startPoint.vectorTo(endPoint);
    }

    public double distanceTo(P point) {
        ExceptionChecker.assertEqual(point.dim(), vector.dim(), new ArithmeticException("dimensionality of point and line do not match"));
        return vector.subtract(startPoint.vectorTo(point).project(vector)).magnitude();
    }

    public double length() {
        return startPoint.distanceTo(endPoint);
    }

    public V toVector() {
        return vector;
    }

    public boolean isOnLine(P point) {
        return startPoint.distanceTo(point) + point.distanceTo(endPoint) == length();
    }
}
