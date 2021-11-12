public class Point2D extends EuclideanPoint<Point2D, math.Vector2D> {
    private static final Point2D ORIGIN = new Point2D(0,0);

    public Point2D(double a, double b, CoordinateSystem2D coordinateSystem) {
        super(new math.Vector2D(
                new MatrixSimple(coordinateSystem.convertTo(CoordinateSystem2D.CARTESIAN).apply(
                        new double[] {a, b}
                )).transpose()
        ));
    }

    public Point2D(double x, double y) {
        this(x, y, CoordinateSystem2D.CARTESIAN);
    }

    public Point2D(double r, double theta, AngleUnit angleUnit) {
        this(r, angleUnit.convertTo(AngleUnit.RADIANS).apply(theta), CoordinateSystem2D.POLAR);
    }
    public Point2D(math.Vector2D vector) {
        super(vector);
    }

    @Override
    public Point2D copy() {
        return new Point2D(pointVector);
    }

    public double getX() {
        return pointVector.getX();
    }

    public double getY() {
        return pointVector.getY();
    }

    public double getR() {
        return pointVector.magnitude();
    }
    public double getTheta() {
        return pointVector.getTheta();
    }
    public double getTheta(AngleUnit angleUnit) {
        return pointVector.getTheta(angleUnit);
    }
    public void setX(double x) {
        pointVector.setX(x);
    }
    public void setY(double y) {
        pointVector.setY(y);
    }
    public void setR(double r) {
        pointVector.scaleTo(r);
    }
    public void setTheta(double thetaRadians) {
        pointVector.setTheta(thetaRadians);
    }
    public void setTheta(double theta, AngleUnit angleUnit) {
        pointVector.setTheta(theta, angleUnit);
    }

    public static Point2D getOrigin() {
        return ORIGIN.copy();
    }

    @Override
    public String toString() {
        return "(" + getX() + ", " + getY() + ")";
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Point2D point2D = (Point2D) o;
        return Double.compare(point2D.getX(), getX()) == 0 &&
                Double.compare(point2D.getY(), getY()) == 0;
    }
}
