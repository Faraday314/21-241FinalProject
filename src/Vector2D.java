import com.sun.istack.internal.NotNull;

import static java.lang.Math.*;

public class Vector2D extends EuclideanVector<Vector2D> {
    private static final Vector2D ZERO_VECTOR = new Vector2D(0,0);

    public Vector2D(MatrixSimple vectorMatrix) {
        super(vectorMatrix);
        ExceptionChecker.assertTrue(vectorMatrix.getNumRows() == 2, new ArithmeticException("Matrix vector is wrong size."));
    }

    public Vector2D(Point2D start, Point2D end) {
        super(start.vectorTo(end).vectorMatrix);
    }
    public Vector2D(Point2D end) {
        this(Point2D.getOrigin(), end);
    }

    public Vector2D(double x, double y) {
        super(x, y);
    }

    public Vector2D(double r, double theta, AngleUnit angleUnit) {
        super(CoordinateSystem2D.POLAR.convertTo(
                CoordinateSystem2D.CARTESIAN
        ).apply(
                new double[] {r, angleUnit.convertTo(AngleUnit.RADIANS).apply(theta)}
        ));
    }

    public static Vector2D fromND(VectorND vector) {
        ExceptionChecker.assertEqual(vector.dim(), 2, new ArithmeticException("Vector dimensionality does not match."));
        return new Vector2D(vector.vectorMatrix.copy());
    }

    @Override
    public Vector2D copy() {
        return new Vector2D(vectorMatrix.copy());
    }

    public double getX() {
        return vectorMatrix.get(0,0);
    }

    public double getY() {
        return vectorMatrix.get(1,0);
    }

    public void setX(double x) {
        vectorMatrix.set(0,0, x);
    }

    public void setY(double y) {
        vectorMatrix.set(0,0, y);
    }

    public double getTheta() {
        return atan2(getY(), getX());
    }
    public double getTheta(AngleUnit angleUnit) {
        return AngleUnit.RADIANS.convertTo(angleUnit).apply(getTheta());
    }

    public void setTheta(double thetaRadians) {
        double magnitude = magnitude();
        setX(magnitude*cos(thetaRadians));
        setY(magnitude*sin(thetaRadians));
    }
    public void setTheta(double theta, @NotNull AngleUnit angleUnit) {
        setTheta(angleUnit.convertTo(AngleUnit.RADIANS).apply(theta));
    }

    public Vector2D rotate(double angle, AngleUnit angleUnit) {
        Vector2D cpy = copy();
        if(!isZeroVector()) {
            double theta = angleUnit.convertTo(AngleUnit.RADIANS).apply(angle);
            double cosTheta = cos(theta);
            double sinTheta = sin(theta);

            double rotX = getX() * cosTheta - getY() * sinTheta;
            double rotY = getX() * sinTheta + getY() * cosTheta;

            cpy.setX(rotX);
            cpy.setY(rotY);
        }
        return cpy;
    }

    public Vector2D rotate(double angle) {
        return rotate(angle, AngleUnit.RADIANS);
    }

    public Vector3D cross(Vector2D vector) {
        return new Vector3D(
                0,
                0,
                getX()*vector.getY() + vector.getX()*getY());
    }

    public static Vector2D zeroVector() {
        return ZERO_VECTOR.copy();
    }
}
