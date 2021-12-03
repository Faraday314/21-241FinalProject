import static java.lang.Math.sqrt;

public class ConstantAccelPiecewiseModel extends ConstantAccel2DModel {

    private final double turnRadius;
    public ConstantAccelPiecewiseModel(MatrixSimple initialState, double turnRadius, double dt) {
        super(initialState, dt);
        this.turnRadius = turnRadius;
        /*
        State should be:
            x pos,
            x vel,
            x accel,
            y pos,
            y vel,
            y accel
         */
    }

    @Override
    public MatrixSimple update() {
        //Start Circle if past X = 0
        if(state.get(0,0) >= 0) {
            Point2D center = Point2D.getOrigin();

            Point2D loc = new Point2D(state.get(0,0), state.get(3,0));
            Point2D newLoc = new Point2D(loc.getX(), sqrt(turnRadius*turnRadius-loc.getX()*loc.getX()));
            Vector2D toCenter = newLoc.vectorTo(center);

            Vector2D vel = new Vector2D(state.get(1,0), state.get(4,0));

            toCenter = toCenter.scaleTo(vel.magnitudeSquared()/turnRadius);
            Vector2D newVel = toCenter.rotate(90, AngleUnit.DEGREES).scaleTo(vel.magnitude());

            state.set(0,0, newLoc.getX());
            state.set(1,0, newVel.getX());
            state.set(2,0, toCenter.getX());
            state.set(3,0, newLoc.getY());
            state.set(4, 0, newVel.getY());
            state.set(5, 0, toCenter.getY());
        }
        return super.update();
    }
}
