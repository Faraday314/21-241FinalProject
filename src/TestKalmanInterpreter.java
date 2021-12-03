import org.jetbrains.annotations.NotNull;
import org.jfree.data.xy.XYSeries;

import java.awt.*;

/**
 * A test interpreter for Kalman filter data. Note: initial state must be (-400,24.89,0,300,0,0) and turning radius must be 300 meters for this to work.
 */
public class TestKalmanInterpreter implements DataInterpreter {
    private final LinearKalmanFilter kf;
    private final XYSeries
            trueValue = new XYSeries("True Value", false),
            measuredValue = new XYSeries("Measured Value", false),
            kalmanEstimation = new XYSeries("Kalman Filter Estimation", false),
            squareError = new XYSeries("Error^2", false);

    //Sample measurement dataset, used for repeatability.
    private final double[] xMeasures = new double[] {-393.66, -375.93, -351.04, -328.96,-299.35,-273.36,	-245.89,	-222.58,	-198.03,	-174.17,	-146.32,	-123.72,	-103.47,	-78.23,	-52.63,	-23.34,	25.96,	49.72,	76.94,	95.38,	119.83,	144.01,	161.84,	180.56,	201.42,	222.62,	239.4,	252.51,	266.26,	271.75,	277.4,	294.12,	301.23,	291.8,	299.89};
    private final double[] yMeasures = new double[] {300.4,	301.78,	295.1,	305.19,	301.06,	302.05,	300,	303.57,	296.33,	297.65,	297.41,	299.61,	299.6,	302.39,	295.04,	300.09,	294.72,	298.61,	294.64,	284.88,	272.82,	264.93,	251.46,	241.27,	222.98,	203.73,	184.1,	166.12,	138.71,	119.71,	100.41,	79.76,	50.62,	32.99,	2.14};

    private int i = 0;
    private final LineChart chart = new LineChart("Particle Path", "X (meters)", "Y (meters)");
    private final LineChart errorChart = new LineChart("Error", "Time (seconds)", "Error^2 (m)");
    private final MatrixSimple observationMatrix;
    private double t = 0;

    public TestKalmanInterpreter() {
        double sigma_x = 3;
        double sigma_y = 3;
        double sigma_accel = 0.15;
        double init_uncertainty = 500;

        VectorND initialState = new VectorND(-400,24.89,0,300,0,0);
        ConstantAccelPiecewiseModel model = new ConstantAccelPiecewiseModel(initialState.toMatrix(), 300, 1);


        MatrixSimple H = new MatrixSimple(new double[][] {
                {1,0,0,0,0,0},
                {0,0,0,1,0,0}
        });
        MatrixSimple R = new MatrixSimple(new double[][] {
                {sigma_x*sigma_x, 0},
                {0, sigma_y*sigma_y}
        });
        MatrixSimple Q = model.calculateProcessNoise(sigma_accel);
        MatrixSimple cov_init = MatrixSimple.identityMatrix(6).multiply(init_uncertainty);
        MatrixSimple init_estimate = MatrixSimple.zeroMatrix(6,1);

        this.kf = new LinearKalmanFilter(model, init_estimate, cov_init, R, H, Q);;
        this.observationMatrix = kf.getObservationMatrix();
    }

    @Override
    public void interpret(MatrixSimple state) {
        MatrixSimple truth = observationMatrix.multiply(state);
        trueValue.add(truth.get(0,0), truth.get(1,0));

        MatrixSimple measurement = (new Vector2D(xMeasures[i], yMeasures[i])).toMatrix();
        measuredValue.add(measurement.get(0,0), measurement.get(1,0));
        i++;

        kf.update(measurement);
        MatrixSimple estimate = observationMatrix.multiply(kf.getEstimate());
        kalmanEstimation.add(estimate.get(0,0), estimate.get(1,0));

        MatrixSimple error = estimate.subtract(truth);
        squareError.add(t, error.toVector().magnitudeSquared());

        t += kf.getModel().dt;
    }

    @Override
    public void finish() {
        chart.addData(trueValue, Color.GREEN);
        chart.addData(measuredValue, Color.BLUE);
        chart.addData(kalmanEstimation, Color.RED);
        chart.render();

        errorChart.addData(squareError, Color.MAGENTA);
        errorChart.render();
    }
}
