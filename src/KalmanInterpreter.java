import org.jetbrains.annotations.NotNull;
import org.jfree.data.xy.XYSeries;

import java.awt.*;

public class KalmanInterpreter implements DataInterpreter {
    private final LinearKalmanFilter kf;
    private final VariancePhysicsModel model;
    private final XYSeries
            trueValue = new XYSeries("True Value", false),
            measuredValue = new XYSeries("Measured Value", false),
            kalmanEstimation = new XYSeries("Kalman Filter Estimation", false),
            squareError = new XYSeries("Error^2", false);
    private final LineChart chart = new LineChart("Particle Path", "X (meters)", "Y (meters)");
    private final LineChart errorChart = new LineChart("Error", "Time (seconds)", "Error^2 (m)");
    private final MatrixSimple observationMatrix;
    private double t = 0;
    public KalmanInterpreter(@NotNull LinearKalmanFilter kf, PhysicsModel model, MatrixSimple measurementVariance) {
        this.kf = kf;
        this.model = new VariancePhysicsModel(model, measurementVariance);
        this.observationMatrix = kf.getObservationMatrix();
    }

    @Override
    public void interpret(MatrixSimple state) {
        MatrixSimple truth = observationMatrix.multiply(state);
        trueValue.add(truth.get(0,0), truth.get(1,0));

        MatrixSimple measurement = observationMatrix.multiply(model.getState());
        measuredValue.add(measurement.get(0,0), measurement.get(1,0));

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
