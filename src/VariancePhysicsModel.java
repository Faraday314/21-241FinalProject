import org.jetbrains.annotations.Contract;
import org.jetbrains.annotations.NotNull;

import java.util.Random;

import static java.lang.Math.sqrt;

public class VariancePhysicsModel extends PhysicsModel {
    protected final PhysicsModel model;
    protected final MatrixSimple variance;
    private final Random r = new Random();
    public VariancePhysicsModel(@NotNull PhysicsModel model, @NotNull MatrixSimple variance) {
        super(model.initialState.copy(), model.dt);
        ExceptionChecker.assertTrue(model.initialState.getNumRows() == variance.getNumRows() && variance.isVector(), new ArithmeticException("Variance matrix is not the same size as the state vector. This is a problem."));
        this.model = model;
        this.variance = variance;
    }

    private double getRandomGaussian(double variance, double mean) {
        return r.nextGaussian()*sqrt(variance) + mean;
    }

    @Contract("_ -> param1")
    private @NotNull MatrixSimple addNoise(@NotNull MatrixSimple state) {
        int vals = state.getNumRows();
        for (int i = 0; i < vals; i++) {
            double var = variance.get(i, 0);
            double mean = state.get(i, 0);
            state.set(i, 0, getRandomGaussian(var, mean));
        }
        return state;
    }

    @Override
    public MatrixSimple getTransitionMatrix() {
        return model.getTransitionMatrix();
    }

    @Override
    public MatrixSimple update() {
        return addNoise(model.update());
    }

    @Override
    public MatrixSimple getState() {
        return addNoise(model.getState().copy());
    }
}
