import org.jetbrains.annotations.NotNull;
import org.jfree.data.xy.XYSeries;

import java.awt.*;

public class PointInterpreter implements DataInterpreter {

    private final XYSeries xySeries = new XYSeries("True Value");
    private final LineChart chart = new LineChart("Particle Path", "X (meters)", "Y (meters)");
    @Override
    public void interpret(@NotNull MatrixSimple state) {
        ExceptionChecker.assertTrue(state.isVector(), new ArithmeticException("Cannot interpret a non-vector."));
        xySeries.add(state.get(0,0), state.get(3,0));
    }

    @Override
    public void finish() {
        chart.addData(xySeries, Color.GREEN);
        chart.render();
    }
}
