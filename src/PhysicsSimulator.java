import java.util.function.Function;

public class PhysicsSimulator {
    private final PhysicsModel model;
    private final DataInterpreter interpreter;
    public PhysicsSimulator(PhysicsModel model, DataInterpreter interpreter) {
        this.model = model;
        this.interpreter = interpreter;
    }

    public void run(double duration) {
        model.reset();
        MatrixSimple state = model.initialState;

        for (double t = 0; t <= duration; t += model.dt) {
            interpreter.interpret(state);
            state = model.update();
        }

        interpreter.finish();
    }
}
