import java.math.BigInteger;
import java.util.Arrays;
import java.util.stream.IntStream;

import static java.lang.Math.*;

public class MainClass {
    public static void main(String[] args) {
        MatrixSimple eig = new MatrixSimple(new double[][]{
                {5, 2, -2},
                {2, 5, -2},
                {-2, -2, 5}
        });

        MatrixSimple eig2 = new MatrixSimple(new double[][]{
                {5, 2, -2},
                {10, 4, -4},
                {-10, -4, 4}
        });

        MatrixSimple markov = new MatrixSimple(new double[][]{
                {0.3, 0.1, 0.6},
                {0.5, 0.3, 0.2},
                {0.8, 0.1, 0.1}
        });

        MatrixSimple eVector1 = new MatrixSimple(new double[][]{
                {0.5773502709817891},
                {0.5773502709817891},
                {0.5773502709817891}
        });

        MatrixSimple noEigenValues = new MatrixSimple(new double[][]{
                {0, -1},
                {1, 0}
        });

        MatrixSimple noEig2 = new MatrixSimple(new double[][]{
                {-1.0 / 4, (1 / 4.0) + (1 / sqrt(2)), -(1.0 / 2) + 1 / (2 * sqrt(2))},
                {(1.0 / 4) - 1 / sqrt(2), -1.0 / 4, -(1.0 / 2) - 1 / (2 * sqrt(2))},
                {-(1.0 / 2) - (1 / (2 * sqrt(2))), -(1.0 / 2) + (1 / (2 * sqrt(2))), 1.0 / 2}
        });

        //System.out.println(Arrays.toString(noEig2.eigenvalues()));
/*
        System.out.println(markov.multiply(eVector1));
        System.out.println();
        System.out.println(eVector1);
        System.out.println();
        System.out.println(eVector1.augment(markov.multiply(eVector1)).nullity());
        System.out.println();
        System.out.println(markov.rref());
*/
        MatrixSimple test = new MatrixSimple(new double[][]{
                {1, 0, -3, 0, 2, -8},
                {0, 1, 5, 0, -1, 4},
                {0, 0, 0, 1, 7, -9},
                {0, 0, 0, 0, 0, 0}
        });

        MatrixSimple test2 = new MatrixSimple(new double[][]{
                {3, 6, 1},
                {1, 2, 1},
                {1, 2, 1},
                {1, 2, 1}
        });

        MatrixSimple test3 = new MatrixSimple(new double[][]{
                {1, -2, 3},
                {1, -3, 5},
                {2, 3, 1},
                {-3, 4, 2}
        });

        //System.out.println(eig2.determinant());
        ComplexNumber abc = new ComplexNumber(10,10);
        ComplexNumber abc2 = new ComplexNumber(10, 0.9352797163456124);
        ComplexNumber def = new ComplexNumber(1,3);
        //System.out.println(ComplexNumber.fromReal(2).pow(abc));
        //System.out.println(ComplexNumber.fromReal(2).pow(abc2));

        //should have 1,2,3 as eigenvalues.
        //markov.diagonalize();

        //ConstantModel model = new ConstantModel((new Vector2D(1,2)).toMatrix(), 0.25);
        VectorND initialState = new VectorND(-400,24.89,0,300,0,0);
        ConstantAccelPiecewiseModel model = new ConstantAccelPiecewiseModel(initialState.toMatrix(), 300,1);
        //ConstantAccel2DModel model = new ConstantAccel2DModel(initialState.toMatrix(), 1);


        //VectorND noiseVector = new VectorND(sigma_x*sigma_x,9,9,sigma_y*sigma_y,9,9);
        //PointInterpreter interpreter = new PointInterpreter();

        TestKalmanInterpreter interpreter = new TestKalmanInterpreter();

        PhysicsSimulator sim = new PhysicsSimulator(model, interpreter);
        sim.run(34);



        //LinearKalmanFilter kf = new LinearKalmanFilter(model, );
/*
        PhysicsSimulator sim = new PhysicsSimulator(model, interpreter, (MatrixSimple state) ->
                new MatrixSimple(state.getRowArray(0), state.getRowArray(3)));
        sim.run(3);*/
    }
}
