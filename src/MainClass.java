import java.util.Arrays;

import static java.lang.Math.sqrt;

public class MainClass {
    public static void main(String[] args) {
        MatrixSimple eig = new MatrixSimple(new double[][]{
                {5,2,-2},
                {2,5,-2},
                {-2,-2,5}
        });

        MatrixSimple eig2 = new MatrixSimple(new double[][]{
                {5,2,-2},
                {10,4,-4},
                {-10,-4,4}
        });

        MatrixSimple markov = new MatrixSimple(new double[][]{
                {0.3,0.1,0.6},
                {0.5,0.3,0.2},
                {0.8,0.1,0.1}
        });

        MatrixSimple eVector1 = new MatrixSimple(new double[][]{
                {0.5773502709817891},
                {0.5773502709817891},
                {0.5773502709817891}
        });

        MatrixSimple noEigenValues = new MatrixSimple(new double[][]{
                {0,-1},
                {1,0}
        });

        MatrixSimple noEig2 = new MatrixSimple(new double[][]{
                {-1.0/4, (1/4.0)+(1/sqrt(2)), -(1.0/2)+1/(2*sqrt(2))},
                {(1.0/4)-1/sqrt(2), -1.0/4, -(1.0/2)-1/(2*sqrt(2))},
                {-(1.0/2)-(1/(2*sqrt(2))), -(1.0/2)+(1/(2*sqrt(2))), 1.0/2}
        });
        System.out.println(Arrays.toString(noEig2.eigenvalues()));
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
                {1,0,-3,0,2,-8},
                {0,1,5,0,-1,4},
                {0,0,0,1,7,-9},
                {0,0,0,0,0,0}
        });

        MatrixSimple test2 = new MatrixSimple(new double[][]{
                {3,6,1},
                {1,2,1},
                {1,2,1},
                {1,2,1}
        });

        MatrixSimple test3 = new MatrixSimple(new double[][]{
                {1,-2,3},
                {1,-3,5},
                {2,3,1},
                {-3,4,2}
        });

        //System.out.println(test.transpose().QRFactorization().Q.round(3));
        //System.out.println();
        //System.out.println(test.transpose().QRFactorization().R.round(3));


/*
        int row1 = 4-1; //elimination row
        int row2 = 1-1;
        int col = 0;

        double a = test3.get(row1,col);
        double b = test3.get(row2, col);

        double r = Math.hypot(a, b);

        double c = b / r;
        double s = -a / r;
        MatrixSimple A = MatrixSimple.givensRotationMatrix(test3.getNumRows(), row2,row1, c,s);

        MatrixSimple test3_new = A.transpose().multiply(test3);


        //System.out.println(test3_new);
        //System.out.println();

        row1 = 4-1; //elimination row
        row2 = 3-1;
        col = 1;

        a = test3_new.get(row1,col);
        b = test3_new.get(row2, col);

        r = Math.hypot(a, b);

        c = b / r;
        s = -a / r;

        MatrixSimple B = MatrixSimple.givensRotationMatrix(test3.getNumRows(), row2,row1, c,s);
*/
        //System.out.println(A.multiply(B.transpose().multiply(test3_new)));
        MatrixSimple f = zero(test3, 3, 0);
        f = zero(f, 2, 0);
        f = zero(f, 3,1);
        f = zero(f, 2, 1);
        //System.out.println(zero(f, 3, 1));

        /*
        f = test3.copy();
        System.out.println(f);
        System.out.println("---------------------------------");
        for (int col = 0; col < f.getNumCols(); col++) {
            for (int row = f.getNumRows()-1; row >= f.getNumRows()-min(test3.getNumRows(), test3.getNumCols())+col; row--) {
                f = zero(f, row, col);
                System.out.println(f);
                System.out.println("---------------------------------");
            }
        }
        System.out.println(f);
         */
        //test2.QRFactorization();

        /*
        MatrixSimple.QRFactorization qr = test.transpose().QRFactorization();

        System.out.println(qr.Q.multiply(qr.R).floatingPointFix(1e3));
        System.out.println();

        for (int col = 0; col < qr.Q.getNumCols(); col++) {
            /*System.out.println(qr.Q.getCol(col).floatingPointFix(1e3));
            System.out.println();
            System.out.println(test.multiply(qr.Q.getCol(col)).floatingPointFix(1e3));
            System.out.println("--------------------------------------------------");*/
            /*for (int col2 = 0; col2 < qr.Q.getNumCols(); col2++) {
                if(col != col2) {
                    double dot = HALMathUtil.floatingPointFix(qr.Q.getCol(col).toVector().dot(qr.Q.getCol(col2).toVector()));
                    if(dot != 0) {
                        System.out.println("AGHH BAD BAD BAD");
                        System.out.println(dot);
                        System.out.println();
                        System.out.println(qr.Q.getCol(col));
                        System.out.println();
                        System.out.println(qr.Q.getCol(col2));
                        System.out.println("-----------------------------");
                    }
                }
            }
        }
        //System.out.println(test.nullSpace().getBasisMatrix().floatingPointFix(1e3));
        /*System.out.println(qr.Q.floatingPointFix(1e3));
        System.out.println();
        System.out.println(qr.R.floatingPointFix(1e3));
        System.out.println();*/

        //System.out.println();
        //System.out.println(test2.rref());


        //MatrixSimple.QRFactorization qr = test2.transpose().QRFactorization();
        //System.out.println(qr.R);
        //System.out.println(MatrixSimple.floatingPointFix(qr.Q.multiply(qr.R),1e3));
        //System.out.println(test2.rank());
        //System.out.println(test.nullSpace());
        //System.out.println(test.augmentVertical(MatrixSimple.identityMatrix(test.getNumCols())).rref());


        //should have 1,2,3 as eigenvalues.
       //markov.diagonalize();

    }

    public static MatrixSimple zero(MatrixSimple m, int row, int col) {
        ExceptionChecker.assertTrue(row >= 1, new ArithmeticException("you're dumb Cole :/"));

        double a = m.get(row,col);
        double b = m.get(row - 1, col);

        double r = Math.hypot(a, b);

        double c = b / r;
        double s = -a / r;
        MatrixSimple A = MatrixSimple.givensRotationMatrix(m.getNumRows(), row-1,row, c,s);

        return A.multiply(m).round(3);
    }
}
