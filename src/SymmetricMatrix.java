import com.sun.istack.internal.NotNull;

public class SymmetricMatrix extends MatrixSimple {
    public SymmetricMatrix(MatrixSimple matrix) {
        super(matrix.vals);
        ExceptionChecker.assertTrue(isSymmetric(), new ArithmeticException("Matrix is not symmetric."));
    }

    @Override
    public EigenData eigen() {
        ExceptionChecker.assertTrue(this.isSquare(), new RuntimeException("Matrix must be square!"));

        QRFactorization qrFactorization = this.QRFactorization();
        MatrixSimple A = qrFactorization.R.multiply(qrFactorization.Q);
        MatrixSimple S = qrFactorization.Q;

        MatrixSimple identity = MatrixSimple.identityMatrix(qrFactorization.Q.getNumRows());

        while(!A.round(7).isDiagonal()) {

            double s = A.get(getNumRows()-1, getNumRows()-1);
            MatrixSimple smult = MatrixSimple.identityMatrix(A.getNumRows()).multiply(s);

            qrFactorization = A.subtract(smult).QRFactorization();

            A = qrFactorization.R.multiply(qrFactorization.Q).add(smult);
            S = S.multiply(qrFactorization.Q);

            //If matrix is already upper triangular
            if(qrFactorization.Q.equals(identity)) {
                for (int diag = 0; diag < A.getNumRows(); diag++) {
                    identity.set(diag, diag, A.get(diag, diag));
                }
                A = identity;
                break;
            }
        }

        return new SpectralDecomposition(A, S);
    }

    public static class SpectralDecomposition extends EigenData {
        protected SpectralDecomposition(MatrixSimple diagEigenvalueMatrix, @NotNull MatrixSimple eigenVectorMatrix) {
            super(diagEigenvalueMatrix, eigenVectorMatrix);
            //TODO ensure eigenvectormatrix is orthogonal (orthonormal cols and unit eigenvectors)
        }
    }
}
