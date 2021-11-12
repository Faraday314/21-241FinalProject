public class VectorSpace {
    private final math.MatrixSimple basisMatrix;
    public VectorSpace(math.MatrixSimple[] vectors) {
        math.ExceptionChecker.assertTrue(vectors.length > 0, new ArithmeticException("A basis must have at least one vector."));

        math.MatrixSimple basisMatrix = vectors[0];
        for (int i = 1; i < vectors.length; i++) {
            basisMatrix = basisMatrix.augment(vectors[i]);
        }
        this.basisMatrix = basisMatrix;
        math.ExceptionChecker.assertTrue(basisMatrix.nullity() == 0, new ArithmeticException("Basis vectors are not linearly independent"));
    }

    public VectorSpace(math.MatrixSimple basisMatrix) {
        this.basisMatrix = basisMatrix;
        math.ExceptionChecker.assertTrue(basisMatrix.nullity() == 0, new ArithmeticException("Basis vectors are not linearly independent"));
    }

    public math.MatrixSimple getBasisMatrix() {
        return basisMatrix;
    }

    public math.MatrixSimple[] getBasisVectors() {
        return basisMatrix.getCols();
    }

    boolean contains(math.MatrixSimple vector) {
        math.ExceptionChecker.assertTrue(vector.isVector(), new ArithmeticException("matrix is not a vector."));

        math.MatrixSimple sol = basisMatrix.augment(vector).rref();
        for (int row = 0; row < sol.getNumRows(); row++) {
            for (int col = 0; col < sol.getNumCols(); col++) {
                if(sol.get(row, col) == 1 && col < basisMatrix.getNumCols()) {
                    break;
                }
                else if(sol.get(row, col) != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    @Override
    public String toString() {
        return basisMatrix.toString();
    }
}
