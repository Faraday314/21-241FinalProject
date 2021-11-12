public class VectorND extends math.EuclideanVector<VectorND> {

    public VectorND(MatrixSimple matrixVector) {
        super(matrixVector);
    }

    public VectorND(double... components) {
        this(new MatrixSimple(components).transpose());
    }

    @Override
    public VectorND copy() {
        return new VectorND(vectorMatrix.copy());
    }

    @Override
    public VectorND set(VectorND vector) {
        return set(vector.vectorMatrix);
    }

    @Override
    public VectorND set(MatrixSimple vectorMatrix) {
        ExceptionChecker.assertTrue(vectorMatrix.isVector(), new ArithmeticException("matrix is not a vector."));
        ExceptionChecker.assertEqual(dim(), vectorMatrix.getNumRows(), new ArithmeticException("dimensionality of vectors does not match."));

        for (int row = 0; row < vectorMatrix.getNumRows(); row++) {
            this.vectorMatrix.set(row, 0, vectorMatrix.get(row, 0));
        }
        return this;
    }
}
