import config_testing.HALMathUtil;
import org.jetbrains.annotations.NotNull;

import java.util.*;

import static java.lang.Math.max;
import static java.lang.Math.pow;

/**
 * A matrix class for doing matrix math.
 *
 * Creation Date: 5/27/20
 *
 * @author Cole Savage, Level Up
 * @since 1.1.0
 * @version 1.0.0
 *
 * @see MatrixSimple
 */
public class MatrixSimple {

    //The double array matrix of values contained in the matrix.
    protected final double[][] vals;

    /**
     * A constructor for matrix where individual arrays of doubles can be entered or a matrix of doubles can be entered.
     *
     * @param vals The list of rows in the matrix.
     * @throws RuntimeException Throws this exception if the entered rows are not all of the same length.
     */
    public MatrixSimple(double[]... vals) {
        if(vals == null) {
            this.vals = new double[1][1];
            return;
        }

        for(double[] row : vals) math.ExceptionChecker.assertTrue(row.length == vals[0].length, new RuntimeException("Entered rows in matrix must all be the same length."));
        this.vals = vals.clone();
    }

    public MatrixSimple(int rows, int cols) {
        this.vals = new double[rows][cols];
    }

    /**
     * A private constructor for matrix used for cloning.
     *
     * @param matrix The matrix to clone.
     */
    private MatrixSimple(MatrixSimple matrix) {
        vals = new double[matrix.getNumRows()][matrix.getNumCols()];
        for (int i = 0; i < matrix.getNumRows(); i++) System.arraycopy(matrix.vals[i], 0, vals[i], 0, vals[i].length);
    }

    public MatrixSimple copy() {
        return new MatrixSimple(this);
    }

    /**
     * Sets an entry in the matrix at the given row and column to the desired value.
     *
     * @param row The row index of the entry to set.
     * @param col The column index of the entry to set.
     * @param value The value to set the entry in the matrix at (row, column) to.
     * @throws RuntimeException Throws this exception if the row index isn't within the range [0, getNumRows()) or the column index isn't within the range [0, getNumCols()).
     */
    public void set(int row, int col, double value) {
        math.ExceptionChecker.assertTrue(row < getNumRows(), new RuntimeException("Cannot set an entry in a row outside the matrix."));
        math.ExceptionChecker.assertTrue(row >= 0, new RuntimeException("Cannot set an entry in a negative row."));
        math.ExceptionChecker.assertTrue(col < getNumCols(), new RuntimeException("Cannot set an entry in a column outside the matrix."));
        math.ExceptionChecker.assertTrue(col >= 0, new RuntimeException("Cannot set an entry in a negative column."));

        vals[row][col] = value;
    }

    /**
     * Gets the value in the matrix at the specified row, column entry.
     *
     * @param row The row index.
     * @param col The column index.
     * @return The entry at the given row and column in the matrix.
     * @throws RuntimeException Throws this exception if the row index isn't within the range [0, getNumRows()) or the column index isn't within the range [0, getNumCols()).
     */
    public double get(int row, int col) {
        math.ExceptionChecker.assertTrue(row < getNumRows(), new RuntimeException("Cannot get an entry in a row outside the matrix."));
        math.ExceptionChecker.assertTrue(row >= 0, new RuntimeException("Cannot get an entry in a negative row."));
        math.ExceptionChecker.assertTrue(col < getNumCols(), new RuntimeException("Cannot get an entry in a column outside the matrix."));
        math.ExceptionChecker.assertTrue(col >= 0, new RuntimeException("Cannot get an entry in a negative column."));

        return vals[row][col];
    }

    /**
     * Gets the row at the specified row index in the matrix.
     *
     * @param row The row index.
     * @return The row at the given row index.
     * @throws RuntimeException Throws this exception if the row index isn't within the range [0, getNumRows()).
     */
    public double[] getRowArray(int row) {
        math.ExceptionChecker.assertTrue(row < getNumRows(), new RuntimeException("Cannot get an entry in a row outside the matrix."));
        math.ExceptionChecker.assertTrue(row >= 0, new RuntimeException("Cannot get an entry in a negative row."));

        return vals[row].clone();
    }

    /**
     * Gets the column at the specified column index in the matrix.
     *
     * @param col The column index.
     * @return The column at the given column index.
     * @throws RuntimeException Throws this exception if the row index isn't within the range [0, getNumCols()).
     */
    public double[] getColArray(int col) {
        math.ExceptionChecker.assertTrue(col < getNumCols(), new RuntimeException("Cannot get an entry in a column outside the matrix."));
        math.ExceptionChecker.assertTrue(col >= 0, new RuntimeException("Cannot get an entry in a negative column."));

        double[] colArray = new double[getNumRows()];
        for (int row = 0; row < getNumRows(); row++) colArray[row] = get(row, col);
        return colArray;
    }

    public MatrixSimple getCol(int col) {
        math.ExceptionChecker.assertTrue(col < getNumCols(), new RuntimeException("Cannot get an entry in a column outside the matrix."));
        math.ExceptionChecker.assertTrue(col >= 0, new RuntimeException("Cannot get an entry in a negative column."));

        MatrixSimple column = MatrixSimple.zeroMatrix(getNumRows(), 1);
        for (int row = 0; row < getNumRows(); row++) column.set(row, 0, get(row, col));
        return column;
    }

    public MatrixSimple getRow(int row) {
        return new MatrixSimple(getRowArray(row));
    }

    public MatrixSimple[] getCols() {
        MatrixSimple[] columns = new MatrixSimple[getNumCols()];
        for (int col = 0; col < getNumCols(); col++) {
            columns[col] = getCol(col);
        }
        return columns;
    }

    public MatrixSimple[] getRows() {
        MatrixSimple[] rows = new MatrixSimple[getNumRows()];
        for (int row = 0; row < getNumRows(); row++) {
            rows[row] = getRow(row);
        }
        return rows;
    }

    /**
     * Gets the total number of rows in the matrix.
     *
     * @return The total number of rows in the matrix.
     */
    public int getNumRows() {
        return vals.length;
    }

    /**
     * Gets the total number of columns in the matrix.
     *
     * @return The total number of columns in the matrix.
     */
    public int getNumCols() {
        if (this.getNumRows() == 0) return 0;
        return vals[0].length;
    }

    /**
     * Sets the row at the given row index to the given row.
     *
     * @param rowNum The row index.
     * @param row The row to add at that row index.
     * @throws RuntimeException Throws this exception if the length of the new row does not match the length of the old row.
     */
    public void setRow(int rowNum, double[] row) {
        math.ExceptionChecker.assertTrue(row.length == getNumCols(), new RuntimeException("New row length does not match matrix row length."));
        vals[rowNum] = row;
    }

    /**
     * Sets the column at the given column index to the given column.
     *
     * @param colNum The column index.
     * @param col The column to add at that column index.
     * @throws RuntimeException Throws this exception if the length of the new column does not match the length of the old column.
     */
    public void setCol(int colNum, double[] col) {
        math.ExceptionChecker.assertTrue(col.length == getNumRows(), new RuntimeException("New column length does not match matrix column length."));
        for (int i = 0; i < getNumRows(); i++) set(i, colNum, col[i]);
    }

    /**
     * Gets whether the matrix is the zero matrix.
     *
     * @return Whether the matrix is the zero matrix.
     */
    public boolean isZeroMatrix() {
        int rows = getNumRows();
        int cols = getNumCols();

        boolean isZeroMatrix = true;
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                isZeroMatrix &= HALMathUtil.round(get(row, col), 7) == 0;
            }
        }
        return isZeroMatrix;
    }

    /**
     * Gets whether the matrix is the identity matrix.
     *
     * @return Whether the matrix is the identity matrix.
     */
    public boolean isIdentityMatrix() {
        if (!isSquare()) return false;

        int rows = getNumRows();
        int cols = getNumCols();

        boolean isIdentity = true;
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                double val = get(row, col);

                if (row == col) isIdentity &= HALMathUtil.round(val, 7) == 1;
                else isIdentity &= HALMathUtil.round(val, 7) == 0;
            }
        }
        return isIdentity;
    }

    /**
     * Gets whether the matrix is a square matrix.
     *
     * @return Whether the matrix is a square matrix.
     */
    public boolean isSquare() {
        return this.getNumRows() == this.getNumCols();
    }

    /**
     * Gets whether the matrix is a diagonal matrix.
     *
     * @return Whether the matrix is a diagonal matrix.
     */
    public boolean isDiagonal() {
        math.ExceptionChecker.assertTrue(isSquare(), new ArithmeticException("Matrix is not square."));
        for (int row = 0; row < getNumRows(); row++) {
            for (int col = 0; col < getNumCols(); col++) {
                if(row != col && get(row, col) != 0) {
                    return false;
                }
            }
        }
        return  true;
    }

    /**
     * Gets whether the matrix is symmetric.
     *
     * @return Whether the matrix is symmetric.
     */
    public boolean isSymmetric() {
        return copy().transpose().equals(this);
    }

    /**
     * Calculates the transpose of the matrix.
     *
     * @return This matrix.
     */
    public MatrixSimple transpose() {
        int rows = getNumRows();
        int cols = getNumCols();

        double[][] matrixTranspose = new double[cols][rows];
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                matrixTranspose[col][row] = vals[row][col];
            }
        }
        return new MatrixSimple(matrixTranspose);
    }

    /**
     * Calculates the trace of this matrix.
     *
     * @return The trace of this matrix.
     */
    public double trace() {
        math.ExceptionChecker.assertTrue(isSquare(), new RuntimeException("Matrix is not square"));

        double trace = 0;
        for (int diag = 0; diag < getNumRows(); diag++) {
            trace += vals[diag][diag];
        }
        return trace;
    }

    public MatrixSimple augment(MatrixSimple matrix) {
        double[][] newVals = new double[max(getNumRows(),matrix.getNumRows())][getNumCols()+matrix.getNumCols()];
        for (int row = 0; row < getNumRows(); row++) {
            if (getNumCols() >= 0) System.arraycopy(vals[row], 0, newVals[row], 0, getNumCols());
        }

        for (int row = 0; row < matrix.getNumRows(); row++) {
            for (int col = 0; col < matrix.getNumCols(); col++) {
                newVals[row][getNumCols()+col] = matrix.vals[row][col];
            }
        }

        return new MatrixSimple(newVals);
    }

    public MatrixSimple augmentVertical(MatrixSimple matrix) {
        double[][] newVals = new double[getNumRows()+matrix.getNumRows()][max(getNumCols(),matrix.getNumCols())];

        for (int row = 0; row < getNumRows(); row++) {
            if (getNumCols() >= 0) System.arraycopy(vals[row], 0, newVals[row], 0, getNumCols());
        }

        for (int row = 0; row < matrix.getNumRows(); row++) {
            for (int col = 0; col < matrix.getNumCols(); col++) {
                newVals[getNumRows()+row][col] = matrix.vals[row][col];
            }
        }

        return new MatrixSimple(newVals);
    }

    /**
     * Calculates the determinant of this matrix.
     *
     * @return The determinant of this matrix.
     * @throws RuntimeException Throws this exception when the matrix is not square.
     */
    public double determinant() {
        math.ExceptionChecker.assertTrue(isSquare(), new RuntimeException("Is not a square matrix"));
        return determinant(this);
    }

    private double determinant(MatrixSimple subMatrix) {

        //if the sub-matrix is 2x2
        if(subMatrix.getNumRows() == 2) {
            return subMatrix.get(0,0)* subMatrix.get(1,1) - subMatrix.get(0,1)*subMatrix.get(1,0);
        }

        //Reminder: The matrix is square so getNumRows == getNumCols
        double[] zeroEntry = new double[subMatrix.getNumRows()];
        double determinantSum = 0;
        for (int col = 0; col < subMatrix.getNumCols(); col++) {
            double selectedVal = subMatrix.vals[0][col];

            MatrixSimple mask = MatrixSimple.onesMatrix(subMatrix.getNumRows(),subMatrix.getNumCols());
            mask.setRow(0, zeroEntry);
            mask.setCol(col, zeroEntry);

            MatrixSimple maskedSubMatrix = subMatrix.mask(mask);
            double[][] subMatrixVals = new double[subMatrix.getNumRows()-1][subMatrix.getNumCols()-1];
            int totalRowsFilled = 0, totalColsFilled;
            for (int subRow = 0; subRow < maskedSubMatrix.getNumRows(); subRow++) {
                if(subRow == 0) continue;

                totalColsFilled = 0;
                for (int subCol = 0; subCol < maskedSubMatrix.getNumCols(); subCol++) {
                    if(subCol != col) {
                        subMatrixVals[totalRowsFilled][totalColsFilled] = maskedSubMatrix.vals[subRow][subCol];
                        totalColsFilled++;
                    }
                }
                totalRowsFilled++;
            }

            MatrixSimple newSubMatrix = new MatrixSimple(subMatrixVals);
            double det = determinant(newSubMatrix);
            if(col % 2 == 0) {
                determinantSum += selectedVal * det;
            }
            else {
                determinantSum -= selectedVal * det;
            }
        }

        return determinantSum;
    }


    public MatrixSimple invert() {
        math.ExceptionChecker.assertTrue(isSquare(), new RuntimeException("Is not a square matrix"));
        MatrixSimple inverseFinderRREF = this.augment(MatrixSimple.identityMatrix(getNumRows())).rref();
        MatrixSimple potentialIdentity = inverseFinderRREF.crop(0, getNumRows(), 0, getNumCols());
        potentialIdentity = potentialIdentity.round(7);
        math.ExceptionChecker.assertTrue(potentialIdentity.isIdentityMatrix(), new RuntimeException("Matrix is not invertible!"));
        return inverseFinderRREF.crop(0, getNumRows(), getNumCols(), 2*getNumCols());
    }

    /**
     * Multiplies this matrix by a constant value.
     *
     * @param value The value to multiply by.
     * @return This matrix.
     */
    public MatrixSimple multiply(double value) {
        double[][] multipliedVals = new double[getNumRows()][getNumCols()];
        for (int row = 0; row < getNumRows(); row++) {
            for (int col = 0; col < getNumCols(); col++) {
                multipliedVals[row][col] = vals[row][col]*value;
            }
        }
        return new MatrixSimple(multipliedVals);
    }

    /**
     * Multiplies this matrix by another matrix.
     *
     * @param matrix The matrix to multiply by.
     * @return This matrix.
     * @throws RuntimeException Throws this exception when this matrix cannot be multiplied by the given matrix.
     */
    public MatrixSimple multiply(MatrixSimple matrix) {
        math.ExceptionChecker.assertTrue(this.getNumCols() == matrix.getNumRows(), new RuntimeException("Cannot multiply given matrices, invalid dimensions."));

        double[][] multipliedMatrix = new double[getNumRows()][matrix.getNumCols()];
        for (int rowNum = 0; rowNum < getNumRows(); rowNum++) {
            double[] row = getRowArray(rowNum);
            for (int colNum = 0; colNum < matrix.getNumCols(); colNum++) {
                double[] col = matrix.getColArray(colNum);
                double dotProduct = 0;
                for (int k = 0; k < row.length; k++) {
                    dotProduct += row[k]*col[k];
                }
                multipliedMatrix[rowNum][colNum] = dotProduct;
            }
        }
        return new MatrixSimple(multipliedMatrix);
    }

    /**
     * Divides this matrix by a constant value.
     *
     * @param value The value to divide by.
     * @return This matrix.
     * @throws ArithmeticException Throws this exception if you try to divide by 0.
     */
    public MatrixSimple divide(double value) {
        math.ExceptionChecker.assertFalse(value == 0, new ArithmeticException("Divide by zero error."));
        return multiply(1.0 / value);
    }

    public MatrixSimple add(MatrixSimple matrix) {
        math.ExceptionChecker.assertTrue(matrix.getNumCols() == this.getNumCols() && matrix.getNumRows() == this.getNumRows(), new RuntimeException("Matrices must be the same size!"));
        double[][] sumVals = new double[getNumRows()][getNumCols()];
        for (int row = 0; row < getNumRows(); row++) {
            for (int col = 0; col < getNumCols(); col++) {
                sumVals[row][col] = this.vals[row][col] + matrix.vals[row][col];
            }
        }

        return new MatrixSimple(sumVals);
    }

    public MatrixSimple subtract(MatrixSimple matrix) {
        return this.add(matrix.multiply(-1));
    }

    /**
     * Masks the matrix, so that only entries in this matrix with corresponding entries > 0 in the mask will remain (other entries are set to 0).
     *
     * @param mask The matrix being used to mask this matrix.
     * @return This matrix.
     * @throws RuntimeException Throws this exception when the size of the mask is not the same size as this matrix.
     */
    public MatrixSimple mask(MatrixSimple mask) {
        double[][] maskedVals = new double[getNumRows()][getNumCols()];
        for (int row = 0; row < getNumRows(); row++) {
            for (int col = 0; col < getNumCols(); col++) {
                if (mask.vals[row][col] > 0) maskedVals[row][col] = vals[row][col];
            }
        }

        return new MatrixSimple(maskedVals);
    }

    public MatrixSimple crop(int fromRow, int toRow, int fromCol, int toCol) {
        double[][] croppedVals = new double[toRow-fromRow][toCol-fromCol];
        for (int row = 0; row < getNumRows(); row++) {
            for (int col = 0; col < getNumCols(); col++) {
                if((row >= fromRow && row < toRow) && (col >= fromCol && col < toCol)) croppedVals[row-fromRow][col-fromCol] = vals[row][col];
            }
        }
        return new MatrixSimple(croppedVals);
    }

    /**
     * Calculates the row echelon form of this matrix using Gaussian elimination.
     *
     * @return This matrix.
     */
    public MatrixSimple rref() {
        LUFactorization luFactorization = this.LUFactorization();
        MatrixSimple rref = luFactorization.U.round(7);

        for (int pivot = Math.min(getNumRows(), getNumCols())-1; pivot >= 0; pivot--) {
            double pivotVal = rref.vals[pivot][pivot];

            if(HALMathUtil.round(pivotVal, 7) != 0) {
                //Divide pivot row by pivot value to make the pivot 1.
                double[] pivotRow = rref.getRowArray(pivot);
                math.FakeNumpy.divide(pivotRow, pivotVal);
                rref.setRow(pivot, pivotRow);

                for (int row = 0; row < pivot; row++) {
                    double elimVal = rref.vals[row][pivot];

                    //Generate elimination matrix
                    MatrixSimple E = MatrixSimple.identityMatrix(getNumRows());
                    E.set(row, pivot, -elimVal);

                    rref = E.multiply(rref);
                }
            }
        }

        boolean reachedZeroRow = false;
        for (int row = 0; row < rref.getNumRows(); row++) {
            boolean isZeroRow = true;
            for (int col = 0; col < rref.getNumCols(); col++) {
                if(HALMathUtil.round(rref.get(row, col), 7) != 0) {
                    isZeroRow = false;
                    break;
                }
            }
            reachedZeroRow |= isZeroRow;
            if(reachedZeroRow && !isZeroRow) {
                double[] nonZeroRow = rref.getRowArray(row);
                double[] zeroRow = rref.getRowArray(row-1);
                rref.setRow(row, zeroRow);
                rref.setRow(row-1, nonZeroRow);
            }
        }

        return rref;
    }

    public VectorSpace nullSpace() {
        QRFactorization qr = transpose().QRFactorization();
        int nullity = nullityFromQR(qr);

        if(nullity == 0) return new VectorSpace(MatrixSimple.zeroMatrix(getNumCols(), 1));

        MatrixSimple N = qr.Q.getCol(getNumCols()-nullity);
        for (int col = getNumCols()-nullity+1; col < getNumCols(); col++) {
            N = N.augment(qr.Q.getCol(col));
        }
        return new VectorSpace(N);
    }

    public VectorSpace columnSpace() {
        MatrixSimple rref = rref();
        List<MatrixSimple> columnSpaceBasis = new ArrayList<>();
        for (int row = 0; row < rref.getNumRows(); row++) {
            for (int col = 0; col < getNumCols(); col++) {
                if(HALMathUtil.floatingPointFix(rref.get(row, col)) == 1) {
                    MatrixSimple basisVector = MatrixSimple.zeroMatrix(getNumRows(), 1);
                    for (int originalRow = 0; originalRow < getNumRows(); originalRow++) {
                        basisVector.set(originalRow, 0, get(originalRow, col));
                    }
                    columnSpaceBasis.add(basisVector);
                    break;
                }
            }
        }

        if(columnSpaceBasis.isEmpty()) return new VectorSpace(new MatrixSimple[] {MatrixSimple.zeroMatrix(getNumCols(), 1)});
        return new VectorSpace(columnSpaceBasis.toArray(new MatrixSimple[0]));
    }

    public VectorSpace rowSpace() {
        return transpose().columnSpace();
    }

    public VectorSpace leftNullSpace() {
        return transpose().nullSpace();
    }

    private LUFactorization LUFactorization(MatrixSimple P) {
        MatrixSimple L = MatrixSimple.identityMatrix(this.getNumRows());
        MatrixSimple U = P.multiply(this);

        for (int pivot = 0; pivot < Math.min(getNumRows(), getNumCols()); pivot++) {
            double pivotVal = U.get(pivot, pivot);

            if(pivotVal == 0) {

                int nonZeroRowIdx = -1;
                for (int row = pivot; row < U.getNumRows(); row++) {
                    if(U.get(row, pivot) != 0) {
                        nonZeroRowIdx = row;
                        break;
                    }
                }
                if(nonZeroRowIdx == -1) {
                    continue;
                }

                MatrixSimple permutation = MatrixSimple.identityMatrix(getNumRows());
                double[] pivotRow = permutation.getRowArray(pivot);
                double[] nonZeroRow = permutation.getRowArray(nonZeroRowIdx);

                permutation.setRow(pivot, nonZeroRow);
                permutation.setRow(nonZeroRowIdx, pivotRow);

                return LUFactorization(permutation.multiply(P));
            }

            for (int row = pivot+1; row < getNumRows(); row++) {

                double elimVal = U.get(row, pivot) / pivotVal;

                //Generate elimination matrix and its inverse
                MatrixSimple E = MatrixSimple.identityMatrix(getNumRows());
                MatrixSimple Einv = MatrixSimple.identityMatrix(getNumRows());
                E.set(row, pivot, -elimVal);
                Einv.set(row, pivot, elimVal);

                U = E.multiply(U);
                L = L.multiply(Einv);
            }
        }

        return new LUFactorization(P, L, U);
    }

    public LUFactorization LUFactorization() {
        return LUFactorization(MatrixSimple.identityMatrix(this.getNumRows()));
    }

    public static MatrixSimple givensRotationMatrix(int size, int modifyRow, int setRow, double c, double s)  {
        math.ExceptionChecker.assertTrue(setRow >= 1, new ArithmeticException("Givens' rotation can't cancel anything on row 0."));
        math.ExceptionChecker.assertTrue(modifyRow < setRow, new ArithmeticException("Givens' rotation secondary modified row must be before the main modified row."));
        MatrixSimple givensRotationMatrix = MatrixSimple.identityMatrix(size);

        givensRotationMatrix.set(modifyRow, modifyRow, c);
        givensRotationMatrix.set(setRow, modifyRow, s);
        givensRotationMatrix.set(modifyRow, setRow, -s);
        givensRotationMatrix.set(setRow, setRow, c);

        return givensRotationMatrix;
    }

    public QRFactorization QRFactorization() {


        int rows = getNumRows();
        int cols = getNumCols();

        MatrixSimple Q = MatrixSimple.identityMatrix(rows);
        MatrixSimple R = copy();

        for (int col = 0; col < cols; col++) {
            for (int row = rows-1; row > col; row--) {
                double a = R.get(row,col);
                if(a == 0) continue; //If the value is already canceled, we're done.

                double b = R.get(row - 1, col);

                double r = Math.hypot(a, b);

                double c = b / r;
                double s = -a / r;

                MatrixSimple G = MatrixSimple.givensRotationMatrix(rows, row-1, row, c,s);

                R = G.multiply(R);
                Q = G.multiply(Q);
            }
        }

        Q = Q.transpose();

        return new QRFactorization(Q, R);
    }

    public double[] eigenvalues() {
        math.ExceptionChecker.assertTrue(this.isSquare(), new RuntimeException("Matrix must be square!"));

        MatrixSimple A = getEigenvalueDiagonalMatrix();

        double[] eigenvalues = new double[A.getNumRows()];
        for (int diag = 0; diag < A.getNumRows(); diag++) {
            eigenvalues[diag] = A.get(diag, diag);
        }

        return eigenvalues;
    }

    private MatrixSimple[] eigenvectorsFromLambda(MatrixSimple eigenvalues) {
        Set<Double> usedEigenvalues = new LinkedHashSet<>();
        List<MatrixSimple> eigenvectors = new ArrayList<>();
        for(int diag = 0; diag < eigenvalues.getNumRows(); diag++) {
            double lambda = eigenvalues.get(diag, diag);
            if(!usedEigenvalues.contains(lambda)) {
                MatrixSimple eigenvector = (this.subtract(MatrixSimple.identityMatrix(this.getNumRows()).multiply(lambda)))
                        .nullSpace().getBasisMatrix(); //The null space of A-(lambda)I
                eigenvectors.add(eigenvector);
            }
            usedEigenvalues.add(lambda);
        }
        return eigenvectors.toArray(new MatrixSimple[0]);
    }

    private MatrixSimple getEigenvalueDiagonalMatrix() {
        QRFactorization qrFactorization = this.QRFactorization();
        MatrixSimple A = qrFactorization.R.multiply(qrFactorization.Q);

        MatrixSimple identity = MatrixSimple.identityMatrix(qrFactorization.Q.getNumRows());

        while(!A.round(7).isDiagonal()) {

            double s = A.get(getNumRows()-1, getNumRows()-1);
            MatrixSimple smult = MatrixSimple.identityMatrix(A.getNumRows()).multiply(s);

            qrFactorization = A.subtract(smult).QRFactorization();

            A = qrFactorization.R.multiply(qrFactorization.Q).add(smult);

            //If matrix is already upper triangular
            if(qrFactorization.Q.equals(identity)) {
                for (int diag = 0; diag < A.getNumRows(); diag++) {
                    identity.set(diag, diag, A.get(diag, diag));
                }
                A = identity;
                break;
            }
        }
        return A;
    }

    public MatrixSimple[] eigenvectors() {
        math.ExceptionChecker.assertTrue(this.isSquare(), new RuntimeException("Matrix must be square!"));
        MatrixSimple eigenvalues = getEigenvalueDiagonalMatrix();
        return eigenvectorsFromLambda(eigenvalues);
    }

    public EigenData eigen() {
        math.ExceptionChecker.assertTrue(this.isSquare(), new RuntimeException("Matrix must be square!"));

        math.ExceptionChecker.assertTrue(this.isSquare(), new RuntimeException("Matrix must be square!"));

        MatrixSimple lambda = getEigenvalueDiagonalMatrix();
        MatrixSimple[] eigenvectors = eigenvectorsFromLambda(lambda);

        math.ExceptionChecker.assertTrue(eigenvectors.length > 0, new ArithmeticException("Something has gone wrong, there were no eigenvectors."));

        MatrixSimple eigenbasis = eigenvectors[0];
        for (int i = 1; i < eigenvectors.length; i++) {
            eigenbasis = eigenbasis.augment(eigenvectors[i]);
        }

        return new EigenData(lambda, eigenbasis);
    }

    public boolean isDiagonalizable() {
        return eigen().eigenbasis.nullity() == 0;
    }

    public Diagonalization diagonalize() {
        math.ExceptionChecker.assertTrue(isSquare(), new ArithmeticException("Matrix is not square."));
        //ExceptionChecker.assertTrue(isDiagonalizable(), new ArithmeticException("Matrix is not diagonalizable."));

        EigenData eigenData = eigen();

        MatrixSimple lambda = MatrixSimple.identityMatrix(eigenData.eigenvalues.length);
        for (int diag = 0; diag < lambda.getNumRows(); diag++) {
            lambda.set(diag, diag, eigenData.eigenvalues[diag]);
        }

        return new Diagonalization(eigenData.eigenbasis, lambda);
    }

    private static int nullityFromQR(QRFactorization qr) {
        int rows = qr.R.getNumRows();
        int cols = qr.R.getNumCols();

        int rank = 0;
        for (int row = rows-1; row >= 0; row--) {
            boolean isZeroRow = true;
            for (int col = 0; col < cols; col++) {
                if(HALMathUtil.round(qr.R.get(row, col),7) != 0) {
                    isZeroRow = false;
                    break;
                }
            }
            if(isZeroRow) rank++;
        }

        return rank;
    }

    /**
     * Calculates the rank of the matrix.
     *
     * @return The rank of the matrix.
     */
    public int rank() {
        return getNumCols()-nullity();
    }

    /**
     * Calculates the nullity of the matrix (uses the rank-nullity theorem).
     *
     * @return The nullity of the matrix.
     */
    public int nullity() {
        return nullityFromQR(transpose().QRFactorization());
    }

    public boolean isVector() {
        return getNumCols() == 1;
    }

    public math.VectorND toVector() {
        math.ExceptionChecker.assertTrue(isVector(), new ArithmeticException("matrix is not a vector."));
        return new math.VectorND(this);
    }

    public MatrixSimple round(int places) {
        math.ExceptionChecker.assertTrue(places >= 0, new ArithmeticException("You cannot have negative decimal places."));

        MatrixSimple rounded = copy();

        int rows = getNumRows();
        int cols = getNumCols();

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                rounded.set(row, col, HALMathUtil.round(get(row, col), places));
            }
        }

        return rounded;
    }

    /**
     * Creates an identity matrix of the given size.
     *
     * @param size The size of the identity matrix to create (used for width and height, as identity matrices are square).
     * @return An identity matrix of the given size.
     * @throws RuntimeException Throws this exception if the given size is 1 or less, as identity matrices must be at least 2x2.
     */

    @NotNull
    public static MatrixSimple identityMatrix(int size) {
        math.ExceptionChecker.assertTrue(size > 0, new RuntimeException("An identity matrix must be at least 1x1"));
        double[][] identityVals = new double[size][size];
        for (int diag = 0; diag < size; diag++) {
            identityVals[diag][diag] = 1;
        }

        return new MatrixSimple(identityVals);
    }

    /**
     * Creates a zero matrix of the given size.
     *
     * @param rows The number of rows in the zero matrix.
     * @param cols The number of columns in the zero matrix.
     * @return A zero matrix of the given size.
     * @throws RuntimeException Throws this exception if the number of rows or columns given was negative.
     */
    @NotNull
    public static MatrixSimple zeroMatrix(int rows, int cols) {
        math.ExceptionChecker.assertTrue(rows > 0, new RuntimeException("Can't have 0 or negative amount of rows."));
        math.ExceptionChecker.assertTrue(cols > 0, new RuntimeException("Can't have 0 or negative amount of columns."));

        return new MatrixSimple(new double[rows][cols]);
    }

    /**
     * Creates a ones matrix of the given size.
     *
     * @param rows The number of rows in the ones matrix.
     * @param cols The number of columns in the ones matrix.
     * @return A ones matrix of the given size.
     * @throws RuntimeException Throws this exception if the number of rows or columns given was negative.
     */
    @NotNull
    public static MatrixSimple onesMatrix(int rows, int cols) {
        math.ExceptionChecker.assertTrue(rows > 0, new RuntimeException("Can't have 0 or negative amount of rows."));
        math.ExceptionChecker.assertTrue(cols > 0, new RuntimeException("Can't have 0 or negative amount of columns."));

        double[][] onesVals = new double[rows][cols];
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                onesVals[row][col] = 1;
            }
        }

        return new MatrixSimple(onesVals);
    }

    public static class LUFactorization {
        public final MatrixSimple P, L, U;
        protected LUFactorization(MatrixSimple P, MatrixSimple L, MatrixSimple U) {
            this.P = P;
            this.L = L;
            this.U = U;
        }
    }

    public static class QRFactorization {
        public final MatrixSimple Q, R;
        protected QRFactorization(MatrixSimple Q, MatrixSimple R) {
            this.Q = Q;
            this.R = R;
        }
    }

    public static class EigenData {
        public final double[] eigenvalues;
        public final MatrixSimple[] eigenvectors;
        public final MatrixSimple eigenbasis;

        private final String eigenvectorString;
        protected EigenData(MatrixSimple diagEigenvalueMatrix, @NotNull MatrixSimple eigenVectorMatrix) {
            math.ExceptionChecker.assertTrue(diagEigenvalueMatrix.round(7).isDiagonal(), new ArithmeticException("Matrix not diagonal"));
            math.ExceptionChecker.assertTrue(eigenVectorMatrix.isSquare(), new ArithmeticException("Matrix not square"));

            this.eigenbasis = eigenVectorMatrix;

            this.eigenvalues = new double[diagEigenvalueMatrix.getNumRows()];
            for (int diag = 0; diag < diagEigenvalueMatrix.getNumRows(); diag++) {
                eigenvalues[diag] = diagEigenvalueMatrix.get(diag, diag);
            }

            eigenvectorString = eigenVectorMatrix.toString();
            eigenvectors = new MatrixSimple[eigenVectorMatrix.getNumCols()];
            for (int col = 0; col < eigenVectorMatrix.getNumCols(); col++) {
                MatrixSimple vector = eigenVectorMatrix.crop(0, eigenVectorMatrix.getNumRows(), col, col+1);
                eigenvectors[col] = vector;
            }
        }

        @Override
        public String toString() {
            return "Eigenvalues "+Arrays.toString(eigenvalues)+"\nEigenvectors:\n"+eigenvectorString;
        }
    }

    //TODO BROKEN
    public static class Diagonalization {
        public final MatrixSimple X, lambda;
        private final MatrixSimple Xinv;
        private Diagonalization(MatrixSimple X, MatrixSimple lambda) {
            math.ExceptionChecker.assertTrue(lambda.isDiagonal(), new ArithmeticException("Matrix must be diagonal"));

            /*
            for (MatrixSimple vector : X.getCols()) {
                double max = -Double.MAX_VALUE;
                for (int row = 0; row < vector.getNumRows(); row++) {
                    max = Math.max(vector.get(row, 0), max);
                }
                for (int row = 0; row < vector.getNumRows(); row++) {
                    vector.set(row, 0, vector.get(row, 0)/max);
                }
                System.out.println(vector);
                System.out.println();
            }

            System.out.println("-------------------------------");

            System.out.println(X);
            System.out.println();
            System.out.println(lambda);
*/
            this.X = X;
            this.lambda = lambda;
            this.Xinv = X.invert();
        }

        public MatrixSimple fastPow(int exponent) {
            math.ExceptionChecker.assertTrue(exponent > 0, new ArithmeticException("matrix exponent must be an integer"));

            MatrixSimple lambdaPow = lambda.copy();
            for (int diag = 0; diag < lambdaPow.getNumRows(); diag++) {
                lambdaPow.set(diag, diag, pow(lambda.get(diag, diag), exponent));
            }

            return X.multiply(lambdaPow).multiply(Xinv);
        }
    }

    @Override
    public String toString() {
        StringBuilder output = new StringBuilder();
        for (double[] row : vals) {
            output.append(Arrays.toString(row));
            output.append('\n');
        }
        output.deleteCharAt(output.length() - 1);
        return output.toString();
    }

    @Override
    public int hashCode() {
        int hashCode = 17;
        for (int row = 0; row < this.getNumRows(); row++) {
            for (int col = 0; col < this.getNumCols(); col++) {
                long bits = Double.doubleToLongBits(this.get(row, col));
                int hash = (int) (bits ^ (bits >>> 32));
                hashCode = 37 * hashCode + hash;
            }
        }
        return hashCode;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof MatrixSimple) {
            MatrixSimple matrix = (MatrixSimple) obj;
            if (this.getNumRows() != matrix.getNumRows() || this.getNumCols() != matrix.getNumCols()) {
                return false;
            }

            boolean equal = true;
            for (int row = 0; row < this.getNumRows(); row++) {
                for (int col = 0; col < this.getNumCols(); col++) {
                    equal &= matrix.get(row, col) == this.get(row, col);
                }
            }
            return equal;
        }
        return false;
    }
}
