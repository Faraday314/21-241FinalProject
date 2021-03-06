import org.jetbrains.annotations.Contract;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.*;

/**
 * A class for doing mathematical operations on arrays. Basically numpy, but less good.
 *
 * @author Cole Savage, Level Up
 * @since 1.0.0
 * @version 1.0.0
 *
 * Creation Date: 10/19/20.
 */
@SuppressWarnings("unused")
public class FakeNumpy {

    private static final double FLOATING_POINT_FIXER_CONSTANT = 1e9;

    /**
     * Private default constructor to make class basically static.
     */
    private FakeNumpy() {}

    /**
     * Finds the maximum value of an array.
     *
     * @param array The input array.
     * @param <T> The element type of the array.
     * @return The maximum of the array.
     */
    public static <T extends Comparable<? super T>> @Nullable T max(T @NotNull [] array) {
        if(array.length == 0) {
            return null;
        }
        return max(array,array.length);
    }

    /**
     * Finds the maximum value of an array.
     *
     * @param array The input array.
     * @param n The maximum index of the elements in the array that will be compared.
     * @param <T> The element type of the array.
     * @return The maximum of the array.
     */
    private static <T extends Comparable<? super T>> T max(T[] array, int n) {
        if(n == 1) {
            return array[0];
        }
        T currentMax = max(array, n-1);
        return array[n-1].compareTo(currentMax) > 0 ? array[n-1] : currentMax;
    }

    /**
     * Finds the maximum value of an array of doubles.
     *
     * @param array The input array.
     * @return The maximum of the array.
     */
    public static double max(double @NotNull [] array) {
        if(array.length == 0) {
            return 0;
        }
        return max(array,array.length);
    }

    /**
     * Finds the maximum value of an array of doubles.
     *
     * @param array The input array.
     * @param n The maximum index of the elements in the array that will be compared.
     * @return The maximum of the array.
     */
    private static double max(double[] array, int n) {
        if(n == 1) {
            return array[0];
        }
        return Math.max(array[n-1], max(array, n-1));
    }

    /**
     * Finds the maximum value of an array of integers.
     *
     * @param array The input array.
     * @return The maximum of the array.
     */
    public static int max(int @NotNull [] array) {
        if(array.length == 0) {
            return 0;
        }
        return max(array,array.length);
    }

    /**
     * Finds the maximum value of an array of integers.
     *
     * @param array The input array.
     * @param n The maximum index of the elements in the array that will be compared.
     * @return The maximum of the array.
     */
    private static int max(int[] array, int n) {
        if(n == 1) {
            return array[0];
        }
        return Math.max(array[n-1], max(array, n-1));
    }

    /**
     * Finds the maximum value of an array of floats.
     *
     * @param array The input array.
     * @return The maximum of the array.
     */
    public static float max(float @NotNull [] array) {
        if(array.length == 0) {
            return 0;
        }
        return max(array,array.length);
    }

    /**
     * Finds the maximum value of an array of floats.
     *
     * @param array The input array.
     * @param n The maximum index of the elements in the array that will be compared.
     * @return The maximum of the array.
     */
    private static float max(float[] array, int n) {
        if(n == 1) {
            return array[0];
        }
        return Math.max(array[n-1], max(array, n-1));
    }

    /**
     * Finds the maximum value of an array of longs.
     *
     * @param array The input array.
     * @return The maximum of the array.
     */
    public static long max(long @NotNull [] array) {
        if(array.length == 0) {
            return 0;
        }
        return max(array,array.length);
    }

    /**
     * Finds the maximum value of an array of longs.
     *
     * @param array The input array.
     * @param n The maximum index of the elements in the array that will be compared.
     * @return The maximum of the array.
     */
    private static long max(long[] array, int n) {
        if(n == 1) {
            return array[0];
        }
        return Math.max(array[n-1], max(array, n-1));
    }

    /**
     * Finds the maximum value of an array of shorts.
     *
     * @param array The input array.
     * @return The maximum of the array.
     */
    public static short max(short @NotNull [] array) {
        if(array.length == 0) {
            return 0;
        }
        return max(array,array.length);
    }

    /**
     * Finds the maximum value of an array of shorts.
     *
     * @param array The input array.
     * @param n The maximum index of the elements in the array that will be compared.
     * @return The maximum of the array.
     */
    private static short max(short[] array, int n) {
        if(n == 1) {
            return array[0];
        }
        return (short) Math.max(array[n-1], max(array, n-1));
    }

    /**
     * Finds the minimum value of an array.
     *
     * @param array The input array.
     * @param <T> The element type of the array.
     * @return The minimum of the array.
     */
    
    public static <T extends Comparable<? super T>> @Nullable T min(T @NotNull [] array) {
        if(array.length == 0) {
            return null;
        }
        return min(array,array.length);
    }

    /**
     * Finds the minimum value of an array.
     *
     * @param array The input array.
     * @param n The maximum index of the elements in the array that will be compared.
     * @param <T> The element type of the array.
     * @return The minimum of the array.
     */
    private static <T extends Comparable<? super T>> T min(T[] array, int n) {
        if(n == 1) {
            return array[0];
        }
        T currentMin = min(array, n-1);
        return array[n-1].compareTo(currentMin) < 0 ? array[n-1] : currentMin;
    }

    /**
     * Finds the minimum value of an array of doubles.
     *
     * @param array The input array.
     * @return The minimum of the array.
     */
    public static double min(double @NotNull [] array) {
        if(array.length == 0) {
            return 0;
        }
        return min(array,array.length);
    }

    /**
     * Finds the minimum value of an array of doubles.
     *
     * @param array The input array.
     * @param n The maximum index of the elements in the array that will be compared.
     * @return The minimum of the array.
     */
    private static double min(double[] array, int n) {
        if(n == 1) {
            return array[0];
        }
        return Math.min(array[n-1], min(array, n-1));
    }

    /**
     * Finds the minimum value of an array of integers.
     *
     * @param array The input array.
     * @return The minimum of the array.
     */
    public static int min(int @NotNull [] array) {
        if(array.length == 0) {
            return 0;
        }
        return min(array,array.length);
    }

    /**
     * Finds the minimum value of an array of integers.
     *
     * @param array The input array.
     * @param n The maximum index of the elements in the array that will be compared.
     * @return The minimum of the array.
     */
    private static int min(int[] array, int n) {
        if(n == 1) {
            return array[0];
        }
        return Math.min(array[n-1], min(array, n-1));
    }

    /**
     * Finds the minimum value of an array of floats.
     *
     * @param array The input array.
     * @return The minimum of the array.
     */
    public static float min(float @NotNull [] array) {
        if(array.length == 0) {
            return 0;
        }
        return min(array,array.length);
    }

    /**
     * Finds the minimum value of an array of floats.
     *
     * @param array The input array.
     * @param n The maximum index of the elements in the array that will be compared.
     * @return The minimum of the array.
     */
    private static float min(float[] array, int n) {
        if(n == 1) {
            return array[0];
        }
        return Math.min(array[n-1], min(array, n-1));
    }

    /**
     * Finds the minimum value of an array of longs.
     *
     * @param array The input array.
     * @return The minimum of the array.
     */
    public static long min(long @NotNull [] array) {
        if(array.length == 0) {
            return 0;
        }
        return min(array,array.length);
    }

    /**
     * Finds the minimum value of an array of longs.
     *
     * @param array The input array.
     * @param n The maximum index of the elements in the array that will be compared.
     * @return The minimum of the array.
     */
    private static long min(long[] array, int n) {
        if(n == 1) {
            return array[0];
        }
        return Math.min(array[n-1], min(array, n-1));
    }

    /**
     * Finds the minimum value of an array of shorts.
     *
     * @param array The input array.
     * @return The minimum of the array.
     */
    public static short min(short @NotNull [] array) {
        if(array.length == 0) {
            return 0;
        }
        return min(array,array.length);
    }

    /**
     * Finds the minimum value of an array of shorts.
     *
     * @param array The input array.
     * @param n The maximum index of the elements in the array that will be compared.
     * @return The minimum of the array.
     */
    private static short min(short[] array, int n) {
        if(n == 1) {
            return array[0];
        }
        return (short) Math.min(array[n-1], min(array, n-1));
    }

    /**
     * Slices an array from startIdx to endIdx (inclusive).
     *
     * @param array The input array.
     * @param startIdx The starting index for the slice.
     * @param endIdx The ending index for the slice.
     * @param <T> The element type of the array.
     * @return The slice of the array that goes from startIdx to endIdx (inclusive).
     */
    @Contract(value = "_, _, _ -> new", pure = true)
    public static <T> T @NotNull [] slice(T[] array, int startIdx, int endIdx) {
        return Arrays.copyOfRange(array, startIdx, endIdx + 1);
    }

    /**
     * Slices a double array from startIdx to endIdx (inclusive).
     *
     * @param array The input array.
     * @param startIdx The starting index for the slice.
     * @param endIdx The ending index for the slice.
     * @return The slice of the array that goes from startIdx to endIdx (inclusive).
     */
    @Contract(value = "_, _, _ -> new", pure = true)
    public static double @NotNull [] slice(double[] array, int startIdx, int endIdx) {
        return Arrays.copyOfRange(array, startIdx, endIdx + 1);
    }

    /**
     * Slices an integer array from startIdx to endIdx (inclusive).
     *
     * @param array The input array.
     * @param startIdx The starting index for the slice.
     * @param endIdx The ending index for the slice.
     * @return The slice of the array that goes from startIdx to endIdx (inclusive).
     */
    @Contract(value = "_, _, _ -> new", pure = true)
    public static int @NotNull [] slice(int[] array, int startIdx, int endIdx) {
        return Arrays.copyOfRange(array, startIdx, endIdx + 1);
    }

    /**
     * Slices a float array from startIdx to endIdx (inclusive).
     *
     * @param array The input array.
     * @param startIdx The starting index for the slice.
     * @param endIdx The ending index for the slice.
     * @return The slice of the array that goes from startIdx to endIdx (inclusive).
     */
    @Contract(value = "_, _, _ -> new", pure = true)
    public static float @NotNull [] slice(float[] array, int startIdx, int endIdx) {
        return Arrays.copyOfRange(array, startIdx, endIdx + 1);
    }

    /**
     * Slices a long array from startIdx to endIdx (inclusive).
     *
     * @param array The input array.
     * @param startIdx The starting index for the slice.
     * @param endIdx The ending index for the slice.
     * @return The slice of the array that goes from startIdx to endIdx (inclusive).
     */
    @Contract(value = "_, _, _ -> new", pure = true)
    public static long @NotNull [] slice(long[] array, int startIdx, int endIdx) {
        return Arrays.copyOfRange(array, startIdx, endIdx + 1);
    }

    /**
     * Slices a short array from startIdx to endIdx (inclusive).
     *
     * @param array The input array.
     * @param startIdx The starting index for the slice.
     * @param endIdx The ending index for the slice.
     * @return The slice of the array that goes from startIdx to endIdx (inclusive).
     */
    @Contract(value = "_, _, _ -> new", pure = true)
    public static short @NotNull [] slice(short[] array, int startIdx, int endIdx) {
        return Arrays.copyOfRange(array, startIdx, endIdx + 1);
    }

    /**
     * Multiplies every element in a double array by a constant double.
     *
     * @param array The input array.
     * @param multiplier The constant to multiply by.
     */
    public static void multiply(double @NotNull [] array, double multiplier) {
        for(int i = 0; i < array.length; i++) {
            array[i] = array[i]*multiplier;
        }
    }

    /**
     * Multiplies every element in an integer array by a constant integer.
     *
     * @param array The input array.
     * @param multiplier The constant to multiply by.
     */
    public static void multiply(int @NotNull [] array, int multiplier) {
        for(int i = 0; i < array.length; i++) {
            array[i] = array[i]*multiplier;
        }
    }

    /**
     * Multiplies every element in a float array by a constant float.
     *
     * @param array The input array.
     * @param multiplier The constant to multiply by.
     */
    public static void multiply(float @NotNull [] array, float multiplier) {
        for(int i = 0; i < array.length; i++) {
            array[i] = array[i]*multiplier;
        }
    }

    /**
     * Multiplies every element in a long array by a constant long.
     *
     * @param array The input array.
     * @param multiplier The constant to multiply by.
     */
    public static void multiply(long @NotNull [] array, long multiplier) {
        for(int i = 0; i < array.length; i++) {
            array[i] = array[i]*multiplier;
        }
    }

    /**
     * Multiplies every element in an integer array by a constant double and rounds it to an integer.
     *
     * @param array The input array.
     * @param multiplier The constant to multiply by.
     */
    public static void multiply(int @NotNull [] array, double multiplier) {
        for(int i = 0; i < array.length; i++) {
            array[i] = (int) Math.round(array[i]*multiplier);
        }
    }

    /**
     * Multiplies every element in a float array by a constant double.
     *
     * @param array The input array.
     * @param multiplier The constant to multiply by.
     */
    public static void multiply(float @NotNull [] array, double multiplier) {
        for(int i = 0; i < array.length; i++) {
            array[i] = array[i]*(float) multiplier;
        }
    }

    /**
     * Multiplies every element in a long array by a constant double.
     *
     * @param array The input array.
     * @param multiplier The constant to multiply by.
     */
    public static void multiply(long @NotNull [] array, double multiplier) {
        for(int i = 0; i < array.length; i++) {
            array[i] = Math.round(array[i]*multiplier);
        }
    }

    /**
     * Divides every element in a double array by a constant double.
     *
     * @param array The input array.
     * @param multiplier The constant to divide by.
     */
    public static void divide(double[] array, double multiplier) {
        ExceptionChecker.assertNotEqual(multiplier,0.0,new ArithmeticException("You can't divide by zero."));
        multiply(array, 1.0/multiplier);
    }

    /**
     * Divides every element in an integer array by a constant integer.
     *
     * @param array The input array.
     * @param multiplier The constant to divide by.
     */
    public static void divide(int[] array, int multiplier) {
        ExceptionChecker.assertNotEqual(multiplier,0,new ArithmeticException("You can't divide by zero."));
        multiply(array, 1.0/multiplier);
    }

    /**
     * Divides every element in a float array by a constant float.
     *
     * @param array The input array.
     * @param multiplier The constant to divide by.
     */
    public static void divide(float[] array, float multiplier) {
        ExceptionChecker.assertNotEqual(multiplier,0.0f,new ArithmeticException("You can't divide by zero."));
        multiply(array, 1.0/multiplier);
    }

    /**
     * Divides every element in a long array by a constant long.
     *
     * @param array The input array.
     * @param multiplier The constant to divide by.
     */
    public static void divide(long[] array, long multiplier) {
        ExceptionChecker.assertNotEqual(multiplier,0L,new ArithmeticException("You can't divide by zero."));
        multiply(array, 1.0/multiplier);
    }

    /**
     * Divides every element in a integer array by a constant double and rounds it to an integer.
     *
     * @param array The input array.
     * @param multiplier The constant to divide by.
     */
    public static void divide(int[] array, double multiplier) {
        ExceptionChecker.assertNotEqual(multiplier,0.0,new ArithmeticException("You can't divide by zero."));
        multiply(array, 1.0/multiplier);
    }

    /**
     * Divides every element in a float array by a constant double.
     *
     * @param array The input array.
     * @param multiplier The constant to divide by.
     */
    public static void divide(float[] array, double multiplier) {
        ExceptionChecker.assertNotEqual(multiplier,0.0,new ArithmeticException("You can't divide by zero."));
        multiply(array, 1.0/multiplier);
    }

    /**
     * Divides every element in a long array by a constant double.
     *
     * @param array The input array.
     * @param multiplier The constant to divide by.
     */
    public static void divide(long[] array, double multiplier) {
        ExceptionChecker.assertNotEqual(multiplier,0.0,new ArithmeticException("You can't divide by zero."));
        multiply(array, 1.0/multiplier);
    }

    /**
     * Takes the absolute value of every element in a double array.
     *
     * @param array The input array.
     * @return The absolute value of the input array.
     */
    @Contract(pure = true)
    public static double @NotNull [] abs(double @NotNull [] array) {
        double[] output = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            output[i] = Math.abs(array[i]);
        }
        return output;
    }

    /**
     * Takes the absolute value of every element in an integer array.
     *
     * @param array The input array.
     * @return The absolute value of the input array.
     */
    @Contract(pure = true)
    public static int @NotNull [] abs(int @NotNull [] array) {
        int[] output = new int[array.length];
        for (int i = 0; i < array.length; i++) {
            output[i] = Math.abs(array[i]);
        }
        return output;
    }

    /**
     * Takes the absolute value of every element in a float array.
     *
     * @param array The input array.
     * @return The absolute value of the input array.
     */
    @Contract(pure = true)
    public static float @NotNull [] abs(float @NotNull [] array) {
        float[] output = new float[array.length];
        for (int i = 0; i < array.length; i++) {
            output[i] = Math.abs(array[i]);
        }
        return output;
    }

    /**
     * Takes the absolute value of every element in a long array.
     *
     * @param array The input array.
     * @return The absolute value of the input array.
     */
    @Contract(pure = true)
    public static long @NotNull [] abs(long @NotNull [] array) {
        long[] output = new long[array.length];
        for (int i = 0; i < array.length; i++) {
            output[i] = Math.abs(array[i]);
        }
        return output;
    }

    public static <T> boolean checkForDuplicates(T @NotNull [] array) {
        Set<T> set = new HashSet<>();
        for(T element : array) {
            if(set.contains(element)) {
                return true;
            }
            set.add(element);
        }
        return false;
    }

    public static <T> T @NotNull [] removeDuplicates(T[] array)
    {
        Set<T> set = new LinkedHashSet<>(Arrays.asList(array));

        List<T> lst = new ArrayList<>(set);

        T[] arrOut = slice(array.clone(),0,lst.size()-1);
        for(int i = 0; i < lst.size(); i++) {
            arrOut[i] = lst.get(i);
        }
        return arrOut;
    }

    public static double mod(double x, double modulus) {
        return (x % modulus + modulus) % modulus;
    }

    @Contract(pure = true)
    public static double @NotNull [] add(double @NotNull [] list1, double @NotNull [] list2) {
        if(list1.length != list2.length) {
            throw new ArithmeticException("Arrays are different sizes, can't be subtracted");
        }

        double[] output = new double[list1.length];
        for(int i = 0; i < list1.length; i++) {
            output[i] = list1[i] + list2[i];
        }
        return output;
    }

    @Contract(pure = true)
    public static int @NotNull [] add(int @NotNull [] list1, int @NotNull [] list2) {
        if(list1.length != list2.length) {
            throw new ArithmeticException("Arrays are different sizes, can't be subtracted");
        }

        int[] output = new int[list1.length];
        for(int i = 0; i < list1.length; i++) {
            output[i] = list1[i] + list2[i];
        }
        return output;
    }

    //list1 - list2
    public static double @NotNull [] subtract(double[] list1, double @NotNull [] list2) {
        double[] list2cpy = list2.clone();
        multiply(list2cpy, -1);
        return add(list2cpy, list1);
    }

    public static int @NotNull [] subtract(int[] list1, int @NotNull [] list2) {
        int[] list2cpy = list2.clone();
        multiply(list2cpy, -1);
        return add(list2cpy, list1);
    }

    public static double @NotNull [] absdiff(double[] list1, double[] list2) {
        return abs(subtract(list1,list2));
    }

    public static int @NotNull [] absdiff(int[] list1, int[] list2) {
        return abs(subtract(list1,list2));
    }

    @Contract(pure = true)
    public static double average(double @NotNull [] list) {
        double sum = 0;
        for(double d : list) {
            sum += d;
        }
        return sum/list.length;
    }

    @Contract(pure = true)
    public static int @NotNull [] round(double @NotNull [] list) {
        int[] output = new int[list.length];
        for (int i = 0; i < list.length; i++) {
            output[i] = (int) Math.round(list[i]);
        }
        return output;
    }

    @Contract("_ -> param1")
    public static double[] floatingPointFix(double @NotNull [] list) {
        for (int i = 0; i < list.length; i++) {
            list[i] = Math.round(list[i]*FLOATING_POINT_FIXER_CONSTANT)/FLOATING_POINT_FIXER_CONSTANT;
        }
        return list;
    }

}