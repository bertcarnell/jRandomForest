package randomForestConnection;

/// <summary>
/// Utilities for the Random Forest class
/// </summary>
public class RandomForestUtilities {
    /// <summary>
    /// Create an array and fill it with a value
    /// </summary>
    /// <typeparam name="T">Type</typeparam>
    /// <param name="length">Length of the array</param>
    /// <param name="value">value to fill the array with</param>
    /// <returns>An array of Type T</returns>
    public static T[] CreateAndFill<T>(int length, T value) {
        T[] result = new T[length];
        for (int i = 0; i < length; i++) {
            result[i] = value;
        }
        return result;
    }
}
