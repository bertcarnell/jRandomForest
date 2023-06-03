package randomForestConnection;

import java.util.List;

public class ConfusionMatrix
{
    /// <summary>
    /// The Confusion Matrix where the rows are the true classes and the columns are the predicted classes
    /// </summary>
    public int[][] confusionMatrix;
    /// <summary>
    /// The error rate of the classifier for each truth class
    /// </summary>
    public double[] error_rate;
    /// <summary>
    /// The class labels for each row and column
    /// </summary>
    public int[] classes;

    /// <summary>
    /// Creates a confusion matrix from the truth classes and predicted classes
    /// </summary>
    /// <param name="truthClass">Truth classes for each observation</param>
    /// <param name="predClass">Predicted classes for each observation</param>
    /// <returns>The confusion matrix</returns>
    public static ConfusionMatrix CreateConfusion(int[] truthClass, int[] predClass)
    {
        if (truthClass.length != predClass.length)
            throw new ArgumentException("Truth and prediction are not equal lengths", "truthClass");
        // check that all classes in predClass are found in truthClass if they aren't, add them
        //   This can happen when creating a confusion matrix for validation options
        List<Integer> availableClasses = truthClass.ToList();
        for (Integer a : predClass) {
            if (!availableClasses.contains(a))
            {
                availableClasses.add(a);
            }
        }
        ConfusionMatrix result = new ConfusionMatrix();
        result.classes = availableClasses.Distinct().ToArray();
        int nclasses = result.classes.length;
        result.confusionMatrix = new int[nclasses][nclasses];
        result.error_rate = new double[nclasses];

        int tempTruth, tempPred;
        for (int i = 0; i < nclasses; i++) // rows
        {
            tempTruth = result.classes[i];
            for (int j = 0; j < nclasses; j++) // cols
            {
                tempPred = result.classes[j];
                result.confusionMatrix[i][j] = 0;
                for (int k = 0; k < truthClass.length; k++) // over the input vectors
                {
                    if (tempTruth == truthClass[k] && tempPred == predClass[k])
                        result.confusionMatrix[i][j] += 1;
                }
            }
        }
        double[] rowsums = new MatrixVector<int>(result.confusionMatrix).rowSum();
        for (int i = 0; i < nclasses; i++)
        {
            // initialize the error_rate to zero
            result.error_rate[i] = 0.0;
            // check for divide by zero and set the error rate
            if (rowsums[i] != 0.0)
                result.error_rate[i] = 1.0 - Convert.ToDouble(result.confusionMatrix[i][i]) / Convert.ToDouble(rowsums[i]);
            // else the error_rate has already been initialized to zero
            // this presents doing this...
            //   else result.error_rate[i] = 0.0;
            // which is never executed and therefore not covered by tests
        }
        return result;
    }
}