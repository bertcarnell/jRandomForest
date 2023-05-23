/*
 * Copyright 2013 Rob Carnell
 * 
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Data;
using System.Linq.Expressions;
using System.Diagnostics;

namespace RandomForestConnection
{
    /// <summary>
    /// RandomForest class
    /// </summary>
    [Serializable]
    public class RandomForest
    {
        #region member variables
        // member variables, modeled after the R class RandomForest
        private MatrixVector<double> m_x, m_xtest;
        private MatrixVector<double> m_votes, m_votes_test, m_xbestsplit, m_error_rate_test, m_proximity_test;
        private MatrixVector<int> m_treemap, m_nodestatus, m_bestvar, m_nodepred, m_inbag;
        private int m_ntree, m_mtry, m_nclass, m_maxcat, m_nrnodes;
        private ConfusionMatrix m_confusion, m_confusion_test;
        private List<int[]> m_xlevels;
        private string m_type;
        private int[] m_predicted, m_y, m_classes, m_ndbigtree,
            m_ncat, m_predicted_test;
        private double[,] m_importance, m_importanceSD, m_localImportance, m_err_rate, m_proximity;
        // refactored member variables
        private double[] m_oob_times, m_pid, m_cutoff;
        private BiDirectionalDictionary<int, int> m_classCorrespondence;
        private Random m_oRandom;
        #endregion member variables

        #region accessors
        /// <summary>
        /// Get the Predicted Classes of the Training Data
        /// </summary>
        public int[] PredictedClass { get { return m_predicted; } }
        /// <summary>
        /// Get the Predicted Class(es) of the unknown sample(s)
        /// </summary>
        public int[] PredictedClassTest { get { return m_predicted_test; } }
        /// <summary>
        /// Get the votes for each class in the unknown sample
        /// </summary>
        public double[] VotesTest { get { return m_votes_test.values; } }
        /// <summary>
        /// Get the confusion matrix for the training data
        /// </summary>
        public ConfusionMatrix ConfusionMatrix
        {
            get
            {
                if (m_confusion.confusionMatrix != null)
                    return m_confusion;
                else
                    return new ConfusionMatrix();
            }
        }
        #endregion accessors
        
        #region helper functions
        private static int[] createOptions(bool addclass, bool ReturnImportance, bool LocalImportance,
            bool ReturnProximity, bool ProximityOnlyOutOfBag, bool DoTrace, bool KeepForest,
            bool SampleWithReplacement, bool Stratify, bool KeepInBagSamples)
        {
            int[] Options = new int[10];
            Options[0] = Convert.ToInt32(addclass);
            Options[1] = Convert.ToInt32(ReturnImportance);
            Options[2] = Convert.ToInt32(LocalImportance);
            Options[3] = Convert.ToInt32(ReturnProximity);
            Options[4] = Convert.ToInt32(ProximityOnlyOutOfBag);
            Options[5] = Convert.ToInt32(DoTrace);
            Options[6] = Convert.ToInt32(KeepForest);
            Options[7] = Convert.ToInt32(SampleWithReplacement);
            Options[8] = Convert.ToInt32(Stratify);
            Options[9] = Convert.ToInt32(KeepInBagSamples);
            return Options;
        }
        #endregion helper functions

        #region constructors
        /// <summary>
        /// 
        /// </summary>
        public RandomForest()
        {
            m_oRandom = new Random(1976);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="seed1"></param>
        /// <param name="seed2"></param>
        public RandomForest(int seed1, int seed2)
        {
            m_oRandom = new Random(Math.Abs(seed1 * seed2));
        }
        #endregion constructors

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="tree"></param>
        /// <param name="mtry"></param>
        /// <param name="bReturnImportance"></param>
        /// <param name="bReturnProximity"></param>
        /// <param name="xtest"></param>
        /// <param name="ytest"></param>
        /// <param name="indexCategoricalVariable"></param>
        /// <param name="classwt"></param>
        public void BalancedClassification(double[,] x, int[] y, int tree = 500, int mtry = -1,
            bool bReturnImportance = true, bool bReturnProximity = true,
            double[,] xtest = null, int[] ytest = null, int[] indexCategoricalVariable = null,
            double[] classwt = null)
        {
            // find the smallest class size and use that as the sample size
            List<int> distincty = y.Distinct().ToList();
            List<int> sampsize = new List<int>(distincty.Count);
            for (int i = 0; i < distincty.Count; i++)
                sampsize.Add(y.Count(a => a == distincty[i]));

            // TODO:  Worry about the order of the sampsize.  Is the distinct() method guaranteed to have the same order?  What about CRandomForest?
            Classification(x, y, tree, mtry, bReturnImportance, bReturnProximity, xtest, ytest, indexCategoricalVariable, classwt, y, sampsize.ToArray());
        }

        /// <summary>
        /// Perform a Random Forest classification
        /// </summary>
        /// <param name="x">The data matrix where the rows are observations and the columns are parameters</param>
        /// <param name="y">The classes of the observations</param>
        /// <param name="ntree">The number of trees to grow</param>
        /// <param name="mtry">The number of variables to choose from at each split</param>
        /// <param name="bReturnImportance">Should importance of predictors be assessed?</param>
        /// <param name="bReturnProximity">Should proximity measure among the rows be calculated?</param>
        /// <param name="xtest">Test data set to be classified</param>
        /// <param name="ytest">Classes of the test data, if known for validation testing</param>
        /// <param name="indexCategoricalVariable">Indices of the predictor variables that are considered categorical (zero based)</param>
        /// <param name="classwt">Prior probabilities for the classes</param>
        /// <param name="strata">The strata for stratified sampling of each observation</param>
        /// <param name="sampsize">The sample sizes for each strata</param>
        public void Classification(double[,] x, int[] y, int ntree=500, int mtry = -1, 
            bool bReturnImportance = true, bool bReturnProximity = true,
            double[,] xtest = null, int[] ytest = null, int[] indexCategoricalVariable = null,
            double[] classwt = null, int[] strata = null, int[] sampsize=null)
        {
            /*
             * R function call
             * 
             * function(x, y=NULL,  xtest=NULL, ytest=NULL, ntree=500,
             * mtry=if (!is.null(y) && !is.factor(y)) max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))),
             * replace=TRUE, classwt=NULL, cutoff, strata,
             * sampsize = if (replace) nrow(x) else ceiling(.632*nrow(x)),
             * nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
             * maxnodes=NULL,importance=FALSE, localImp=FALSE, nPerm=1,
             * proximity, oob.prox=proximity,
             * norm.votes=TRUE, do.trace=FALSE,
             * keep.forest=!is.null(y) && is.null(xtest), corr.bias=FALSE,
             * keep.inbag=FALSE, ...)
             */

            // data should have rows that are shots and columns that are metrics
            // TODO:  worry about nulls in data

            // Variables needed up front - not following the original R code
            if (x == null)
                throw new Exception("X may not be null");
            m_x = new MatrixVector<double>(x);
            int n = m_x.nrow;
            int p = m_x.ncol;
            if (y != null)
            {
                m_classes = y.Distinct().ToArray();
                m_nclass = m_classes.Length;
            }
            else
                throw new ArgumentException("Expected a y vector");

            #region Class Label Translation
            //y and ytest must be id's in integer order starting at 1 (see rf.c at around line 172
                /* Count number of cases in each class. */
                /*zeroInt(classFreq, nclass);
                for (n = 0; n < nsample; ++n) classFreq[cl[n] - 1]++;*/
            // therefore, we must keep track of the original classIDs (database) and the classNumbers (1:nclass) for the randomForest
            m_classCorrespondence = new BiDirectionalDictionary<int, int>(m_nclass);
            for (int i = 0; i < m_nclass; i++)
            {
                // key is classID, value is class number
                m_classCorrespondence.Add(m_classes[i], i + 1);
            }
            // replace classes with the translated classes
            m_classes = m_classCorrespondence.getValueWithKey(m_classes);
            // m_classes should now be the first m_nclass integers
            Debug.Assert(m_classes.Sum() == m_nclass * (m_nclass + 1) / 2);
            // replace y with translated classes
            if (y != null)
            {
                y = m_classCorrespondence.getValueWithKey(y);
                // y should now have nothing in it greater than m_nclass
                Debug.Assert(Array.TrueForAll(y, (a) => (a <= m_nclass)));
            }

            // replace ytest with translated classes
            if (ytest != null)
            {
                ytest = m_classCorrespondence.getValueWithKey(ytest);
                // y should now have nothing in it greater than m_nclass
                Debug.Assert(Array.TrueForAll(ytest, (a) => (a <= m_nclass)));
            }
            #endregion

            // Variables from the function call
            m_ntree = ntree;
            if (mtry == -1)
                m_mtry = Convert.ToInt32(Math.Floor(Math.Sqrt(Convert.ToDouble(p))));
            else
                m_mtry = mtry;
            int nodesize = 1;
            // if strata is null and sampsize is null, make a single sampsize
            if (strata == null && sampsize == null)
            {
                sampsize = new int[1];
            }
            // else if strata and sapsize are not null, check that they are equal length
            else if (strata != null && sampsize != null)
            {
                if (strata.Distinct().ToArray().Length != sampsize.Length)
                    throw new ArgumentException("The length of samplesize must equal the number of unique strata");
            }
            // else if strata is not null and sampsize is null
            else if (strata != null && sampsize == null)
            {
                throw new ArgumentException("Strata are defined, but not the samplesize for each.");
            }
            // else if strata are null and sampsize is not null, the strata will be replaced with y and sample sizes checked later
            int? maxnodes = null;
            int nPerm = 1;
            bool bSampleWithReplacement = true;
            bool bLocalImportance = false;
            bool bProximityOnlyOutOfBag = bReturnProximity;
            bool bNormalizeVotes = true;
            bool bDoTrace = false;
            bool bKeepForest = true;
#pragma warning disable 0219
            bool bCorrectBiasForRegression = false; // not used... // TODO:  fix this by using it
#pragma warning restore 0219
            bool bKeepInBagSamples = false;
            
            // initialize some variables
            m_cutoff = null;
            if (sampsize != null && sampsize.Length == 1)
                sampsize[0] = (bSampleWithReplacement) ? n : Convert.ToInt32(Math.Ceiling(0.632 * Convert.ToDouble(n)));

            // other local variables
            int ipi = 0;
            int itestdat, ilabelts;
            int ntest = 0;
            //int error_test = 0;
            MatrixVector<double> prox, impout, impsd, impmat, errtr;
            MatrixVector<int> counttr;
            int nsample, nsum, nt, max_nodes;
            double[] cwt, threshold;
            bool bLabelts, bStratify, bAddClass, bClassRF, bTestdat;
            int[] Options, xdim, outcl;

            // begin function body

            /*
             * addclass <- is.null(y)
             * classRF <- addclass || is.factor(y) // y is a factor
             */
            bAddClass = false;
            bClassRF = true;

            /*
             * if (!classRF && y.Distinct().Count() <= 5)
             *      warning("The response has five or fewer unique values.  Are you sure you want to do regression?")
             * if (classRF && !addclass && length(unique(y)) < 2)
             *      stop("Need at least two classes to do classification.")
             */
            // TODO:  add regression capability 

            if (bClassRF && !bAddClass && y.Distinct().Count() < 2)
                throw new Exception("need at least two classes to perform classification");
            
            /*
             * n <- nrow(x) // already done
             * p <- ncol(x) // already done
             * if (n == 0) stop("data (x) has 0 rows")
             */
            if (n == 0) 
                throw new ArgumentException("Data has 0 rows", "x");

            /*
             * x.row.names <- rownames(x)
             * x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)
             * ## overcome R's lazy evaluation:
             * keep.forest <- keep.forest
             * testdat <- !is.null(xtest)
             */
            bTestdat = !(xtest == null);
            if (bTestdat) // if test data is supplied
            {
                /*
                 * if (ncol(x) != ncol(xtest))
                 *     stop("x and xtest must have same number of columns")
                 * ntest <- nrow(xtest)
                 * xts.row.names <- rownames(xtest)
                 */
                if (p != xtest.GetLength(1))
                    throw new ArgumentException("x and xtest must have same number of columns", "xtest");
                ntest = xtest.GetLength(0);
                m_xtest = new MatrixVector<double>(xtest);
                if (ytest != null && ytest.Length != ntest)
                    throw new ArgumentException("ytest must have the same number of rows as xtest", "ytest");
            }

            /*
             * ## Make sure mtry is in reasonable range.
             * if (mtry < 1 || mtry > p)
             *     warning("invalid mtry: reset to within valid range")
             * mtry <- max(1, min(p, round(mtry)))
             */
            m_mtry = Convert.ToInt32(Math.Max(1, Math.Min(Convert.ToDouble(p), Convert.ToDouble(m_mtry))));
            if (y != null)
            {
                /*
                 * if (length(y) != n) stop("length of response must be the same as predictors")
                 * addclass <- FALSE
                 */
                if (n != y.Length)
                    throw new Exception("length of response must be the same as predictors");
                bAddClass = false;
            }
            else
            {
                /*
                 * if (!addclass) addclass <- TRUE
                 * y <- factor(c(rep(1, n), rep(2, n)))
                 * x <- rbind(x, x)
                 */
                if (!bAddClass)
                    bAddClass = true;
                m_y = new int[n + n];
                for (int i = 0; i < n; i++)
                {
                    m_y[i] = 1;
                    m_y[i+n] = 2;
                }
                MatrixVector<double> temp_x = new MatrixVector<double>(m_x.nrow * 2, m_x.ncol);
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < p; j++)
                    {
                        temp_x[i, j] = m_x[i, j];
                        temp_x[i + n, j] = m_x[i, j];
                    }
                }
                m_x = temp_x;
            }


            /*
             * ## Check for NAs.
             * if (any(is.na(x))) stop("NA not permitted in predictors")
             * if (testdat && any(is.na(xtest))) stop("NA not permitted in xtest")
             * if (any(is.na(y))) stop("NA not permitted in response")
             * if (!is.null(ytest) && any(is.na(ytest))) stop("NA not permitted in ytest")
             */

            /*
             * if (is.data.frame(x)) {
             */
            // Note:  x is always a dataframe in the C# implmentation used for classification
            if (x.Rank == 2)
            {
                /*
                 * mylevels <- function(x) if (is.factor(x)) levels(x) else 0
                 * xlevels <- lapply(x, mylevels)
                 * ncat <- sapply(xlevels, length)
                 * ## Treat ordered factors as numerics.
                 * ncat <- ifelse(sapply(x, is.ordered), 1, ncat)
                 * x <- data.matrix(x)
                 */
                m_ncat = RandomForestUtilities.CreateAndFill(p, 1); // all numeric variables
                m_xlevels = new List<int[]>(p);
                for (int i = 0; i < p; i++)
                {
                    m_xlevels.Add(new int[1] { 0 });
                }
                if (indexCategoricalVariable != null)
                {
                    if (indexCategoricalVariable.Length > p)
                        throw new ArgumentException("Cannot have more categorical indicators than variables", "indexCategoricalVariable");
                    for (int i = 0; i < indexCategoricalVariable.Length; i++)
                    {
                        if (indexCategoricalVariable[i] >= p)
                            throw new ArgumentException(string.Format("The {0} element of indexCategoricalVariable is a {1} (zero based).  This is not allowed.  must be >= 0 and < the number of variables", i, indexCategoricalVariable[i]), "indexCategoricalVariable");
                        int [] temp_test = new int[n];
                        m_xlevels[indexCategoricalVariable[i]] = Array.ConvertAll<double, int>(m_x.GetColumn(indexCategoricalVariable[i]).Distinct().ToArray(), new Converter<double, int>(a => Convert.ToInt32(a)));
                        m_ncat[indexCategoricalVariable[i]] = m_xlevels[indexCategoricalVariable[i]].Length;
                    }
                }

                if (bTestdat)
                {
                    /*
                     * if(!is.data.frame(xtest))
                     *    stop("xtest must be data frame if x is")
                     * xfactor <- which(sapply(xtest, is.factor))
                     * if (length(xfactor) > 0) {
                     *    for (i in xfactor) {
                     *       if (any(! levels(xtest[[i]]) %in% xlevels[[i]]))
                     *          stop("New factor levels in xtest not present in x")
                     *       xtest[[i]] <- factor(xlevels[[i]][match(xtest[[i]], xlevels[[i]])], levels=xlevels[[i]])
                     *    }
                     * }
                     * xtest <- data.matrix(xtest)
                     */
                    // Assume that xtest has categorical variables in the same location as x and that xtest is a data.frame
                }
            }
            else
            {
                /*
                 * ncat <- rep(1, p)
                 * xlevels <- as.list(rep(0, p))
                 */
                m_ncat = RandomForestUtilities.CreateAndFill<int>(p, 1);
                m_xlevels = new List<int[]>(p);
                for (int i = 0; i < p; i++)
                {
                    m_xlevels.Add(new int[1] { 0 });
                }
            }

            /*
             * maxcat <- max(ncat)
             * if (maxcat > 32)
             *    stop("Can not handle categorical predictors with more than 32 categories.")
             */
            m_maxcat = m_ncat.Max();
            if (m_maxcat > 32)
                throw new Exception("Can not handle categorical predictors with more than 32 categories.");

            if (bClassRF)
            {
                /*
                 * nclass <- length(levels(y)) // already set
                 * ## Check for empty classes:
                 * if (any(table(y) == 0)) stop("Can't have empty classes in y.")
                 */
                // Not necessary to perform the above check because the classes are determined from y
                if (ytest != null)
                {
                    /*
                     * if (!is.factor(ytest)) stop("ytest must be a factor")
                     * if (!all(levels(y) == levels(ytest)))
                     *    stop("y and ytest must have the same levels")
                     */
                    int[] ytest_distinct = ytest.Distinct().ToArray();
                    for (int i = 0; i < ytest_distinct.Length; i++)
                        if (!m_classes.Contains(ytest_distinct[i]))
                            throw new Exception("y and ytest must have the same levels");
                }
                if (m_cutoff == null)
                {
                    /*
                     * cutoff <- rep(1 / nclass, nclass)
                     */
                    if (m_nclass != 0)
                        m_cutoff = RandomForestUtilities.CreateAndFill(m_nclass, 1.0 / Convert.ToDouble(m_nclass));
                    else
                        throw new NotFiniteNumberException("division by zero nclass");
                }
                else
                {
                    /*
                     * if (sum(cutoff) > 1 || sum(cutoff) < 0 || !all(cutoff > 0) ||
                     *    length(cutoff) != nclass) {
                     *      stop("Incorrect cutoff specified.")
                     *    }
                     * if (!is.null(names(cutoff))) {
                     *    if (!all(names(cutoff) %in% levels(y))) {
                     *        stop("Wrong name(s) for cutoff")
                     *    }
                     *    cutoff <- cutoff[levels(y)]
                     * }
                     */
                    if (m_cutoff.Sum() > 1.0 || m_cutoff.Sum() < 0.0 || m_cutoff.Min() <= 0 || m_cutoff.Length != m_nclass)
                        throw new Exception("Incorrect cutoff specified.");
                }
                if (classwt != null)
                {
                    /*
                     * if (length(classwt) != nclass)
                     *    stop("length of classwt not equal to number of classes")
                     * ## If classwt has names, match to class labels.
                     * if (!is.null(names(classwt))) {
                     *    if (!all(names(classwt) %in% levels(y))) {
                     *        stop("Wrong name(s) for classwt")
                     *    }
                     *    classwt <- classwt[levels(y)]
                     * }
                     * if (any(classwt <= 0)) stop("classwt must be positive")
                     * ipi <- 1
                     */
                    if (classwt.Length != m_nclass)
                        throw new Exception("Length of classwt not equal to number of classes");
                    if (classwt.Min() <= 0)
                        throw new Exception("classwt must be positive");
                    ipi = 1;
                }
                else
                {
                    /*
                     * classwt <- rep(1, nclass)
                     * ipi <- 0
                     */
                    classwt = RandomForestUtilities.CreateAndFill(m_nclass, 1.0);
                    ipi = 0;
                }
            }
            else
            {
                /*
                 * addclass <- FALSE
                 */
                bAddClass = false;
            }

            /*
             * if (missing(proximity)) proximity <- addclass
             * if (proximity == null)
             *    proximity = addclass;
             */
            // Note:  proximity will never be missing in this implementation, it has a default of true

            if (bReturnProximity) 
            {
                /*
                 * prox <- matrix(0.0, n, n)
                 * proxts <- if (testdat) matrix(0, ntest, ntest + n) else double(1)
                 */
                prox = new MatrixVector<double>(n, n, 0.0);
                if (bTestdat)
                {
                    m_proximity_test = new MatrixVector<double>(ntest, ntest + n, 0.0);
                }
                else
                {
                    m_proximity_test = new MatrixVector<double>(1, 1, 0.0);
                }
            } 
            else 
            {
                /*
                 * prox <- proxts <- double(1)
                 */
                prox = new MatrixVector<double>(1, 1, 1.0);
                m_proximity_test = new MatrixVector<double>(1, 1, 1.0);
            }

            if (bLocalImportance) 
            {
                /*
                 * importance <- TRUE
                 * impmat <- matrix(0, p, n)
                 */
                bReturnImportance = true;
                impmat = new MatrixVector<double>(p,n, 0.0);
            } 
            else
            {
                /*
                 * impmat <- double(1)
                 */
                impmat = new MatrixVector<double>(1,1, 1.0);
            }

            if (bReturnImportance) 
            {
                /*
                 * if (nPerm < 1) nPerm <- as.integer(1) else nPerm <- as.integer(nPerm)
                 */
                if (nPerm < 1)
                {
                    nPerm = 1;
                }
                if (bClassRF) 
                {
                    /*
                     * impout <- matrix(0.0, p, nclass + 2)
                     * impSD <- matrix(0.0, p, nclass + 1)
                     */
                    impout = new MatrixVector<double>(p, m_nclass + 2, 0.0);
                    impsd = new MatrixVector<double>(p, m_nclass + 1, 0.0);
                } 
                else 
                {
                    /*
                     * impout <- matrix(0.0, p, 2)
                     * impSD <- double(p)
                     * names(impSD) <- x.col.names
                     */
                    impout = new MatrixVector<double>(p, 2, 0.0);
                    impsd = new MatrixVector<double>(p, 1, 0.0);
                }
            } 
            else 
            {
                /*
                 * impout <- double(p)
                 * impSD <- double(1)
                 */
                impout = new MatrixVector<double>(p, 1, 0.0);
                impsd = new MatrixVector<double>(1, 1, 0.0);
            }

            /*
             * nsample <- if (addclass) 2 * n else n
             * Stratify <- length(sampsize) > 1
             * if ((!Stratify) && sampsize > nrow(x)) stop("sampsize too large")
             * if (Stratify && (!classRF)) stop("sampsize should be of length one")
             */
            nsample = (bAddClass) ? 2 * n : n;
            bStratify = sampsize.Length > 1;
            if ((!bStratify) && sampsize[0] > n) 
                throw new Exception("sampsize too large");
            if (bStratify && (!bClassRF)) 
                throw new Exception("sampsize should be of length one");
            if (bClassRF) 
            {
                if (bStratify) 
                {
                    /*
                     * if (missing(strata)) strata <- y
                     * if (!is.factor(strata)) strata <- as.factor(strata)
                     * nsum <- sum(sampsize)
                     * if (length(sampsize) > nlevels(strata)) stop("sampsize has too many elements.")  // already done
                     * if (any(sampsize <= 0) || nsum == 0) stop("Bad sampsize specification")
                     * ## If sampsize has names, match to class labels.
                     * if (!is.null(names(sampsize))) {
                     *    sampsize <- sampsize[levels(strata)]
                     * }
                     * if (any(sampsize > table(strata)))
                     *   stop("sampsize can not be larger than class frequency")
                     */
                    if (strata == null)
                        strata = y;
                    else if (strata.Length != n)
                        throw new ArgumentException("strata must have length of observations");
                    nsum = sampsize.Sum();
                    if (sampsize.Length > strata.Distinct().Count())
                        throw new Exception("sampsize has too many elements.");
                    if (sampsize.Min() <= 0 || nsum == 0)
                        throw new Exception("Bad sampsize specification");
                    // TODO:  How do I know that the order of the distinct Strata will be the same order as the sampsize?
                    int[] distinctStrata = strata.Distinct().ToArray();
                    int[] stratasize = new int[distinctStrata.Length];
                    for(int i = 0; i < distinctStrata.Length; i++)
                    {
                        stratasize[i] = strata.Count((a) => (a == distinctStrata[i]));
                        if (sampsize[i] > stratasize[i])
                            throw new ArgumentException("sampsize for a strata or the overall sampsize can not be larger than class frequency");
                        // since this is in the bStratify block, sampsize cannot be length 1
                    }

                } 
                else 
                {
                    /*
                     * nsum <- sampsize
                     */
                    nsum = sampsize[0];
                }
                /*
                 * nrnodes <- 2 * trunc(nsum / nodesize) + 1
                 */
                if (nodesize != 0)
                    m_nrnodes = 2 * (nsum / nodesize) + 1; // does not need truncate because wil do integer division
                else
                    throw new NotFiniteNumberException("Division by zero nodesize");
            } 
            else
            {
                /*
                 * ## For regression trees, need to do this to get maximal trees.
                 * nrnodes <- 2 * trunc(sampsize/max(1, nodesize - 4)) + 1
                 */
                m_nrnodes = 2 * (sampsize[0] / Math.Max(1, nodesize - 4)) + 1;
            }
            if (maxnodes != null) 
            {
                /*
                 * ## convert # of terminal nodes to total # of nodes
                 * maxnodes <- 2 * maxnodes - 1
                 * if (maxnodes > nrnodes) warning("maxnodes exceeds its max value.")
                 * nrnodes <- min(c(nrnodes, max(c(maxnodes, 1))))
                 */
                maxnodes = 2 * maxnodes - 1;
                m_nrnodes = Math.Min(m_nrnodes, Math.Max((int) maxnodes, 1));
            }
            /*
             * ## Compiled code expects variables in rows and observations in columns.
             * x <- t(x)
             * storage.mode(x) <- "double"
             */
            m_x.transpose();
            Debug.Assert(m_x.ncol == n && m_x.nrow == p);
            if (bTestdat) 
            {
                /*
                 * xtest <- t(xtest)
                 * storage.mode(xtest) <- "double"
                 */
                m_xtest.transpose();
                Debug.Assert(m_xtest.nrow == p);
                if (ytest == null) 
                {
                    /*
                     * ytest <- labelts <- 0
                     */
                    ytest = new int[1];
                    ytest[0] = 0;
                    bLabelts = false;
                } 
                else 
                {
                    /*
                     * labelts <- TRUE
                     */
                    bLabelts = true;
                }
            } 
            else 
            {
                /*
                 * xtest <- double(1)
                 * ytest <- double(1)
                 * ntest <- 1
                 * labelts <- FALSE
                 */
                m_xtest = new MatrixVector<double>(1, 1); 
                m_xtest[0, 0] = 0.0;
                ytest = new int[1]; 
                ytest[0] = 0;
                ntest = 1;
                bLabelts = false;
            }
            /*
             * nt <- if (keep.forest) ntree else 1
             */
            nt = (bKeepForest) ? m_ntree : 1;

            if (bClassRF)
            {
                /*
                 * cwt <- classwt
                 * threshold <- cutoff
                 * error.test <- if (labelts) double((nclass+1) * ntree) else double(1)
                 */
                cwt = classwt;
                threshold = m_cutoff;
                m_error_rate_test = (bLabelts) ? new MatrixVector<double>((m_nclass + 1), m_ntree) : new MatrixVector<double>(1,1);
                
                xdim = new int[2];
                xdim[0] = p;
                xdim[1] = n;
                Options = createOptions(bAddClass, bReturnImportance, bLocalImportance,
                    bReturnProximity,bProximityOnlyOutOfBag,bDoTrace,bKeepForest,bSampleWithReplacement,
                    bStratify,bKeepInBagSamples);
                //outclts = as.integer(numeric(ntest)),
                m_predicted_test = new int[ntest];
                //outcl = integer(nsample),
                outcl = new int[nsample];
                //counttr = integer(nclass * nsample),
                counttr = new MatrixVector<int>(m_nclass, nsample); // related to m_votes
                //ndbigtree = integer(ntree),
                m_ndbigtree = new int[m_ntree];
                //nodestatus = integer(nt * nrnodes),
                m_nodestatus = new MatrixVector<int>(m_nrnodes, nt);
                //bestvar = integer(nt * nrnodes),
                m_bestvar = new MatrixVector<int>(m_nrnodes, nt);
                //treemap = integer(nt * 2 * nrnodes),
                m_treemap = new MatrixVector<int>(2, m_nrnodes, nt);
                //nodepred = integer(nt * nrnodes),
                m_nodepred = new MatrixVector<int>(m_nrnodes, nt);
                //xbestsplit = double(nt * nrnodes),
                m_xbestsplit = new MatrixVector<double>(m_nrnodes, nt);
                //testdat = as.integer(testdat),
                itestdat = Convert.ToInt32(bTestdat);
                //countts = double(nclass * ntest),
                m_votes_test = new MatrixVector<double>(m_nclass, ntest);
                //labelts = as.integer(labelts),
                ilabelts = Convert.ToInt32(bLabelts);
                //inbag = if (keep.inbag) matrix(integer(n * ntree), n) else integer(n),
                m_inbag = (bKeepInBagSamples) ? new MatrixVector<int>((n * m_ntree), n) : new MatrixVector<int>(1, n);
                errtr = new MatrixVector<double>(m_nclass + 1, m_ntree);
                if (!bStratify)
                {
                    strata = new int[1];
                    strata[0] = 1;
                }
                // else strata is passed in, else, strata is teh same as y

                #region classRF Function Call
                CsRandomForest.rf.classRF(
                    m_x.values, //x = x,
                    xdim, // xdim = as.integer(c(p, n)),
                    new CsRandomForest.CArray<int>(y), // y = as.integer(y),
                    ref m_nclass, // nclass = as.integer(nclass),
                    new CsRandomForest.CArray<int>(m_ncat), // ncat = as.integer(ncat),
                    ref m_maxcat, // maxcat = as.integer(maxcat),
                    sampsize, // sampsize = as.integer(sampsize),
                    strata, // strata = if (Stratify) as.integer(strata) else integer(1),
                    Options, // Options = as.integer(c(addclass, importance, localImp, proximity, oob.prox, do.trace, keep.forest, replace, Stratify, keep.inbag)),
                    ref m_ntree, // ntree = as.integer(ntree),
                    ref m_mtry, // mtry = as.integer(mtry),
                    ref ipi, // ipi = as.integer(ipi),
                    classwt, // classwt = as.double(cwt),
                    threshold, // cutoff = as.double(threshold),
                    ref nodesize, // nodesize = as.integer(nodesize),
                    outcl, // outcl = integer(nsample),
                    counttr.values, //counttr = integer(nclass * nsample),
                    prox.values,  // prox = prox,
                    impout.values, // impout = impout,
                    impsd.values, // impSD = impSD,
                    impmat.values, // impmat = impmat,
                    ref m_nrnodes, // nrnodes = as.integer(nrnodes),
                    new CsRandomForest.CArray<int>(m_ndbigtree), //ndbigtree = integer(ntree),
                    new CsRandomForest.CArray<int>(m_nodestatus.values), //nodestatus = integer(nt * nrnodes),
                    new CsRandomForest.CArray<int>(m_bestvar.values), //bestvar = integer(nt * nrnodes),
                    new CsRandomForest.CArray<int>(m_treemap.values), //treemap = integer(nt * 2 * nrnodes),
                    new CsRandomForest.CArray<int>(m_nodepred.values), //nodepred = integer(nt * nrnodes),
                    new CsRandomForest.CArray<double>(m_xbestsplit.values), //xbestsplit = double(nt * nrnodes),
                    new CsRandomForest.CArray<double>(errtr.values), //errtr = double((nclass+1) * ntree),
                    ref itestdat, //testdat = as.integer(testdat),
                    m_xtest.values, //xts = as.double(xtest),
                    ytest, //clts = as.integer(ytest),
                    ref ntest, //nts = as.integer(ntest),
                    m_votes_test.values, //countts = double(nclass * ntest),
                    m_predicted_test, //outclts = as.integer(numeric(ntest)),
                    ref ilabelts, //labelts = as.integer(labelts),
                    m_proximity_test.values, // proxts = proxts,
                    new CsRandomForest.CArray<double>(m_error_rate_test.values), //errts = error.test,
                    m_inbag.values, //inbag = if (keep.inbag) matrix(integer(n * ntree), n) else integer(n),
                    m_oRandom
                    );
                #endregion

                #region Class Label Translation
                // change the y and outcl values back to the original classes for the confusion matrix
                if (y != null)
                    y = m_classCorrespondence.getKeyWithValue(y);
                if (m_y != null)
                    m_y = m_classCorrespondence.getKeyWithValue(m_y);
                if (outcl != null)
                    outcl = m_classCorrespondence.getKeyWithValue(outcl);
                if (ytest != null && bTestdat && bLabelts)
                    ytest = m_classCorrespondence.getKeyWithValue(ytest);
                if (m_classes != null)
                    m_classes = m_classCorrespondence.getKeyWithValue(m_classes);
                if (m_predicted != null)
                    m_predicted = m_classCorrespondence.getKeyWithValue(m_predicted);
                if (m_predicted_test != null && bTestdat)
                    m_predicted_test = m_classCorrespondence.getKeyWithValue(m_predicted_test);
                #endregion

                #region output R object
                if (bKeepForest)
                {
                    /*
                     * ## deal with the random forest outputs
                     * max.nodes <- max(rfout$ndbigtree)
                     * treemap <- aperm(array(rfout$treemap, dim = c(2, nrnodes, ntree)),
                     *                c(2, 1, 3))[1:max.nodes, , , drop=FALSE]
                     */
                    max_nodes = m_ndbigtree.Max();
                    int[] dimensionOrder = {1, 0, 2};
                    m_treemap.permute(dimensionOrder);
                    Debug.Assert(m_treemap.nrow == m_nrnodes && m_treemap.ncol == 2 && m_treemap.dim3 == nt);
                }
                if (!bAddClass)
                {
                    /*
                     * ## Turn the predicted class into a factor like y.
                     * out.class <- factor(rfout$outcl, levels=1:nclass, label=levels(y))
                     * names(out.class) <- x.row.names
                     * con <- table(observed = y, predicted = out.class)[levels(y), levels(y)]
                     * con <- cbind(con, class.error = 1 - diag(con)/rowSums(con))
                     */ 

                    m_confusion = ConfusionMatrix.CreateConfusion(y, outcl);
                }
                /*
                 * out.votes <- t(matrix(rfout$counttr, nclass, nsample))[1:n, ]
                 * oob.times <- rowSums(out.votes)
                 * if (norm.votes) out.votes <- t(apply(out.votes, 1, function(x) x/sum(x)))
                 * dimnames(out.votes) <- list(x.row.names, levels(y))
                 * class(out.votes) <- c(class(out.votes), "votes")
                 */
                m_votes = new MatrixVector<double>(Array.ConvertAll<int, double>(counttr.values, new Converter<int, double>(a => Convert.ToDouble(a))).ToArray(), counttr.nrow, counttr.ncol);
                m_votes.transpose();
                Debug.Assert(m_votes.nrow == nsample && m_votes.ncol == m_nclass);
                m_oob_times = m_votes.rowSum();
                if (bNormalizeVotes)
                    m_votes.normalizeByRow();

                if (bTestdat)
                {
                    /*
                     * out.class.ts <- factor(rfout$outclts, levels=1:nclass, label=levels(y)) // already done
                     * names(out.class.ts) <- xts.row.names
                     * out.votes.ts <- t(matrix(rfout$countts, nclass, ntest))
                     * dimnames(out.votes.ts) <- list(xts.row.names, levels(y))
                     * if (norm.votes)
                     *   out.votes.ts <- t(apply(out.votes.ts, 1, function(x) x/sum(x)))
                     * class(out.votes.ts) <- c(class(out.votes.ts), "votes")
                     * if (labelts) 
                     * {
                     *   testcon <- table(observed = ytest, predicted = out.class.ts)[levels(y), levels(y)]
                     *   testcon <- cbind(testcon, class.error = 1 - diag(testcon)/rowSums(testcon))
                     * }
                     */
                    m_votes_test.transpose();
                    Debug.Assert(m_votes_test.nrow == ntest && m_votes_test.ncol == m_nclass);
                    if (bNormalizeVotes)
                    {
                        m_votes_test.normalizeByRow();
                        m_votes_test.transpose();
                    }
                    if (bLabelts)
                    {
                        m_confusion_test = ConfusionMatrix.CreateConfusion(ytest, m_predicted_test);
                    }
                }
                /*
                 * cl <- match.call()
                 * cl[[1]] <- as.name("randomForest")
                 * out <- list(call = cl,
                 * type = if (addclass) "unsupervised" else "classification",
                 * predicted = if (addclass) NULL else out.class,
                 * err.rate = if (addclass) NULL else t(matrix(rfout$errtr,nclass+1, ntree,dimnames=list(c("OOB", levels(y)), NULL))),
                 */
                m_type = (bAddClass) ? "unsupervised" : "classification";
                m_predicted = (bAddClass) ? null : outcl;
                Debug.Assert(errtr.nrow == m_nclass + 1 && errtr.ncol == m_ntree, "Matrix was reshaped");
                m_err_rate = (bAddClass) ? null : errtr.getMatrix(); // Add names

                /*
                 * confusion = if (addclass) NULL else con, // Already done
                 * votes = out.votes,  // Already done and normalized if required
                 * oob.times = oob.times,  // Already Done
                 * classes = levels(y), // Already Done
                 * importance = if (importance) 
                 *     matrix(rfout$impout, p, nclass+2, dimnames = list(x.col.names, c(levels(y), "MeanDecreaseAccuracy", "MeanDecreaseGini")))
                 *     else matrix(rfout$impout, ncol=1, dimnames=list(x.col.names, "MeanDecreaseGini")),
                 * importanceSD = if (importance)
                 *     matrix(rfout$impSD, p, nclass + 1, dimnames = list(x.col.names, c(levels(y), "MeanDecreaseAccuracy")))
                 *     else NULL,
                 */
                if (bReturnImportance)
                {
                    Debug.Assert(impout.nrow == p && impout.ncol == m_nclass + 2, "Matrix has been reshaped");
                    m_importance = impout.getMatrix(); // TODO:  Add names
                    Debug.Assert(impsd.nrow == p && impsd.ncol == m_nclass + 1, "Matrix has been reshaped");
                    m_importanceSD = impsd.getMatrix(); // TODO:  Add names
                }
                else
                {
                    Debug.Assert(impout.ncol == 1, "Matrix has been reshaped");
                    m_importance = impout.getMatrix(); // TODO:  Add names
                    m_importanceSD = null;
                }

                /*
                 * localImportance = if (localImp)  matrix(rfout$impmat, p, n, dimnames = list(x.col.names,x.row.names)) else NULL,
                 */
                if (bLocalImportance)
                {
                    Debug.Assert(impmat.nrow == p && impmat.ncol == n, "Matrix has been reshaped");
                    m_localImportance = impmat.getMatrix(); // TODO:  Add names
                }
                else
                    m_localImportance = null;

                /*
                 * proximity = if (proximity) matrix(rfout$prox, n, n, dimnames = list(x.row.names, x.row.names)) else NULL,
                 */
                if (bReturnProximity)
                {
                    Debug.Assert(prox.nrow == n && prox.ncol == n, "Matrix has been reshaped");
                    m_proximity = prox.getMatrix(); // TODO:  Add names
                }

                /*
                 * ntree = ntree, // already done
                 * mtry = mtry, // already done
                 */
                if (!bKeepForest) // free this memory
                {
                    m_ndbigtree = null;
                    m_nodestatus = null;
                    m_bestvar = null;
                    m_treemap = null;
                    m_nodepred = null;
                    m_xbestsplit = null;
                    m_pid = null;
                    m_cutoff = null;
                    m_ncat = null;
                    m_xlevels = null;
                }
                else
                {
                    /*
                     * forest = if (!keep.forest) NULL else { // already done
                     *    list(ndbigtree = rfout$ndbigtree, // already done
                     *       nodestatus = matrix(rfout$nodestatus, nc = ntree)[1:max.nodes,, drop=FALSE], // already done
                     *       bestvar = matrix(rfout$bestvar, nc = ntree)[1:max.nodes,, drop=FALSE], // already done
                     *       treemap = treemap, // already done
                     *       nodepred = matrix(rfout$nodepred, nc = ntree)[1:max.nodes,, drop=FALSE], // already done
                     *       xbestsplit = matrix(rfout$xbestsplit, nc = ntree)[1:max.nodes,, drop=FALSE], // already done
                     *       pid = rfout$classwt,
                     *       cutoff=cutoff, // already done
                     *       ncat=ncat, // already done
                     *       maxcat = maxcat, // already done
                     *       nrnodes = max.nodes, 
                     *       ntree = ntree, // already done
                     *       nclass = nclass, // already done
                     *       xlevels=xlevels) // already done
                     *   },
                     */
                    m_pid = classwt;
                    m_nrnodes = Convert.ToInt32(maxnodes);
                }

                /*
                 * y = if (addclass) NULL else y,
                 * test = if(!testdat) NULL else list(
                 *    predicted = out.class.ts,  // already done
                 *    err.rate = if (labelts) t(matrix(rfout$errts, nclass+1, ntree, dimnames=list(c("Test", levels(y)), NULL))) else NULL, // already done
                 *      confusion = if (labelts) testcon else NULL,
                 *      votes = out.votes.ts,
                 *      proximity = if(proximity) matrix(rfout$proxts, nrow=ntest,
                 *      dimnames = list(xts.row.names, c(xts.row.names,
                 *      x.row.names))) else NULL),
                 *      inbag = if (keep.inbag) rfout$inbag else NULL)
                 */
                m_y = (bAddClass) ? null : y;
                if (!bTestdat)
                {
                    // free memory
                    m_predicted_test = null;
                    m_error_rate_test = null;
                    m_confusion_test = null;
                    m_votes_test = null;
                    m_proximity_test = null;
                    m_inbag = null;
                }
                else
                {
                    if (!bLabelts)
                        m_error_rate_test = null;
                    if (!bLabelts)
                        m_confusion_test = null;
                    if (!bReturnProximity)
                        m_proximity_test = null;
                    if (!bKeepInBagSamples)
                        m_inbag = null;
                }
                #endregion
            }
            else // !classRF 
            {
                throw new NotImplementedException(); // TODO
                #region regression
                /*
                rfout <- .C("regRF",
                            x,
                            as.double(y),
                            as.integer(c(n, p)),
                            as.integer(sampsize),
                            as.integer(nodesize),
                            as.integer(nrnodes),
                            as.integer(ntree),
                            as.integer(mtry),
                            as.integer(c(importance, localImp, nPerm)),
                            as.integer(ncat),
                            as.integer(maxcat),
                            as.integer(do.trace),
                            as.integer(proximity),
                            as.integer(oob.prox),
                            as.integer(corr.bias),
                            ypred = double(n),
                            impout = impout,
                            impmat = impmat,
                            impSD = impSD,
                            prox = prox,
                            ndbigtree = integer(ntree),
                            nodestatus = matrix(integer(nrnodes * nt), ncol=nt),
                            leftDaughter = matrix(integer(nrnodes * nt), ncol=nt),
                            rightDaughter = matrix(integer(nrnodes * nt), ncol=nt),
                            nodepred = matrix(double(nrnodes * nt), ncol=nt),
                            bestvar = matrix(integer(nrnodes * nt), ncol=nt),
                            xbestsplit = matrix(double(nrnodes * nt), ncol=nt),
                            mse = double(ntree),
                            keep = as.integer(c(keep.forest, keep.inbag)),
                            replace = as.integer(replace),
                            testdat = as.integer(testdat),
                            xts = xtest,
                            ntest = as.integer(ntest),
                            yts = as.double(ytest),
                            labelts = as.integer(labelts),
                            ytestpred = double(ntest),
                            proxts = proxts,
                            msets = double(if (labelts) ntree else 1),
                            coef = double(2),
                            oob.times = integer(n),
                            inbag = if (keep.inbag)
                            matrix(integer(n * ntree), n) else integer(1),
                            DUP=FALSE,
                            PACKAGE="randomForest")[c(16:28, 36:41)]
                ## Format the forest component, if present.
                if (keep.forest) {
                    max.nodes <- max(rfout$ndbigtree)
                    rfout$nodestatus <-
                        rfout$nodestatus[1:max.nodes, , drop=FALSE]
                    rfout$bestvar <-
                        rfout$bestvar[1:max.nodes, , drop=FALSE]
                    rfout$nodepred <-
                        rfout$nodepred[1:max.nodes, , drop=FALSE]
                    rfout$xbestsplit <-
                        rfout$xbestsplit[1:max.nodes, , drop=FALSE]
                    rfout$leftDaughter <-
                        rfout$leftDaughter[1:max.nodes, , drop=FALSE]
                    rfout$rightDaughter <-
                        rfout$rightDaughter[1:max.nodes, , drop=FALSE]
                }
                cl <- match.call()
                cl[[1]] <- as.name("randomForest")
                ## Make sure those obs. that have not been OOB get NA as prediction.
                ypred <- rfout$ypred
                if (any(rfout$oob.times < 1)) {
                    ypred[rfout$oob.times == 0] <- NA
                }
                out <- list(call = cl,
                            type = "regression",
                            predicted = structure(ypred, names=x.row.names),
                            mse = rfout$mse,
                            rsq = 1 - rfout$mse / (var(y) * (n-1) / n),
                            oob.times = rfout$oob.times,
                            importance = if (importance) matrix(rfout$impout, p, 2,
                            dimnames=list(x.col.names,
                                          c("%IncMSE","IncNodePurity"))) else
                                matrix(rfout$impout, ncol=1,
                                       dimnames=list(x.col.names, "IncNodePurity")),
                            importanceSD=if (importance) rfout$impSD else NULL,
                            localImportance = if (localImp)
                            matrix(rfout$impmat, p, n, dimnames=list(x.col.names,
                                                       x.row.names)) else NULL,
                            proximity = if (proximity) matrix(rfout$prox, n, n,
                            dimnames = list(x.row.names, x.row.names)) else NULL,
                            ntree = ntree,
                            mtry = mtry,
                            forest = if (keep.forest)
                            c(rfout[c("ndbigtree", "nodestatus", "leftDaughter",
                                      "rightDaughter", "nodepred", "bestvar",
                                      "xbestsplit")],
                              list(ncat = ncat), list(nrnodes=max.nodes),
                              list(ntree=ntree), list(xlevels=xlevels)) else NULL,
                            coefs = if (corr.bias) rfout$coef else NULL,
                            y = y,
                            test = if(testdat) {
                                list(predicted = structure(rfout$ytestpred,
                                     names=xts.row.names),
                                     mse = if(labelts) rfout$msets else NULL,
                                     rsq = if(labelts) 1 - rfout$msets /
                                                (var(ytest) * (n-1) / n) else NULL,
                                     proximity = if (proximity)
                                     matrix(rfout$proxts / ntree, nrow = ntest,
                                            dimnames = list(xts.row.names,
                                            c(xts.row.names,
                                            x.row.names))) else NULL)
                            } else NULL,
                            inbag = if (keep.inbag)
                            matrix(rfout$inbag, nrow(rfout$inbag),
                                   dimnames=list(x.row.names, NULL)) else NULL)
            } */


                /*****************/

                //double[] x = data.??
                /*double[] x = new double[10]; // TODO
                int[] dimx = new int[2];
                dimx[0] = data.GetLength(0); // TODO:  is this the correct dimension?
                dimx[1] = data.GetLength(1); // TODO:  is this the correct dimension?
                int[] cl = classes;
                int ncl = classes.Length;
                int[] cat = new int[10]; // TODO
                int[] Options = new int[10]; // TODO
                //int ntree = 0; // comes from the function call
                int nvar = 0; // TODO
                int ipi = 0; // TODO
                double? cut = null; // TODO
                int outcl = 0; // TODO
                int[] counttr = new int[10]; // TODO
                double[] prox = new double[10]; // TODO
                double[] imprt = new double[10]; // TODO
                double? impsd = null; // TODO
                double[] impmat = new double[10]; // TODO
                int? nrnodes = null; // TODO
                int? ndbigtree  = null; // TODO
                int? nodestatus = null; // TODO
                int? bestvar = null; // TODO
                int? treemap = null; // TODO
                int? nodeclass = null; // TODO
                double? xbestsplit = null; // TODO
                double? errtr = null; // TODO
                int? testdat = null; // TODO
                double? xts = null; // TODO
                int? clts = null; // TODO
                int? nts = null; // TODO
                double? countts = null; // TODO
                int? outclts = null; // TODO
                int? labelts = null; // TODO
                double? proxts = null; // TODO
                double? errts = null; // TODO
                int? inbag = null; // TODO
                Interface.classRF(x, dimx, cl, ncl, cat, maxcat,
                    sampsize, strata, Options, ntree, nvar,
                    ipi, classwt, cut, nodesize,
                    ref outcl, ref counttr, prox,
                    imprt, impsd, impmat, nrnodes,
                    ndbigtree, nodestatus, bestvar, treemap,
                    nodeclass, xbestsplit, errtr,
                    testdat, xts, clts, nts, countts,
                    outclts, labelts, proxts, errts,
                    inbag);*/
                #endregion regression
            }
        }
    }
}
