%--- embed.m not found. Showing help for tsne instead. ---
%
% tsne - t-Distributed Stochastic Neighbor Embedding
%    This MATLAB function returns a matrix of two-dimensional embeddings of
%    the high-dimensional rows of X.
%
%    Syntax
%      Y = tsne(X)
%      Y = tsne(X,Name,Value)
%      [Y,loss] = tsne(___)
%
%    Input Arguments
%      X - Data points
%        n-by-m matrix
%
%    Name-Value Arguments
%      Algorithm - tsne algorithm
%        'barneshut' (default) | 'exact'
%      CacheSize - Size of Gram matrix in megabytes
%        1e3 (default) | positive scalar | "maximal"
%      Distance - Distance metric
%        'euclidean' (default) | 'seuclidean' | 'fasteuclidean' |
%        'fastseuclidean' | 'cityblock' | 'chebychev' | 'minkowski' |
%        'mahalanobis' | 'cosine' | 'correlation' | 'spearman' | 'hamming' |
%        'jaccard' | function handle
%      Exaggeration - Size of natural clusters in data
%        4 (default) | scalar value 1 or greater
%      NumDimensions - Dimension of the output Y
%        2 (default) | positive integer
%      NumPCAComponents - PCA dimension reduction
%        0 (default) | nonnegative integer
%      Perplexity - Effective number of local neighbors of each point
%        30 (default) | positive scalar
%      Standardize - Flag to normalize input data
%        false (default) | true
%      InitialY - Initial embedded points
%        1e-4*randn(N,NumDimensions) (default) |
%        n-by-NumDimensions real matrix
%      LearnRate - Learning rate for optimization process
%        500 (default) | positive scalar
%      NumPrint - Iterative display frequency
%        20 (default) | positive integer
%      Options - Optimization options
%        structure containing the fields 'MaxIter', 'OutputFcn', and 'TolFun'
%      Theta - Barnes-Hut tradeoff parameter
%        0.5 (default) | scalar from 0 through 1
%      Verbose - Iterative display
%        0 (default) | 1 | 2
%
%    Output Arguments
%      Y - Embedded points
%        n-by-NumDimensions matrix
%      loss - Kullback-Leibler divergence
%        nonnegative scalar
%
%    Examples
%      openExample('stats/VisualizeFisherIrisDataExample')
%      openExample('stats/VisualizeFisherIrisDataWithOptionsExample')
%      openExample('stats/PlotResultsWithNaNInputDataExample')
%      openExample('stats/ObtainTSNELossExample')
%
%    See also pca, pdist, knnsearch, statset, gscatter
%
%    Introduced in Statistics and Machine Learning Toolbox in R2017a