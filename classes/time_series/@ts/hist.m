% HIST  Histogram.
%    HIST is not recommended. Use HISTOGRAM or HISTCOUNTS instead.
%  
%    N = HIST(Y) bins the elements of Y into 10 equally spaced containers
%    and returns the number of elements in each container.  If Y is a
%    matrix, HIST works down the columns.
% 
%    N = HIST(C) returns the category counts for the categorical array C.
%    For a categorical matrix, HIST works down the columns of Y and returns
%    a matrix of counts with one column for each column of Y and one row for
%    each category.
% 
%    N = HIST(Y,M), where M is a scalar, uses M bins.
% 
%    N = HIST(Y,X), where X is a vector, returns the distribution of Y among
%    bins with centers specified by X. The first bin includes data between
%    -inf and the first center and the last bin includes data between the
%    last bin and inf. Note: Use HISTC if it is more natural to specify bin
%    edges instead.
% 
%    N = HIST(C,CATS) returns counts for the categories specified by CATS.
%    CATS is a categorical array, string array, or a cell array of character
%    vectors.
% 
%    [N,X] = HIST(...) also returns in X the position of the bin centers for
%    numeric data, or the categories corresponding to N for categorical
%    data.
% 
%    HIST(...) without output arguments produces a histogram bar plot of the
%    results. For numeric data, the bar edges on the first and last bins may
%    extend to cover the min and max of the data unless a matrix of data is
%    supplied.
% 
%    HIST(AX,...) plots into AX instead of GCA.
% 
%    Class support for inputs Y: double, single, categorical
%                             X: double, single, categorical, string array
% 
%    See also HISTOGRAM, HISTCOUNTS, MODE, COUNTCATS, CATEGORIES.
%
%    Documentation for hist
%       doc hist
%
%    Other uses of hist
%
%       ts/hist
%