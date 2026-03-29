%  `ts` - Time Series Object
% 
%  **Description:**
%    The `ts` class represents a time series object in MATLAB, designed to
%    store and manipulate time series data. It inherits from the `gogetter`
%    class.
% 
%  **Properties:**
%    - `varnames` (:obj:`cellstr`): Names of the variables in the database.
%    - `data` (:obj:`numeric`): Time series data with dimensions nobs x nvars x npages.
%    - `description` (:obj:`cellstr` or :obj:`{''}`): Comments on each variable.
% 
%  **Dependent Properties:**
%    - `NumberOfObservations` (:obj:`numeric`): Number of observations in the time series.
%    - `NumberOfPages` (:obj:`numeric`): Number of pages (third dimension) of the time series.
%    - `NumberOfVariables` (:obj:`numeric`): Number of variables in the time series.
%    - `finish` (:obj:`char` or :obj:`serial date`): End time of the time series.
%    - `frequency` (:obj:`char`): Frequency of the time series.
%    - `start` (:obj:`char` or :obj:`serial date`): Start time of the time series.
%    - `date_numbers` (:obj:`numeric`): Numeric representation of dates.
%    - `start_date_number` (:obj:`numeric`): Numeric representation of the start date.
%    - `dates` (:obj:`datetime`): Date objects representing the time series dates.
%    - Hidden Property: `firstDateObj` (:obj:`serial date`): Hidden property representing the first date.
% 
%  **Methods:**
%    - **Constructor:**
%        - `ts`: Construct a time series object.
% 
%    - **Setter Method:**
%        - `set.start`: Sets the start date of the ts object.
% 
%    - **Visualization Methods:**
%        - `allmean`: Compute mean for all data.
%        - `apply`: Apply a function to the data.
%        - `bsxfun`: Binary Singleton Expansion Function.
%        - `chebyshev_box`: Create a Chebyshev box plot.
%        - `describe`: Display summary statistics.
%        - `disp`: Display time series information.
%        - `expanding`: Expanding window operation.
%        - `fanchart`: Create a fan chart.
%        - `hpfilter`: Apply the Hodrick-Prescott filter.
%        - `interpolate`: Interpolate missing values.
%        - `ma_filter`: Moving average filter.
%        - `moments`: Compute moments of the data.
%        - `npdecomp`: Non-parametric decomposition.
%        - `pdecomp`: Parametric decomposition.
%        - `regress`: Perform regression analysis.
%        - `rolling`: Rolling window operation.
%        - `spectrum`: Compute the spectrum of the data.
%        
%    - **Statistics Methods:**
%        - `corr`: Compute correlation coefficient.
%        - `corrcoef`: Compute correlation coefficient matrix.
%        - `cumprod`: Cumulative product of the data.
%        - `cumsum`: Cumulative sum of the data.
%        - `jbtest`: Perform the Jarque-Bera test for normality.
%        - `kurtosis`: Compute kurtosis of the data.
%        - `mean`: Compute mean of the data.
%        - `median`: Compute median of the data.
%        - `mode`: Compute mode of the data.
%        - `prctile`: Compute percentiles of the data.
%        - `rmse`: Compute root mean squared error.
%        - `skewness`: Compute skewness of the data.
%        - `std`: Compute standard deviation of the data.
%        - `sum`: Compute sum of the data.
%        - `var`: Compute variance of the data.
%        
%    - **Graphing Methods:**
%        - `bar`: Create a bar graph.
%        - `barh`: Create a horizontal bar graph.
%        - `boxplot`: Create a box plot.
%        - `hist`: Create a histogram.
%        - `plot`: Create a line plot.
%        - `plotyy`: Create a plot with y-axes on both sides.
%        - `plot_real_time`: Create a real-time plot.
%        
%    - **Calculus Methods:**
%        - `abs`: Compute absolute values.
%        - `acos`: Compute inverse cosine.
%        - `acosh`: Compute inverse hyperbolic cosine.
%        - `acot`: Compute inverse cotangent.
%        - `acoth`: Compute inverse hyperbolic cotangent.
%        - `asin`: Compute inverse sine.
%        - `asinh`: Compute inverse hyperbolic sine.
%        - `atan`: Compute inverse tangent.
%        - `atanh`: Compute inverse hyperbolic tangent.
%        - `cos`: Compute cosine.
%        - `cosh`: Compute hyperbolic cosine.
%        - `cot`: Compute cotangent.
%        - `coth`: Compute hyperbolic cotangent.
%        - `cov`: Compute covariance matrix.
%        - `eq`: Element-wise equality comparison.
%        - `exp`: Compute exponential.
%        - `ge`: Element-wise greater than or equal to comparison.
%        - `gt`: Element-wise greater than comparison.
%        - `le`: Element-wise less than or equal to comparison.
%        - `log`: Compute natural logarithm.
%        - `lt`: Element-wise less than comparison.
%        - `max`: Element-wise maximum.
%        - `min`: Element-wise minimum.
%        - `minus`: Element-wise subtraction.
%        - `mpower`: Matrix power.
%        - `mrdivide`: Matrix right division.
%        - `mtimes`: Matrix multiplication.
%        - `ne`: Element-wise inequality comparison.
%        - `plus`: Element-wise addition.
%        - `power`: Element-wise power.
%        - `rdivide`: Element-wise right division.
%        - `sin`: Compute sine.
%        - `sinh`: Compute hyperbolic sine.
%        - `times`: Element-wise multiplication.
%        - `uminus`: Unary minus.
%        
%    - **Lookaround Methods:**
%        - `head`: Display the first few observations.
%        - `subsasgn`: Subscripted assignment.
%        - `subsref`: Subscripted reference.
%        - `tail`: Display the last few observations.
%        - `values`: Extract data values.
%        - `double`: Convert to
%
%    Documentation for ts
%       doc ts
%
%