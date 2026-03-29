%--- help for ts/convert ---
%
%  CONVERT Converts data in different frequencies.
% 
%  Usage:
% 
%    newdb = convert(db, newFrequency)
% 
%    newdb = convert(db, newFrequency, squash)
% 
%  Inputs
% 
%  - db: Time series data structure (ts class) with the following fields:
% 
%    - db.dates: Dates of the time series data
% 
%    - db.data: Data values corresponding to the dates
% 
%  - newFrequency: Desired frequency for the converted time series. Supported frequencies are:
% 
%    - 'D': Daily
%    - 'W': Weekly
%    - 'M': Monthly
%    - 'Q': Quarterly
%    - 'H': Half-yearly
%    - 'Y': Yearly
% 
%  - squash (optional): Aggregation function used to combine data points
%    within each time period. It can be either a function handle or a
%    character array representing a function name. If not provided, the
%    default aggregation function is 'mean'. NOTE : The squash function will
%    be applied ONLY when going from higher to lower frequency
%    (aggregation). Otherwise, the higher frequency is written as is in the
%    dates corresponding to the observations, and nan elsewhere.
% 
%  Outputs:
% 
%  - newdb: Converted time series data structure (ts class) with the following fields:
%    - newdb.dates: Dates of the converted time series
%    - newdb.data: Data values corresponding to the dates
% 
%    Example:
%    dates = rd(2022, 1, 1);
%    data = randn(365, 1);
%    db = ts(dates, data);
%    newdb = convert(db, 'M', @sum);
% 
%    See also MEAN, SUM.
% 
%  Notes:
% 
%  - The function converts the data in the time series db into the specified
%    newFrequency by aggregating the data using the provided squash function.
%  - If the newFrequency is the same as the original frequency of the time
%    series, the function returns the input db unchanged.
%  - The function supports different aggregation functions specified by the
%    squash argument. The default aggregation function is 'mean'.
%  - The input db must be a time series represented by the ts class, with
%    dates and data values provided.
%  - The function works with daily, weekly, monthly, quarterly, half-yearly
%    and yearly frequencies.
%  - For frequencies higher than daily (e.g., weekly, monthly), the function
%    aggregates the data within each time period based on the specified squash
%    function.
%
%    Other uses of convert
%
%       matlab.net.http.HeaderField/convert
%       matlab.net.http.io.StringConsumer/convert
%       rise_dates.dates/convert
%