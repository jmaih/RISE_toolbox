%--- index.m not found. Showing help for rsindex instead. ---
%
% rsindex - Relative Strength Index (RSI)
%    This MATLAB function calculates the Relative Strength Index (RSI) from
%    the series of closing stock prices.
%
%    Syntax
%      index = rsindex(Data)
%      index = rsindex(___,Name,Value)
%
%    Input Arguments
%      Data - Data with closing prices
%        matrix | table | timetable
%
%    Name-Value Arguments
%      WindowSize - Moving window size for relative strength index
%        14 (default) | positive integer
%
%    Output Arguments
%      index - Relative strength index
%        matrix | table | timetable
%
%    Examples
%      openExample('finance/CalculateRelativeStrengthIndexForDataSeriesExample')
%
%    See also timetable, table, negvolidx, posvolidx
%
%    Introduced in Financial Toolbox before R2006a
%    Documentation for rsindex
%       doc rsindex
%
%