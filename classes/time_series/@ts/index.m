%--- help for rsindex ---
%
% RSINDEX Relative Strength Index (RSI).
% 
%  Syntax: 
% 
%    index = rsindex(Data)
%    index = rsindex(Data,WindowSize)
% 
%  Description:
% 
%    RSINDEX calculates the Relative Strength Index (RSI) from the series of
%    closing stock prices. By default, RSI values are based on a 14-period window.
% 
%  Input Argument:
% 
%    Data    - A vector, table, or timetable. For matrix input, Data is an 
%              M-by-1 vector of closing prices. Timetables and tables with M 
%              rows contain variables named 'Close' (case insensitive).
% 
%  Optional Input Argument:
% 
%    WindowSize          - Positive integer scalar indicating the moving window
%                          size for relative strength index. The default is 14.
% 
%  Output Argument:
% 
%    index   - Relative StrengthIndex with the same number of rows (M) and 
%              type as the input data.
% 
%  Note: 
%    The RS factor is calculated by dividing the average of the gains by the
%    average of the losses within a specified period.
% 
%          RS = (average gains) / (average losses)
% 
%    Also, the first value of RSI, RSI(1), is a NaN in order to preserve the
%    dimensions of CLOSEP.
% 
%  Example:   
%               load SimulatedStock.mat
%               index = rsindex(TMW)
%               index = rsindex(TMW,14)
% 
%    See also NEGVOLIDX, POSVOLIDX.
% 
%    Reference: Murphy, John J., Technical Analysis of the Futures Market,
%               New York Institute of Finance, 1986, pp. 295-302
%
%    Reference page in Doc Center
%       doc rsindex
%
%    Other functions named rsindex
%
%       fints/rsindex
%