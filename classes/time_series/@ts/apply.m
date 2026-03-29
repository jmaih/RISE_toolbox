%--- help for ts/apply ---
%
%  APPLY - Applies a unary function to each element of a time series database
% 
%  Syntax:
%    db = apply(db, fhandle)
% 
%  Inputs:
%    db - Time series database (ts object)
%    fhandle - Unary function handle that accepts a scalar input and returns a scalar output
% 
%  Outputs:
%    db - Resulting time series database after applying the unary function to each element
% 
%  Description:
%    The 'apply' function applies a unary function to each element of a time series database.
%    It performs the operation element-wise and returns a new time series database with the transformed values.
% 
%  Examples:
%    % Create a time series database
%    dates = rd('2000-01-01');
%    data = rand(5, 3);
%    varnames = {'Var1', 'Var2', 'Var3'};
%    db = ts(dates, data, varnames);
% 
%    % Define a unary function
%    fhandle = @(x) log(x);
% 
%    % Apply the function to each element of the time series database
%    transformed_db = apply(db, fhandle);
% 
%    % Display the resulting time series database
%    disp(transformed_db);
% 
%  See also:
%    ts, unary_operation
%