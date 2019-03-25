%--- help for ts/transform ---
%
%  Apply different transformations following the definitions of Haver
% 
%  ::
% 
%     db = transform(db, type);
% 
%  Args:
%     db (ts object): time series object
%     type: the type of the transformation. Available transformations are
% 
%        - 1 or 'level': Untransformed data - level (default)
%        - 2 or 'pct_ch_ar': Percentage change at compound annual rate (% chg)
%        - 3 or 'pct_ch': Period to period percentage change (% change)
%        - 4 or 'yoy_pct_ch': Year over year percentage change (yr/yr % chg)
%        - 5 or 'diff': Period to period difference change
%        - 6 or 'yoy_diff': Year over Year difference change
%        - 7 or 'log_ch_ar': %Log Change - compound annual rate (ann l-chg)
%        - 8 or 'log_ch': %Log Change - period to period (log-change)
%        - 9 or 'yoy_log_ch': Year to year log change (yr-yr l-chg)
% 
%  Returns:
%     :
% 
%     - **db** (ts object): Transformed time series
% 
%
%    Other functions named transform
%
%       abstvar/transform    sym/transform
%