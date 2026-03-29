% ISFINITE True for finite elements.
%    ISFINITE(X) returns an array that contains 1's where
%    the elements of X are finite and 0's where they are not.
%    For example, ISFINITE([pi NaN Inf -Inf]) is [1 0 0 0].
% 
%    For any X, exactly one of ISFINITE(X), ISINF(X), or ISNAN(X)
%    is 1 for each element.
% 
%    See also ISNAN, ISINF.
%
%    Documentation for isfinite
%       doc isfinite
%
%    Other uses of isfinite
%
%       calendarDuration/isfinite    duration/isfinite    symbolic/isfinite
%       codistributed/isfinite       gpuArray/isfinite    ts/isfinite
%       datetime/isfinite            sym/isfinite
%