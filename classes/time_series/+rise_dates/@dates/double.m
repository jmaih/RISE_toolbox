% DOUBLE Convert to double precision.
%    DOUBLE(X) returns the double precision value for X.
%    If X is already a double precision array, DOUBLE has
%    no effect.
% 
%    DOUBLE is called for the expressions in FOR, IF, and WHILE loops
%    if the expression isn't already double precision.  DOUBLE should
%    be overloaded for all objects where it makes sense to convert it
%    into a double precision value.
% 
%    See also SINGLE, DATATYPES, ISFLOAT, ISNUMERIC.
%
%    Documentation for double
%       doc double
%
%    Other uses of double
%
%       calendarDuration/double    opaque/double
%       categorical/double         rise_dates.dates/double
%       codistributed/double       string/double
%       dataset/double             sym/double
%       datetime/double            tabular/double
%       duration/double            ts/double
%       gpuArray/double
%