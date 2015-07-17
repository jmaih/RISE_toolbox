function db=cumsum(db,dim)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if nargin<2
    dim=1;
end

db=ts(db.date_numbers,cumsum(db.data,dim),db.varnames);
end