function db=cumprod(db,dim)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

if nargin<2
    dim=1;
end

db=ts(db.date_numbers,cumprod(db.data,dim),db.varnames);
end