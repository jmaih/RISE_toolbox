function data=bsxfun(db,fun,b)
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


data=bsxfun(fun,db.data,b);

data=ts(db.start,data);

end