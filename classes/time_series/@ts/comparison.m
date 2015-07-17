function flag=comparison(func,db1,db2)
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

func=str2func(['@(x,y)',func,'(x,y)']);
if isa(db1,'ts')
    db1=db1.data;
end
if isa(db2,'ts')
    db2=db2.data;
end
flag=func(db1,db2);
end