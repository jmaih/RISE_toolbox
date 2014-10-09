function [P,retcode]=doubling_solve(A,B,C,options)
% doubling_solve solves the linear equation X=A*X*B+C
%
% Syntax
% -------
% ::
%   [P,retcode]=doubling_solve(A,B,C)
%   [P,retcode]=doubling_solve(A,B,C,options)
%
% Inputs
% -------
% - A :
% - B :
% - C :
% - options :
% 
% Outputs
% --------
% - P :
% - retcode :
%
% Description
% ------------
% 
% Examples
% ---------
%
% See also: 

if nargin==0
	if nargout>1
		error([mfilename,':: number of output arguments cannot exceed 1 if there are no inputs'])
	end
	P=fix_point_iterator();
	return
end

if nargin<4
    options=[];
    if nargin<3
        error([mfilename,':: at least 3 arguments should be provided'])
    end
end
if isempty(B),B=A';end
if isempty(A),A=B';end

P0=C;
symmetric=isequal(A,B');
Gl=A;
if ~symmetric
    Gr=B;
end

[P,~,retcode]=fix_point_iterator(@iterator,P0,options);

if retcode
    retcode=280+retcode;
end

    function [P,F0]=iterator(P0)
        if symmetric
            P=P0+Gl*P0*Gl';
        else
            P=P0+Gl*P0*Gr;
            Gr=Gr*Gr;
        end
        Gl=Gl*Gl;
        F0=P-P0;
   end
end