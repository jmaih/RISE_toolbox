function [P,retcode,good]=doubling_solve(A,B,C,options)
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
% More About
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
	P=struct('lyapunov_double_stationary',false,...
        'lyapunov_double_stationary_threshold',5000);
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

if isempty(options)
    options=doubling_solve();
end

P0=C;
symmetric=isequal(A,B');
Gl=A;
if ~symmetric
    Gr=B;
end

[P,~,retcode]=fix_point_iterator(@iterator,P0,options);

if retcode
    if options.lyapunov_double_stationary
        do_stationary_variables_only()
    end
else
    good=true(size(P,1),1);
end

    function do_stationary_variables_only()
        % select the stationary guys and proceed: Looks like we can concentrate
        % on the variances, so take the diagonal. Seems to save time
        %----------------------------------------------------------------------
        good=diag(P)<options.lyapunov_double_stationary_threshold;
        if any(good)
            P20=P(good,good);
            Gl=Gl(good,good);
            if ~symmetric
                Gr=Gr(good,good);
            end
            [P2,~,retcode]=fix_point_iterator(@iterator,P20,options);
            P(:)=nan;
            P(good,good)=P2;
        end
        if retcode
            retcode=280+retcode;
        end
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