function [M1,M2,RM2i,q,U,S,retcode]=null_and_column_spaces(R,verbose)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<2
    verbose=false;
end
[q,cols]=size(R);

DegreesOfFreedom=cols-q;
if DegreesOfFreedom<0
    error([mfilename,':: More Restrictions than shocks, no solution'])
else
	if verbose
		if DegreesOfFreedom==0
		    disp('There is a unique solution')
		else
		    disp('There are many solutions')
		end
	end
end

% Avoid complaints
%-----------------
R=full(R);

[U,S,V]=svd(R);

% this can be treacherous if q==1
if q>1
    dS=diag(S);
    failed= ~isempty(dS) && dS(end)<sqrt(eps);
else
    failed=~any(R);
end
if failed
    % Avoid calling SVD twice by computing the rank first
    %-----------------------------------------------------
    disp([mfilename,':: Restrictions are not independent or they are not ',...
        'affected by the shocks selected. You may want to check your ',...
        'restrictions or the shocks that are meant to explain those ',...
        'restrictions'])
    M1=[];M2=[];RM2i=[];U=[];S=[];
    retcode=702;
    return
end
retcode=0;

M2=V(:,1:q); % another basis which is different but gives the same end results is M22=null(M1') 

M1=V(:,q+1:cols); % == null(R)

RM2i=(R*M2)\eye(q); 

%M2RM2=M2/(R*M2); % which gives the same results as R'/(R*R'), the WZ result
% but is different from R\eye(q), the solution that Ida would prefer...
end