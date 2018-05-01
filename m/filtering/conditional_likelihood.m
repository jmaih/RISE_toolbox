function [loglik,v_iF_v,lik]=conditional_likelihood(v,iF,dF,pp)
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

v_iF_v=v'*iF*v;

loglik=-.5*(pp*log(2*pi)+log(dF)+v_iF_v);

if nargout>2

    lik=exp(loglik);
	
end

end
