function [A,retcode]=theoretical_autocorrelations(obj,varargin)%,resolve_flag
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

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    A=theoretical_autocovariances(obj);
    A=struct('autocorr_ar',A.autocov_ar);
    return
end

nobj=numel(obj);
if nobj>1
    A=cell(1,nobj);
    retcode=cell(1,nobj);
    for iobj=1:nobj
        [A{iobj},retcode{iobj}]=theoretical_autocorrelations(obj(iobj),varargin{:});
    end
    return
end
obj=set(obj,varargin{:});
autocorr_ar=obj.options.autocorr_ar;
[A,retcode]=theoretical_autocovariances(obj,'autocov_ar',autocorr_ar);%,resolve_flag

if ~retcode
    stdevmat=sqrt(diag(A(:,:,1)));
    bad=stdevmat<=0;
    if any(bad)
        disp(obj.endogenous.name(bad))
        warning([mfilename,':: have zero or negative theoretical standard deviation'])
    end
    stdevmat=stdevmat*stdevmat';
    for ii=1:autocorr_ar+1
        A(:,:,ii)=A(:,:,ii)./stdevmat;
    end
end
