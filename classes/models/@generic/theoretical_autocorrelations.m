function [A,retcode]=theoretical_autocorrelations(obj,varargin)%,resolve_flag
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

if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        A=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
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

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

A=theoretical_autocovariances(dsge.empty);

A=disp_defaults(A);

d={
    'autocorr_ar',A.autocov_ar,@(x)num_fin_int(x),...
    'autocorr_ar must be a finite and positive integer'
    };
    
end
