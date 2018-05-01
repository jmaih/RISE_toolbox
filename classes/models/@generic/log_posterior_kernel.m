function [log_post,log_lik,log_prior,Incr,retcode,obj]=log_posterior_kernel(obj,param)
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


nobj=numel(obj);

if nargin<2
    
    param=[];
    
end

if nobj==0
    
    if nargout>1
        
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
        
    end
    
    log_post=cell(0,4);
    
    return
    
elseif nobj>1
    
    nout=nargout;
    
    vout=cell(1,nout);
    
    bigvout=vout;
    
    for iobj=1:nobj
        
        [vout{1:nout}]=log_posterior_kernel(obj(iobj),param);
        
        for io=1:nout
            
            bigvout{io}=[bigvout{io},vout{io}];
            
        end
        
    end
    
    if nout
        
        log_post=bigvout{1};
        
        if nout>1
            
            log_lik=bigvout{2};
            
            if nout>2
                
                log_prior=bigvout{3};
                
                if nout>3
                    
                    Incr=bigvout{4};
                    
                    if nout>4
                        
                        retcode=bigvout{5};
                        
                        if nout>5
                            
                            obj=bigvout{6};
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        return
        
    end
    
end
    
    % estim_hyperparams=obj.estim_hyperparams;
    log_prior=nan(1,2);
    
    log_lik=-obj.options.estim_penalty;
    
    log_post=log_lik;
    
    Incr=[];
    
    likelihood_func=obj.routines.likelihood;
    
    [log_prior(1),retcode]=log_prior_density(obj,param);
    
    if ~retcode
        
        [log_lik,Incr,retcode,obj]=likelihood_func(param,obj);
        % under dsge-var, log_lik can be a vector such that the first element
        % is the var-dsge and the second element is the dsge. retcode also can
        % be a two-element vector
        if ~retcode(1) % pick the right retcode
            
            [log_prior(2),retcode(1)]=log_prior_density(obj,'endogenous');
            
            if ~retcode(1)
                
                log_post=log_lik+sum(log_prior);
                
                if log_post<=-obj.options.estim_penalty
                    
                    retcode(1)=306; % unlikely parameter vector
                    
                end
                
            end
            
        end
        
    end

if obj.options.debug
    
    disp(['log_lik ',num2str(log_lik)])
    
    disp(['log_prior ',num2str(log_prior(1))])
    
    disp(['log_endo_prior ',num2str(log_prior(2))])
    
    disp(['log_post ',num2str(log_post)])
    
    utils.error.decipher(retcode)
    
    if obj.estimation_under_way
        
        error([mfilename,':: no debugging of this function should occur during estimation, at least I do not see why'])
        
    else
        
        keyboard
        
    end
    
end

end

