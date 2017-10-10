function varargout=var_likelihood(params,obj)
% var_likelihood -- computes the likelihood of a VAR model.
%
% More About
% ------------
%
% - This function has he same inputs and outpus as
% VARTOOLS/VAR_LIKELIHOOD_STRETCH and VARTOOLS/VAR_LIKELIHOOD_DIRECT . It
% calls one of the function or both, depending on option
%   - **estim_likelihood_exp_mode** ['direct'|'debug'|{'stretch'}]:
%
% Examples
% ---------
%
% See also: VARTOOLS/VAR_LIKELIHOOD_STRETCH, VARTOOLS/VAR_LIKELIHOOD_DIRECT

nout=nargout;

if nargin==0
    
    varargout=cell(1,nout);
    
    varargout{1}=struct('estim_likelihood_exp_mode','stretch');
    
else
    
    switch obj.options.estim_likelihood_exp_mode
        
        case 'stretch'
            
            [varargout{1:nout}]=vartools.var_likelihood_stretch(params,obj);
        
        case 'direct'
            
            [varargout{1:nout}]=vartools.var_likelihood_direct(params,obj);
        
        case 'debug'
            
            tic
            
            [varargout1{1:nout}]=vartools.var_likelihood_stretch(params,obj);
            
            new_time=toc;
            
            tic
            
            [varargout2{1:nout}]=vartools.var_likelihood_direct(params,obj);
            
            old_time=toc;
            
            varargout=varargout1;
            
            if varargout1{3}||varargout2{3}||...
                    abs(varargout1{1}-varargout2{1})>1e-8
                
                keyboard
            
            else
                
                fprintf('timing stretch algo: %0.8f, timing direct algo: %0.8f\n',new_time,old_time)
            
            end
            
        otherwise
            
            error(['unknown option "',parser.any2str(obj.options.estim_likelihood_exp_mode),...
                '". Valid options are ''stretch'', ''direct'' and ''debug'''])
    
    end
    
end

end