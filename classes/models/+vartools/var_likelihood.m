function varargout=var_likelihood(params,obj)%
% H1 line
%
% Syntax
% -------
% ::
%   [LogLik,Incr,retcode,obj]=var_likelihood(params,obj,choice)
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
% See also: store_probabilities save_filters
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