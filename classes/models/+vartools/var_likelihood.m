function varargout=var_likelihood(params,obj,choice)%
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
if nargin<3
    choice='new';
end
nout=nargout;
switch choice
    case 'new'
        [varargout{1:nout}]=vartools.var_likelihood_new(params,obj);
    case 'old'
        [varargout{1:nout}]=vartools.var_likelihood_old(params,obj);
    case 'debug'
        tic
        [varargout1{1:nout}]=vartools.var_likelihood_new(params,obj);
        new_time=toc;
        tic
        [varargout2{1:nout}]=vartools.var_likelihood_old(params,obj);
        old_time=toc;
        varargout=varargout1;
        if varargout1{3}||varargout2{3}||...
                abs(varargout1{1}-varargout2{1})>1e-8
            keyboard
        else
            fprintf('timing new: %0.8f, timing old: %0.8f\n',new_time,old_time)
        end
    otherwise
        error(['unknown option "',choice,'"'])
end

end