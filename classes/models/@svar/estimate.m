function varargout=estimate(obj,varargin)
% ESTIMATE - estimates the parameters of a svar or rfvar model
%
% More About
% ------------
%
% - ESTIMATE is the same as RISE_GENERIC/ESTIMATE except for the fact that
% it adds the following additional options
%   - **estim_analytical_post_mode** [false|{true}]: Computes the posterior
%   mode and does the simulation of the posterior analytically.
%   - **estim_likelihood_exp_mode** ['direct'|'debug'|{'stretch'}]: When
%   the likelihood is very low, the likelihood function may set it to
%   infinity very quickly because of the exponentiation. The **stretch**
%   mode transforms the exponential and thereby extends the search region.
%   The **direct** mode computes the likelihood directly but inadvertently
%   truncates the search space. The **debug** mode computes both
%   alternatives and pauses when the direct mode returns an infinite
%   likelihood.
%
% Examples
% ---------
%
% See also: RISE_GENERIC/ESTIMATE

nout=nargout;
varargout=cell(1,nout);
if isempty(obj)
    est_opts=estimate@rise_generic(obj);
    % option for computing the posterior mode analytically if possible
    est_opts.estim_analytical_post_mode=true;
    est_opts=utils.miscellaneous.mergestructures(est_opts,...
        vartools.var_likelihood());
    if nout
        if nout~=1
            error('number of output arguments can only be 1 in this case')
        end
        varargout=cell(1,1);
        varargout{1}=est_opts;
    else
        disp(est_opts)
    end
else
    [varargout{1:nout}]=estimate@rise_generic(obj,varargin{:});
end
end
