function [smpl,fsmpl,accept_rate,start] = mh_sampling_engine(start,nsamples,varargin)
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

% MH_SAMPLING_ENGINE Generate Markov chain using Metropolis-Hasting algorithm
%   SMPL = MH_SAMPLING_ENGINE(START,NSAMPLES,'pdf',PDF,'proppdf',PROPPDF,'proprnd',PROPRND)
%   draws NSAMPLES random samples from a target stationary distribution PDF
%   using the Metropolis-Hasting algorithm. 
%       - START is a stucture with fields x0 (colum vector or matrix
%       containing the start value of the Markov chain) and f0 (the logpdf
%       of x0). The number of columns of x0 determines the number of
%       independent markov chains
%       - NSAMPLES is an integer  specifying the number of samples to be
%       generated. 
%       - PDF, PROPPDF, and PROPRND are function handles created using @.
%       PROPPDF defines the proposal distribution density and PROPRND
%       defines the random number generator for the proposal distribution.
%       PDF and PROPRND take one argument as an input and this argument has
%       the same type and size as START. PROPPDF takes two arguments as
%       inputs and both arguments have the same type and size as START. 
%       - SMPL is a column vector or matrix containing the samples. 
%       - If log density function is preferred, 'pdf' and 'proppdf' can be
%       replaced with 'logpdf' and 'logproppdf'. The density functions used
%       in Metropolis-Hasting algorithm are not necessarily normalized.  
%       - If proppdf or logproppdf q(x|y) satisfies q(x|y) = q(y|x), i.e.
%       the proposal distribution is symmetric, MH_SAMPLING_ENGINE implements Random
%       Walk Metropolis-Hasting sampling.  
%       - If proppdf or logproppdf q(x|y) satisfies q(x|y) = q(x), i.e. the
%       proposal distribution is independent of current values, MH_SAMPLING_ENGINE
%       implements Independent Metropolis-Hasting sampling.
%
%   The proposal distribution q(x|y) gives the probability density for
%   choosing x as the next point when y is the current point.
%
%   SMPL = MH_SAMPLING_ENGINE(...,'burnin',K) generate a Markov chain with values
%   between the starting point and the K-th point omitted, but keep points
%   after that. K is a non-negative integer. The default value of K is 0.
%
%   SMPL = MH_SAMPLING_ENGINE(...,'thin',M) generate a markov chain with M-1 out
%   of M values omitted in the generated sequence. M is a positive integer.
%   The default value is 1.
%
%   [SMPL,ACCEPT] = MH_SAMPLING_ENGINE(...) also returns ACCEPT as the acceptance
%   rate of the proposed distribution. ACCEPT is a scalar if a single chain
%   is generated and is a vector if multiple chains are generated.
%
%  Example:
%    Estimate the second order moment of a Gamma distribution using
%    Independent Metropolis-Hasting sampling:
%        alpha = 2.43;
%        beta = 1;
%        pdf = @(x) gampdf(x,alpha,beta);     %target distribution
%        proppdf = @(x,y) gampdf ...          % proposal pdf
%                   (x,floor(alpha),floor(alpha)/alpha);
%        proprnd = @(x) sum ...               % proposal random sampler
%                   (exprnd(floor(alpha)/alpha,floor(alpha),1));
%        nsamples = 5000;
%        smpl = mh_sampling_engine(1,nsamples,'pdf',pdf,'proprnd', ...
%                   proprnd,'proppdf',proppdf);
%        xxhat = cumsum(smpl.^2)./(1:nsamples)';
%        plot(1:nsamples,xxhat)
%
%    Generate random samples from N(0,1) using Random Walk
%    Metropolis-Hasting sampling:
%    delta = .5;
%    npar=10;
%    pdf = @(x) mean(normpdf(x),1);                     %target distribution
%    % proppdf = @(x,y) unifpdf(y-x,-delta,delta);% proposal/Jumping pdf
%    proprnd = @(x)x + 2*rand(npar,size(x,2))*delta - delta;   % proposal random sampler
%    nsamples = 10000;
%    nchains=4;
%    [x,fx,accept] = mh_sampling_engine(rand(npar,nchains),nsamples,...
%     'pdf',pdf,'proprnd',proprnd,'burnin',300);

% modified from Mathworks' mhsample

% parse the information in the name/value pairs
pnames_defaults = {
    'pdf' ,[],@(x)isa(x,'function_handle') && nargin(x)==1
    'logpdf',[],@(x)isa(x,'function_handle') && nargin(x)==1
    'proppdf',[],@(x)isa(x,'function_handle') && nargin(x)==2
    'logproppdf',[],@(x)isa(x,'function_handle') && nargin(x)==2
    'proprnd', [],@(x)isa(x,'function_handle') && nargin(x)==1
    'burnin',0,@(x)isa(x,'double') && isscalar(x) && x==round(x) && x>=0
    'thin',1,@(x)isa(x,'double') && isscalar(x) && x==round(x) && x>=1
    'waitbar_update',[],@(x)isa(x,'function_handle')
    };
if nargin==0
    smpl=pnames_defaults(:,1:3);
    return
end

assert(isstruct(start) && ...
    isfield(start,'x0') && ...
    isfield(start,'f0') && ...
    size(start.x0,2)==size(start.f0,2) && ...
    size(start.f0,1)==1,'bad initial input "start"');

[pdf,logpdf,proppdf,logproppdf, proprnd,burnin,thin,waitbar_update]=...
    utils.miscellaneous.parse_arguments(pnames_defaults,varargin{:});

sym=isempty(proppdf)||isempty(logproppdf);

assert(~(isempty(pdf)&& isempty(logpdf)),...
    'pdf or logpdf for target distribution has to be provided.')

assert(~isempty(proprnd),...
    'proprnd, the random number generator for the proposal distribution has to be provided.')

% error checks for the functional handles
if ~isempty(pdf)
    checkFunErrs('pdf',pdf,start.x0);
end;
if ~isempty(logpdf)
    checkFunErrs('logpdf',logpdf,start.x0);
end;
if ~isempty(proppdf)
    checkFunErrs('proppdf',proppdf,start.x0);
end;
if ~isempty(logproppdf)
    checkFunErrs('logproppdf',logproppdf,start.x0);
end;
if ~isempty(proprnd)
    checkFunErrs('proprnd',proprnd,start.x0);
end;

%error checks for burnin and thin
if (burnin<0) || burnin~=round(burnin)
    error('Bad burnin value (must be a positive integer)');
end
if (thin<=0)|| thin~=round(thin)
    error('Bad thin value (must be a positive integer)');
end

% log density is preferred for numerical stability
if ~isempty(pdf) && isempty(logpdf)
    logpdf = @(x) mylog(pdf(x));
end;
if ~isempty(proppdf) && isempty(logproppdf)
    logproppdf = @(x,y) mylog(proppdf(x,y));
end;
if ~sym
    if any (logpdf(start.x0)==-Inf | logproppdf(start.x0,start.x0) == -Inf)
        error('Bad initial conditions logpdf(x0)=-inf or logproppdf(x0,x0)=-inf')
    end
else
    if any (logpdf(start.x0)==-Inf)
        error('Bad initial conditions logpdf(x0)=-inf')
    end
end;
outclass = superiorfloat(start.x0); % single or double

% Put the replicates dimension second.
[npar,nchain] = size(start.x0);

smpl = zeros([npar,nsamples,nchain],outclass);
fsmpl = zeros([1,nsamples,nchain],outclass);

accept =zeros(1,nchain,outclass);
% Metropolis-Hasting Algorithm.
U = log(rand(nsamples*thin+burnin,nchain));
q_x0_given_y = 0; q_y_given_x0 = 0;
% redo this in case the user makes a mistake like giving the negative of
% the log posterior instead of the log posterior directly. 
start.f0=logpdf(start.x0);
start.funevals=start.funevals+nchain;
for idraw = 1-burnin:nsamples*thin
    y = proprnd(start.x0); % sample from proposal dist'n
    if ~sym
        q_x0_given_y = logproppdf(start.x0,y);
        q_y_given_x0 = logproppdf(y,start.x0);
    end
    fy=logpdf(y);
    start.funevals=start.funevals+nchain;
    % this is a generic formula.
    rho = (fy+q_x0_given_y)-(start.f0+q_y_given_x0);
    
    % Accept or reject the proposal.
    acc = min(rho,0)>=U(idraw+burnin,:);
    start.x0(:,acc) = y(:,acc); % preserves x's shape.
    start.f0(acc)=fy(acc);
    accept = accept+acc;
    if idraw>0 && mod(idraw,thin)==0; % burnin and thin
        smpl(:,idraw/thin,:) = start.x0;
        fsmpl(1,idraw/thin,:) = start.f0;
    end
    if ~isempty(waitbar_update)
        new_messages={
            sprintf('Acceptance rate %0.4f',accept/(idraw+burnin))
            sprintf('thinning factor %0.0f',thin)
            };
        waitbar_update(new_messages{:});
    end
end

% Accept rate can be used to optimize the choice of scale parameters in
% random walk MH sampler. See for example Roberts, Gelman and Gilks (1997).
accept_rate = accept/(nsamples*thin+burnin);

end
%-------------------------------------------------
function  y = mylog(x)
% my log function is to define to avoid the warnings.
y = -Inf(size(x));
y(x>0) = log(x(x>0));
end
%----------------------------------------------------
function checkFunErrs(type,fun,param)
%CHECKFUNERRS Check for errors in evaluation of user-supplied function
if isempty(fun), return; end

try
    switch type
        case {'logproppdf','proppdf'}
            out=fun(param,param);
        case {'logpdf','pdf','proprnd'}
            out=fun(param);
    end
catch ME
    throwAsCaller(ME);
end
switch type
    case 'logproppdf'
        if any(isnan(out))||any(isinf(out) & out>0 )
            error('non finite log proposal density')
        end;
    case 'proppdf'
        if any(~isfinite(out))
            error('non finite proposal density');
        end;
        if any(out<0)
            error('negative proposal density');
        end;
    case 'logpdf'
        if any(isnan(out))||any(isinf(out) & out>0 )
            error('non-finite logpdf');
        end;
    case 'pdf'
        if any(~isfinite(out))
            error('non finite target density');
        end;
        if any(out<0)
            error('negative target density');
        end;
    case 'proprnd'
        if any(~isfinite(out))
            error('non finite proposal random draw');
        end;
end

end
