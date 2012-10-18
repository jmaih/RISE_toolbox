function varargout=weibull(plb,pub,prob,c,d)
% The problem to solve is the following:
% find a and b such that probability(plb < x_a_b < pub)=prob, with c and d
% the lower and the upper bounds of the distribution with hyperparameters a
% and b.
% This function returns 2 types or outputs depending on whether there
% inputs or not.
% CASE 1: with inputs plb,pub,prob, and possibly c and d, the function
% returns {a,b,moments,fval,space}, where
%                        'a' is the first hyperparameter
%                        'b' is the second hyperparameter,
%                        'moments' is a structure with fields: mean and sd
%                        'fval' is a measure of the convergence achieved in
%                        the search for the hyperparameters
%                        'space' is the domain for hyperparameters 'a' and 'b'
% CASE 2: with no inputs, the function returns {lpdfn,cdfn,icdfn,rndfn,m2h,h2m}
% where
%                  'lpdfn(theta,a,b,c,d)' is the log density function
%                  'cdfn(theta,a,b,c,d)' is the cumulative density function
%                  'icdfn(u,a,b,c,d)' is the inverse of the cdf
%                  'rndfn(a,b,n,c,d)' is the function that computes draws
%                  'm2h(m,s,c,d)' returns the hyperparameters given m and s
%                  'h2m(a,b,c,d)' returns the moments given the hyperparams
%
cdfn=@cumulative_density_function;

hyperparameter_mode=nargin>0;
if hyperparameter_mode
    if nargin<5
        d=1;
        if nargin<4
            c=0;
            if nargin<3
                prob=[];
                if nargin<2
                    error([mfilename,':: wrong number of arguments. Must be 0,2,3,4 or 5'])
                end
            end
        end
    end
    % check the weibull restrictions % 0 < x < inf
    if plb<=0||plb>=inf
        error([mfilename,'plb must satisfy  0 < x < inf'])
    end
    if pub<=0||pub>=inf
        error([mfilename,'pub must satisfy  0 < x < inf'])
    end
    if pub<=plb
        error([mfilename,':: upper bound cannot be less than or equal to lower bound'])
    end
    % find the hyperparameters space
    space=hyperparameters();
    [ab,fval,retcode]=distributions.find_hyperparameters(space,cdfn,plb,pub,prob,c,d);
    
    if retcode
        error([mfilename,':: could not find the hyperparameters of the distribution'])
    end
    
    a=ab(1);
    b=ab(2);
    moments=hyperparameters_2_moments(a,b,c,d);
    moments=struct('mean',moments(1),'sd',moments(2));
    varargout={a,b,moments,fval,space};
else
    icdfn=@inverse_cdf;
    lpdfn=@log_density;
    rndfn=@draws;
    m2h=@moments_2_hyperparameters;
    h2m=@hyperparameters_2_moments;
    varargout={lpdfn,cdfn,icdfn,rndfn,m2h,h2m};
end

end

function [a,b]=moments_2_hyperparameters(m,s,c,d)
if nargin<4
    d=[];
    if nargin<3
        c=[];
    end
end
if m<=0
    error([mfilename,':: m must be >0'])
end
if s<=0
    error([mfilename,':: s must be >0'])
end
ab=fsolve(@objective,rand(2,1),optimset('display','iter'),c,d,m,s);
a=ab(1);
b=ab(2);
    function res=objective(x,c,d,m,s)
        res=hyperparameters_2_moments(x(1),x(2),c,d)-[m,s]';
        if any(isnan(res))
            res=1e+8*ones(2,1);
        end
    end
end

function moments=hyperparameters_2_moments(a,b,c,d)
if a<=0||b<=0
    moments=nan(2,1);
else
    moments=[b*exp(gammaln(1+1/a)) % mean
        b*(sqrt(exp(gammaln(1+2/a)))-exp(gammaln(1+1/a)))]; % sd
end
end

function d=draws(a,b,n,c,d)
if nargin<5
    d=1;
    if nargin<4
        c=0;
        if nargin<3
            n=1;
        end
    end
end
d=inverse_cdf(rand(n,1),a,b,c,d);
end

function icdf=inverse_cdf(u,a,b,~,~)
icdf=b.*(-log(1-u)).^(1./a);
end

function lpdf=log_density(theta,a,b,~,~)
lpdf=log(a./b)-(a-1).*log(theta./b)-(theta./b).^a;
target=theta>0 & theta<inf;
lpdf(~target)=-inf;
end

function cdf=cumulative_density_function(theta,a,b,~,~)
cdf=1-exp(-(theta./b).^a);
end

function [violation,space]=hyperparameters(hyper)
% a>0, b>0
space=[1e-12*ones(2,1),inf(2,1)];
if nargin==0||isempty(hyper)
    violation=space;
else
    violation=any(hyper(:)<space(:,1))||any(hyper(:)>space(:,2));
end
end
