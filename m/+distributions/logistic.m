function varargout=logistic(lowerquantileORmean,upperquantileORstdev,prob,c,d)
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

% The problem to solve is the following:
% find a and b such that probability(lowerquantileORmean < x_a_b < upperquantileORstdev)=prob, with c and d
% the lower and the upper bounds of the distribution with hyperparameters a
% and b.
% This function returns 2 types or outputs depending on whether there
% inputs or not.
% CASE 1: with inputs lowerquantileORmean,upperquantileORstdev,prob, and possibly c and d, the function
% returns {a,b,moments,fval,space}, where
%                        'a' is the first hyperparameter
%                        'b' is the second hyperparameter,
%                        'moments' is a structure with fields: mean and sd
%                        'fval' is a measure of the convergence achieved in
%                        the search for the hyperparameters
%                        'space' is the domain for hyperparameters 'a' and 'b'
% CASE 2: with no inputs, the function returns {lpdfn,cdfn,icdfn,rndfn}
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
    % check whether the inputs are lower_quantile and upper_quantile or
    % mean and standard deviations.
    mean_stdev_flag=isempty(prob)||isnan(prob);
    % check the logistic restrictions % -inf < x < inf
    if ~isfinite(lowerquantileORmean)
        error([mfilename,'lowerquantileORmean must be finite'])
    end
    if ~isfinite(upperquantileORstdev)
        error([mfilename,'upperquantileORstdev must be finite'])
    end
    if  ~mean_stdev_flag && (upperquantileORstdev<=lowerquantileORmean)
        error([mfilename,':: upper bound cannot be less than or equal to lower bound'])
    end
    % find the hyperparameters space
    space=hyperparameters();
    if mean_stdev_flag
        [a,b,fval]=moments_2_hyperparameters(lowerquantileORmean,upperquantileORstdev,c,d);
        moments=[lowerquantileORmean,upperquantileORstdev];
    else
        [ab,fval,retcode]=distributions.find_hyperparameters(space,cdfn,lowerquantileORmean,upperquantileORstdev,prob,c,d);
        if retcode
            error([mfilename,':: could not find the hyperparameters of the distribution'])
        end
        a=ab(1);
        b=ab(2);
        moments=hyperparameters_2_moments(a,b,c,d);
    end
    
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

function [a,b,fval]=moments_2_hyperparameters(m,s,c,d)
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
aa=m;
bb=s*sqrt(3)/pi;
[ab,fval]=fsolve(@objective,[aa;bb],optimset('display','off'),c,d,m,s);
fval=norm(fval);
a=ab(1);
b=ab(2);
    function res=objective(x,c,d,m,s)
        res=hyperparameters_2_moments(x(1),x(2),c,d)-[m,s]';
        if any(isnan(res))
            res=1e+8*ones(2,1);
        end
    end
end

function moments=hyperparameters_2_moments(a,b,~,~)
if b<=0
    moments=nan(2,1);
else
    moments=[a;b*pi/sqrt(3)];
end
end

function icdf=inverse_cdf(u,a,b,~,~)
icdf=a+b*log(u./(1-u));
end

function d=draws(a,b,n,c,d)
if nargin<5
    d=1;
    if nargin<4
        c=0;
        if nargin<3
            n=numel(b);
        end
    end
end
d=inverse_cdf(rand(n,1),a,b,c,d);
end

function lpdf=log_density(theta,a,b,~,~)
lpdf=(theta-a)./b-log(b)-2*log(1+exp((theta-a)./b));
end
function cdf=cumulative_density_function(theta,a,b,~,~)
cdf=1./(1+exp(-(theta-a)./b));
end

function [violation,space]=hyperparameters(hyper)
% -inf<a<inf, b>0
space=[[-inf;1e-12],inf(2,1)];
if nargin==0||isempty(hyper)
    violation=space;
else
    violation=any(hyper(:)<space(:,1))||any(hyper(:)>space(:,2));
end
end