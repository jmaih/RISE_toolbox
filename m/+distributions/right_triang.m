function varargout=right_triang(lowerquantileORmean,upperquantileORstdev,prob,~,~)
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
%                  'lpdfn(theta,a,b)' is the log density function
%                  'cdfn(theta,a,b)' is the cumulative density function
%                  'icdfn(u,a,b)' is the inverse of the cdf
%                  'rndfn(a,b,n)' is the function that computes draws
%                  'm2h(m,s)' returns the hyperparameters given m and s
%                  'h2m(a,b)' returns the moments given the hyperparams
%
cdfn=@cumulative_density_function;

hyperparameter_mode=nargin>0;
if hyperparameter_mode
    if nargin<3
        prob=[];
        if nargin<2
            error([mfilename,':: wrong number of arguments. Must be 0,2,3,4 or 5'])
        end
    end
    if upperquantileORstdev<=lowerquantileORmean
        error([mfilename,':: upper bound cannot be less than or equal to lower bound'])
    end
    % check whether the inputs are lower_quantile and upper_quantile or
    % mean and standard deviations.
    mean_stdev_flag=isempty(prob)||isnan(prob);
    % find the hyperparameters space
    space=hyperparameters();
    if mean_stdev_flag
        [a,b,fval]=moments_2_hyperparameters(lowerquantileORmean,upperquantileORstdev);
        moments=[lowerquantileORmean,upperquantileORstdev];
    else
        [ab,fval,retcode]=distributions.find_hyperparameters(space,cdfn,lowerquantileORmean,upperquantileORstdev,prob);
        if retcode
            error([mfilename,':: could not find the hyperparameters of the distribution'])
        end
        a=ab(1);
        b=ab(2);
        moments=hyperparameters_2_moments(a,b);
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

function [a,b,fval]=moments_2_hyperparameters(m,s,~,~)
if s<=0
    error([mfilename,':: s must be >0'])
end
s2=s*sqrt(2);
a=m-2*s2;
b=a+3*s2;
fval=0;
end

function moments=hyperparameters_2_moments(a,b,~,~)
mm=a;
Ex=(a+mm+b)/3;
Vx=(a.^2+b.^2+mm.^2-a.*b-a.*mm-b.*mm)/18;
moments=[Ex
    sqrt(Vx)];
end

function d=draws(a,b,n,~,~)
if nargin<3
    n=numel(b);
end
d=inverse_cdf(rand(n,1),a,b);
end

function icdf=inverse_cdf(u,a,b,~,~)
mm=a;
icdf=b-sqrt((1-u).*(b-a).*(b-mm));
end

function lpdf=log_density(theta,a,b,~,~)
mm=a; % x is right to the point with the highest density
lpdf=log(2*(b-theta)./((b-a).*(b-mm)));
target=theta>=a & theta<b;
lpdf(~target)=-inf;
end

function cdf=cumulative_density_function(theta,a,b,~,~)
if ~isequal(size(a),size(b))
    error('hyperparameters should be of the same size')
end
mm=a;
cdf=1-(b-theta).^2./((b-a).*(b-mm));
cdf(theta<a)=0;
cdf(theta>b)=1;
end

function [violation,space]=hyperparameters(hyper)
% a>0, b>0
space=[-inf(2,1),inf(2,1)];
if nargin==0||isempty(hyper)
    violation=space;
else
    violation=any(hyper(:)<space(:,1))||any(hyper(:)>space(:,2));
end
end

