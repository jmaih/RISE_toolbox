function varargout=truncated_normal(lowerquantileORmean,upperquantileORstdev,prob,c,d)
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
lpdfn=@log_density;

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
    if isempty(c)||isnan(c),c=0;end
    if isempty(d)||isnan(d),d=1;end
    % check whether the inputs are lower_quantile and upper_quantile or
    % mean and standard deviations.
    mean_stdev_flag=isempty(prob)||isnan(prob);
    % check the truncated normal restrictions % c <= x <= d
    if ~mean_stdev_flag && (lowerquantileORmean<=c||lowerquantileORmean>=d)
        error([mfilename,'lowerquantileORmean must satisfy  c <= x <= d'])
    end
    if ~mean_stdev_flag && (upperquantileORstdev<=c||upperquantileORstdev>=d)
        error([mfilename,'upperquantileORstdev must satisfy  c <= x <= d'])
    end
    if ~mean_stdev_flag && (upperquantileORstdev<=lowerquantileORmean)
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
    rndfn=@draws;
    m2h=@moments_2_hyperparameters;
    h2m=@hyperparameters_2_moments;
    varargout={lpdfn,cdfn,icdfn,rndfn,m2h,h2m};
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
        aa=rand;
        bb=rand;
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

    function moments=hyperparameters_2_moments(a,b,c,d)
        if b<=0
            moments=nan(2,1);
        else
		    mu=a;
		    sig=b;
		    a1=(c-mu)/sig;
		    a2=(d-mu)/sig;
		    PHI1=.5*(1+erf(a1/sqrt(2)));
		    PHI2=.5*(1+erf(a2/sqrt(2)));
		    D=PHI2-PHI1;
		    
		    pdfn=@(theta,a,b,c,d)exp(lpdfn(theta,a,b,c,d));
		    phi1=pdfn(c,a,b,c,d);
		    phi2=pdfn(d,a,b,c,d);
		    mu_trunc=mu-sig*(phi2-phi1)/D;
		    aa1=1;
		    if isfinite(a1)
		        aa1=a1;
		    end
		    aa2=1;
		    if isfinite(a2)
		        aa2=a2;
		    end
		    var_trunc=sig^2*(1-(aa2*phi2-aa1*phi1)/D-((phi2-phi1)/D)^2);
		    moments=[mu_trunc;sqrt(var_trunc)];
        end
    end
	
    function [violation,space]=hyperparameters(hyper)
        % -inf<a<inf, b>0
        space=[c,d;
            1e-12,inf];
        if nargin==0||isempty(hyper)
            violation=space;
        else
            violation=any(hyper(:)<space(:,1))||any(hyper(:)>space(:,2));
        end
    end

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
if isempty(c)||isnan(c),c=0;end
if isempty(d)||isnan(d),d=1;end
d=inverse_cdf(rand(n,1),a,b,c,d);
end

function icdf=inverse_cdf(u,a,b,c,d)
if nargin<5
    d=1;
    if nargin<4
        c=0;
    end
end
if isempty(c)||isnan(c),c=0;end
if isempty(d)||isnan(d),d=1;end
mu=a;
sig=b;
a1=(c-mu)/sig;
a2=(d-mu)/sig;
PHI1=.5*(1+erf(a1/sqrt(2)));
PHI2=.5*(1+erf(a2/sqrt(2)));
D=PHI2-PHI1;
icdf=a+b.*sqrt(2).*erfinv(2*(D.*u+PHI1)-1);
end

function lpdf=log_density(theta,a,b,c,d)
if isempty(c)||isnan(c),c=0;end
if isempty(d)||isnan(d),d=1;end
mu=a;
sig=b;
a1=(c-mu)/sig;
a2=(d-mu)/sig;
PHI1=.5*(1+erf(a1/sqrt(2)));
PHI2=.5*(1+erf(a2/sqrt(2)));
D=PHI2-PHI1;
lpdf=-log(b*sqrt(2*pi))-.5*((theta-a)./b).^2-log(D);
target=theta>=c&theta<=d;
lpdf(~target)=-inf;
end

function cdf=cumulative_density_function(theta,a,b,c,d)
if isempty(c)||isnan(c),c=0;end
if isempty(d)||isnan(d),d=1;end
a1=(c-a)./b;
a2=(d-a)./b;
PHI1=.5*(1+erf(a1/sqrt(2)));
PHI2=.5*(1+erf(a2/sqrt(2)));
D=PHI2-PHI1;
target=theta>=c&theta<=d;
cdf=(.5*(1+erf(((theta-a)./b)/sqrt(2)))-PHI1)./D;
cdf(~target)=nan;
end


