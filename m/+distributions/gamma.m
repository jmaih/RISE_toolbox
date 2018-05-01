function varargout=gamma(lowerquantileORmean,upperquantileORstdev,prob,c,d)
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

matlab_gamma_definition=true;

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
    % check the gamma restrictions % 0 < x < inf
    if  ~mean_stdev_flag && (lowerquantileORmean<=0||lowerquantileORmean>=inf)
        error([mfilename,'lowerquantileORmean must satisfy  0 < x < inf'])
    end
    if  ~mean_stdev_flag && (upperquantileORstdev<=0||upperquantileORstdev>=inf)
        error([mfilename,'upperquantileORstdev must satisfy  0 < x < inf'])
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
        aa=(m/s)^2;
        if matlab_gamma_definition
            bb=m/aa;
        else
            bb=aa/m;
        end
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
        if a<=0||b<=0
            moments=nan(2,1);
        else
            if matlab_gamma_definition
                moments=[a*b;sqrt(a)*b];
            else
                moments=[a/b;sqrt(a)/b];
            end
        end
    end

    function lpdf=log_density(theta,a,b,c,d)
        if matlab_gamma_definition
            lpdf=-a.*log(b)-gammaln(a)+(a-1).*log(theta)-theta./b;
        else
            lpdf=a.*log(b)-gammaln(a)+(a-1).*log(theta)-b.*theta;
        end
        target=theta>0;
        lpdf(~target)=-inf;
    end
    function cdf=cumulative_density_function(theta,a,b,c,d)
        if matlab_gamma_definition
            cdf=gamcdf(theta,a,b);
        else
            cdf=gamcdf(theta,a,1/b);
        end
    end
    function icdf=inverse_cdf(u,a,b,c,d)
        if nargin<5
            d=1;
            if nargin<4
                c=0;
            end
        end
        if matlab_gamma_definition
            icdf=gaminv(u,a,b);
        else
            icdf=gaminv(u,a,1/b);
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
        d=inverse_cdf(rand(n,1),a,b,c,d);
    end

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