function [obj,draw]=draw_parameters(obj,type)
% this function draws a parameter and assigns it.

% this function might not work well if we have markov chains with more than
% 2 states. the recenter function is not enough in that case. One also
% needs to find a way to make sure that all the drawn elements in the
% dirichlet parameters, say, sum to 1. Also, the use of the recenter
% function in this case might be called into question since only the
% parameters that fall outside their boundaries are "recentered"
persistent cc mode n lb ub

if isempty(n)
    n=numel(obj.estimated_parameters);
    lb=vertcat(obj.estimated_parameters.lb);
    ub=vertcat(obj.estimated_parameters.ub);
end

switch type
    case 'mode'
        if isempty(cc)
            [cc,pp]=chol(obj.estimation.posterior_maximization.vcov);
            if pp
                error([mfilename,':: covariance matrix of estimated parameters not positive definite'])
            else
                cc=transpose(cc);
            end
            mode=obj.estimation.posterior_maximization.mode;
        end
        draw=mode+cc*randn(n,1);
    case 'posterior'
        tmp=what('simulation_folder');
        tmp=tmp.mat;
        % select the mat file to draw from with equal probability
        choice=pick_one_randomly(numel(tmp));
        % now select the parameter vector
        test=load(['simulation_folder\',tmp{choice}]);
		choice=pick_one_randomly(size(test.Params,2));
        draw=test.Params(:,choice);
    otherwise
        error([mfilename,':: currently only the mode is implemented'])
end

draw=recenter(draw,lb,ub);
obj=assign_estimates(obj,draw);

function choice=pick_one_randomly(N)

if N<=0||~isequal(N,ceil(N))
	error([mfilename,':: input must be a positive integer'])
end
probs=1/N*ones(1,N);
cprob=[0,cumsum(probs)];
cprob(end)=1;
choice=find(cprob>rand,1,'first')-1;
