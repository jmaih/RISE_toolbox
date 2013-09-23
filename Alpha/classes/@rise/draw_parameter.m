function [draw,obj]=draw_parameter(obj,simulation_folder)
% outputed are the draw and the rise object in which the draw has been
% pushed.
if isempty(obj)
    draw=struct();
    return
end

if nargin<2
    simulation_folder=obj.folders_paths.simulations;
end

pnames={obj.estimation.priors.name};
if ischar(simulation_folder) && ismember(simulation_folder,{'mode','prior'})
        % this case might not work well if we have markov chains with more
        % than 2 states. the recenter function is not enough in that case.
        % One also needs to find a way to make sure that all the drawn
        % elements in the dirichlet parameters, say, sum to 1. Also, the
        % use of the recenter function in this case might be called into
        % question since only the parameters that fall outside their
        % boundaries are "recentered" 
        n=numel(obj.estimation.priors);
        lb=vertcat(obj.estimation.priors.lower_bound);
        ub=vertcat(obj.estimation.priors.upper_bound);
        switch simulation_folder
            case 'mode'
                [cc,pp]=chol(obj.estimation.posterior_maximization.vcov);
                if pp
                    error([mfilename,':: covariance matrix of estimated parameters not positive definite'])
                end
                xmode=obj.estimation.posterior_maximization.mode;
                draw=xmode+transpose(cc)*randn(n,1);
            case 'prior'
                distribs={obj.estimation.priors.prior_distrib};
                udistrib=unique(distribs);
                draw=nan(numel(distribs),1);
                for idist=1:numel(udistrib)
                    loc=strcmp(udistrib{idist},distribs);
                    a=vertcat(obj.estimation.priors(loc).a);
                    b=vertcat(obj.estimation.priors(loc).b);
                    [~,~,~,rndfn]=distributions.(udistrib{idist});
                    draw(loc)=rndfn(a,b);
                end
        end
        draw=recenter(draw,lb,ub);
else
    is_saved_to_disk=ischar(simulation_folder) && ~strcmp(simulation_folder,'mode');
    if is_saved_to_disk
        W = what(simulation_folder);
        W=W.mat;
        locs=find(strncmp('chain_',W,6));
        if isempty(locs)
            error([mfilename,':: no simulations found'])
        end
        W=strrep(W(locs),'.mat','');
    elseif isstruct(simulation_folder)
        W=fieldnames(simulation_folder);
    else
        error('wrong specification of input')
    end
    % select the mat file to draw from with equal probability
    %--------------------------------------------------------
    choice=pick_one_randomly(numel(W));
    this_matrix=W{choice};
    if is_saved_to_disk
        tmp=load([simulation_folder,filesep,this_matrix]);
    else
        tmp=simulation_folder.(this_matrix);
    end
    % now select the parameter vector
    %--------------------------------
    choice=pick_one_randomly(size(tmp.Params,2));
    draw=tmp.Params(:,choice);
end

draw={pnames,draw};
if nargout>1
    obj=set(obj,'parameters',draw);
end

function choice=pick_one_randomly(N)

if N<=0||~isequal(N,ceil(N))
	error([mfilename,':: input must be a positive integer'])
end
probs=1/N*ones(1,N);
cprob=[0,cumsum(probs)];
cprob(end)=1;
choice=find(cprob>rand,1,'first')-1;
