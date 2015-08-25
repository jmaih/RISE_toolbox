function [draw,obj]=draw_parameter(obj,simulation_folder)
% draw_parameter - random parameter draws for RISE model objects. 
%
% Syntax
% -------
% ::
%
%   [draw,obj]=draw_parameter(obj,simulation_folder)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar]: RISE model object
%
% - **simulation_folder** [char|struct]:
%   - char(1): simulation folder : the stored elements should be structures
%   with fields:
%       - **x** : new form
%       - **Params** : legacy
%   - char(2): ['mode'|'prior']: draw from the prior distribution or from a
%   multivariate normal distribution around the mode.
%
% Outputs
% --------
%
% - **draw** [cell]: the first entry is the names of the estimated
% parameters and the second is a vector of drawn parameters. The whole cell
% can be pushed in to a model as obj=set(obj,'parameters',draw).
%
% - **obj** [rise|dsge|rfvar|svar]: RISE model object in which the drawn
% parameter has been pushed.
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

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
                % replace the dirichlet with the beta
                %-------------------------------------
                distribs=strrep(distribs,'dirichlet','beta');
                udistrib=unique(distribs);
                draw=nan(numel(distribs),1);
                for idist=1:numel(udistrib)
                    % the dirichlets are drawn separately
                    %-------------------------------------
                    if strcmp(udistrib{idist},'dirichlet')
                        continue
                    end
                    loc=strcmp(udistrib{idist},distribs);
                    a=vertcat(obj.estimation.priors(loc).a);
                    b=vertcat(obj.estimation.priors(loc).b);
                    [~,~,~,rndfn]=distributions.(udistrib{idist});
                    draw(loc)=rndfn(a,b);
                end
                for id=1:numel(obj.estim_dirichlet)
                    loc=obj.estim_dirichlet(id).location;
                    draw(loc)=obj.estim_dirichlet(id).rndfn(1);
                end
        end
        draw=utils.optim.recenter(draw,lb,ub);
else
    is_saved_to_disk=ischar(simulation_folder);
    if is_saved_to_disk
        W = what(simulation_folder);
        W=W.mat;
        if isempty(W)
            error([mfilename,':: no simulations found'])
        end
        % check for legacy
        %------------------
        locs=find(strncmp('chain_',W,6));
        is_legacy=~isempty(locs);
        if is_legacy
            W=W(locs);
        end
        %--------------------------
        W=strrep(W,'.mat','');
        N=numel(W);
        choice=pick_one_randomly(N);
        this_matrix=W{choice};
        tmp=load([simulation_folder,filesep,this_matrix]);
        if is_legacy
            choice=pick_one_randomly(size(tmp.Params,2));
            draw=tmp.Params(:,choice);
        else
            if ~isfield(tmp,'x')
                error('wrong format for the stored objects to draw from')
            end
            nn=numel(tmp);
            id=pick_one_randomly(nn);
            draw=tmp.(id).x;
        end
    elseif isstruct(simulation_folder)
        N=numel(simulation_folder);
        choice=pick_one_randomly(N);
        draw=simulation_folder(choice).x;
    else
        error('wrong specification of input')
    end
end

% The user may want to set the parameters directly himself and do it is a
% good idea to untransform the parameters for him
x=unstransform_estimates(obj,draw{2});
draw={pnames,x};
if nargout>1
    % assignin estimates will be faster than calling
    % obj=set(obj,'parameters',draw); 
    obj=assign_estimates(obj,x);
end

function choice=pick_one_randomly(N)

if N<=0||~isequal(N,ceil(N))
	error([mfilename,':: input must be a positive integer'])
end
probs=1/N*ones(1,N);
cprob=[0,cumsum(probs)];
cprob(end)=1;
choice=find(cprob>rand,1,'first')-1;
