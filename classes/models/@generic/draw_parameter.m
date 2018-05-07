function [draw,obj]=draw_parameter(obj,simulation_folder)
% draw_parameter - random parameter draws for RISE model objects.
%
% ::
%
%
%   [draw,obj]=draw_parameter(obj,simulation_folder)
%
% Args:
%
%    - **obj** [rise|dsge|rfvar|svar]: RISE model object
%
%    - **simulation_folder** [char|struct]:
%      - char(1): simulation folder : the stored elements should be structures
%      with fields:
%          - **x** : parameter vectors
%          - **f** : value of minus(log posterior kernel)
%      - char(2): ['mode'|'prior']: draw from the prior distribution or from a
%      multivariate normal distribution around the mode.
%
% Returns:
%    :
%
%    - **draw** [cell]: the first entry is the names of the estimated
%    parameters and the second is a vector of drawn parameters. The whole cell
%    can be pushed in to a model as obj=set(obj,'parameters',draw).
%
%    - **obj** [rise|dsge|rfvar|svar]: RISE model object in which the drawn
%    parameter has been pushed.
%
% Note:
%
% Example:
%
%    See also:

if isempty(obj)
    
    draw=cell(0,4);
    
    return
    
end

if nargin<2

    simulation_folder=[];

end

if isempty(simulation_folder)

    simulation_folder=obj.folders_paths.simulations;

end

pnames={obj.estimation.priors.name};

need_untransform=true;

if ischar(simulation_folder) && ismember(simulation_folder,{'mode','prior'})
    
    switch simulation_folder
        
        case 'mode'
            
            lb=vertcat(obj.estimation.priors.lower_bound);
            
            ub=vertcat(obj.estimation.priors.upper_bound);
            
            xmode=obj.estimation.posterior_maximization.mode;
            
            vcov=obj.estimation.posterior_maximization.vcov;
            % Draw from the short side and then untransform. The dirichlet
            % distributions should be fine.
            [obj,xmode,lb_short,ub_short,vcov]=transform_parameters(obj,xmode,lb,ub,vcov);
            
            [cc,pp]=chol(vcov);
            
            if pp
            
                error([mfilename,':: covariance matrix of estimated parameters not positive definite'])
            
            end
            
            n=numel(xmode);
            
            draw=xmode+transpose(cc)*randn(n,1);
        
            draw=utils.optim.recenter(draw,lb_short,ub_short);
        
        case 'prior'
            
            need_untransform=false;
            
            distribs={obj.estimation.priors.prior_distrib};
            
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
    
else
    
    is_saved_to_disk=ischar(simulation_folder);
    
    if is_saved_to_disk
        
        W = what(simulation_folder);
        
        W=W.mat;
        
        if isempty(W)
        
            error([mfilename,':: no simulations found'])
        
        end
        
        W=strrep(W,'.mat','');
        
        N=numel(W);
        
        choice=randi(N);
        
        this_matrix=W{choice};
        
        tmp=load([simulation_folder,filesep,this_matrix]);
        
        if isfield(tmp,'pop')
        
            tmp=tmp.pop;
        
        end
        
        if ~isfield(tmp,'x')
        
            error('wrong format for the stored objects to draw from')
        
        end
        
        nn=numel(tmp);
        
        id=randi(nn);
    
        draw=tmp(id).x;
    
    elseif isstruct(simulation_folder)
        
        N=numel(simulation_folder);
        
        choice=randi(N);
    
        draw=simulation_folder(choice).x;
    
    else
        
        error('wrong specification of input')

    end
    
end

if need_untransform
    % The user may want to set the parameters directly himself and do it is a
    % good idea to untransform the parameters for him
    x=unstransform_parameters(obj,draw);

else
    % we need to apply linear restrictions if necessary otherwise the
    % linear restrictions are not enforced for the prior.
    x=back_and_forth(draw);

end

draw={pnames,x};

if nargout>1
    % assignin estimates will be faster than calling
    % obj=set(obj,'parameters',draw);
    obj=assign_estimates(obj,x);

end

    function x=back_and_forth(draw)
        
        if isempty(obj(1).linear_restrictions_data)
        
            obj=setup_restrictions(obj);
        
        end
        
        draw=obj(1).linear_restrictions_data.a2tilde_func(draw);
    
        x=obj(1).linear_restrictions_data.a_func(draw);

    end

end
