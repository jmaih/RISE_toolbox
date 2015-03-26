function simdiags=simulation_diagnostics(obj,simulation_folder)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    simdiags=struct();
    return
end

if nargin<2
    simulation_folder=[obj.options.results_folder,filesep,'simulations'];
end
is_saved_to_disk=ischar(simulation_folder);
if is_saved_to_disk
    W = what(simulation_folder);
    W=W.mat;
    locs=find(strncmp('chain_',W,6));
    if isempty(locs)
        error([mfilename,':: no simulations found'])
    end
    W=W(locs);
elseif isstruct(simulation_folder)
    W=fieldnames(simulation_folder);
else
    error('wrong specification of input')
end

% first determine the number of chains and the number  of matrices in each
% chain.
number_of_parallel_chains=obj.options.mcmc_number_of_parallel_chains;

incmnt=500;
npar=numel(obj.estimation.priors);
recursive_mean=nan(npar,incmnt,number_of_parallel_chains);
recursive_variance=nan(npar,incmnt,number_of_parallel_chains);
AllParams=nan(npar,incmnt);
offset=0;
theta_mean=0;
V0=0;
for pc=1:number_of_parallel_chains
    matrices=regexp(W,['(?<!\w)chain_',sprintf('%0.0f',pc),'_\d+(?!\w)'],'match');
    matrices=[matrices{:}];
    number_of_matrices=numel(matrices);
    iter=0;
    for m=1:number_of_matrices
        this_matrix=['chain_',sprintf('%0.0f',pc),'_',sprintf('%0.0f',m)];
        if is_saved_to_disk
            tmp=load([simulation_folder,filesep,this_matrix]);
        else
            tmp=simulation_folder.(this_matrix);
        end
        Params=tmp.Params;
        nvals=size(Params,2);
        if offset+nvals>size(AllParams,2);
            incmnt=max(incmnt,nvals);
            AllParams(:,offset+(1:incmnt))=nan;
        end
        AllParams(:,offset+(1:nvals))=Params;
        offset=offset+nvals;
        for ii=1:nvals
            iter=iter+1;
            [theta_mean,V0]=utils.moments.recursive(theta_mean,V0,Params(:,ii),iter);
            if iter==size(recursive_mean,2)
                recursive_mean(:,end+(1:incmnt),:)=nan;
                recursive_variance(:,end+(1:incmnt),:)=nan;
            end
            recursive_mean(:,iter,pc)=theta_mean;
            recursive_variance(:,iter,pc)=diag(V0);
        end
    end
end
AllParams=AllParams(:,1:offset);
recursive_mean=recursive_mean(:,1:iter,:);
recursive_variance=recursive_variance(:,1:iter,:);
[R,B,W]=potential_scale_reduction(AllParams);
clear AllParams
simdiags=struct();
for iparam=1:numel(obj.estimation.priors)
    pname=obj.estimation.priors(iparam).name;
    simdiags.univariate.recursive_mean.(pname)=squeeze(recursive_mean(iparam,:,:));
    simdiags.univariate.recursive_variance.(pname)=squeeze(recursive_variance(iparam,:,:));
    simdiags.univariate.PSRF.(pname)=R(iparam,:);
    simdiags.univariate.between_variance.(pname)=B(iparam,:);
    simdiags.univariate.within_variance.(pname)=W(iparam,:);
end
clear recursive_mean recursive_variance
