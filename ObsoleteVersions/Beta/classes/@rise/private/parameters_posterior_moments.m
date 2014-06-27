function [theta_mean,V0,quantiles]=parameters_posterior_moments(simulation_folder)
W = what(simulation_folder);
W=W.mat;
locs=find(strncmp('chain_',W,6));
if isempty(locs)
    error([mfilename,':: no simulations found'])
end
W=W(locs);
number_of_matrices=numel(W);
theta_mean=0;
V0=0;
iter=0;
out_nargs=nargout;
for m=1:number_of_matrices
    tmp=load([simulation_folder,filesep,W{m}]);
    Params=tmp.Params;
    nvals=size(Params,2);
    for ii=1:nvals
        iter=iter+1;
        if out_nargs<2
            theta_mean=rise_moments.recursive_moments(theta_mean,V0,Params(:,ii),iter);
        else
            [theta_mean,V0]=rise_moments.recursive_moments(theta_mean,V0,Params(:,ii),iter);
        end
    end
end

% quantiles
if out_nargs>2
    Nsim=iter;
    npar=size(theta_mean,1);
    quantiles=cell(npar,1);
    quant=[round((50:-5:5)/100*Nsim)',round((50:5:95)/100*Nsim)']';
    for ii=1:npar
        all_vals=nan(Nsim,1);
        iter=0;
        for m=1:number_of_matrices
            tmp=load([simulation_folder,filesep,W{m}]);
            Params=tmp.Params(ii,:);
            nvals=size(Params,2);
            all_vals(iter+(1:nvals))=Params(:);
            iter=iter+nvals;
        end
        quantiles{ii}=all_vals(quant);
    end
end
