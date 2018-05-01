function [theta_mean,theta_median,V0]=parameters_moments(simulation_folder)
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
number_of_matrices=numel(W);
theta_mean=0;
V0=0;
iter=0;
% approximated median calculated as the median of the medians
%------------------------------------------------------------
theta_median=cell(2,number_of_matrices);
for m=1:number_of_matrices
    if is_saved_to_disk
        tmp=load([simulation_folder,filesep,W{m}]);
    else
        tmp=simulation_folder.(W{m});
    end
    Params=tmp.Params;
    minus_logpost_params=tmp.minus_logpost_params;
    [~,tags]=sort(minus_logpost_params);
    nvals=size(Params,2);
    md=tags(ceil(0.5*nvals));
    theta_median(:,m)={minus_logpost_params(md);Params(:,md)};
    for ii=1:nvals
        iter=iter+1;
        [theta_mean,V0]=utils.moments.recursive(theta_mean,V0,Params(:,ii),iter);
    end
end
[~,tags]=sort(cell2mat(theta_median(1,:)));
md=tags(ceil(0.5*number_of_matrices));
theta_median=theta_median{2,md};