function [islands,ids]=clear_duplicates(islands,lb,ub,newdeal)
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

if nargin<4
    newdeal=false;
end
norm_=norm(ub-lb);
[number_of_parameters,MaxNodes]=size(islands);
ids=false(1,MaxNodes);
for i1=1:MaxNodes
    agent_ii=islands(:,i1);
    for jj=i1+1:MaxNodes
        dd=utils.optim.distance(agent_ii,islands(:,jj));
        if dd/norm_<1e-2
            ids(jj)=true;
            if newdeal
                uu=rand(number_of_parameters,1);
                islands(:,jj)=islands(:,jj)+.5*uu+...
                    .25*uu.*randn(number_of_parameters,1);
                islands(:,jj)=utils.optim.recenter(islands(:,jj),lb,ub);
            else
                % select the chromosomes to change randomly
                change=randperm(number_of_parameters);
                change=change(1:round(number_of_parameters/3));
                islands(change,jj)=lb(change)+...
                    rand(numel(change),1).*(ub(change)-lb(change));
            end
        end
    end
end
ids=find(ids);

%% Reliable: 
% bbo_3
% gampc_2: sensitive to crossover prob .1 seems to work ok
% bee : may/may not require as many probes as advodated by the gampc guys
% cmsa: understands the game pretty well, but impoverishment still needs
% tuning...
% studga: holding its own

% %%
% [energy,lb,ub,fname]=cec_2011_problems(112);
% Options=struct('MaxFunEvals',50000*3,...
%     'verbose',50,'penalty',1e+12,...
%     'MaxNodes',90,...
%     'MaxIter',10000);
% samode(energy,[],[],lb,ub,Options);
%{
	C:\Users\Junior\Dropbox\SandboxFiles\Code\MarkovSwitchingDsge\optimizers\clever_algorithms
	C:\Users\Junior\Dropbox\SandboxFiles\Code\MarkovSwitchingDsge\optimizers\cec_2011_recompiled_problems
	C:\Users\Junior\Dropbox\SandboxFiles\doc\GlobalOptimization\CEC_2011\CEC_2011_Spl_Session\CEC_2011_Spl_Session\Probs_1_to_8
	C:\Users\Junior\Dropbox\SandboxFiles\Code
%}

