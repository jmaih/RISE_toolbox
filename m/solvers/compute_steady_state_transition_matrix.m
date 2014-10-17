function [Q,retcode]=compute_steady_state_transition_matrix(trans_mat_func,...
    ss_0,pp0,def_0,exo_nbr)

% compute_steady_state_transition_matrix computes the transition matrix a
% the steady state
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

s0=1;s1=1;
xss=zeros(exo_nbr,1);
ys_=ss_0;
[Q,retcode]= utils.code.evaluate_transition_matrices(trans_mat_func,...
    ys_,xss,ss_0,pp0,[],def_0,s0,s1);
end

% function [Q,retcode]=compute_steady_state_transition_matrix(obj,ss_0,pp)
% s0=1;s1=1;
% def_=obj.solution.definitions{1};
% x=zeros(sum(obj.exogenous.number),1);
% ys_=ss_0(:,1);
% [Q,retcode]= utils.code.evaluate_transition_matrices(...
%     obj.routines.transition_matrix,...
%     ys_,x,ss_0(:,1),pp(:,1),[],def_,s0,s1);
% end
