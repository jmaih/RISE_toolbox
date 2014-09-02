function [st,Q,PAI,retcode]=choose_state(st,Q,PAI,y)
% st: state
% Q: transition matrix or ...
% PAI: current updated probabilities
% y: data for current period
retcode=0;
if isempty(st)||isnan(st)
    endogenous_switching=~isempty(Q{2});
    Q0=Q{1};
    % update probabilities
    %---------------------
    if endogenous_switching
        % in the endogenous probability case the configuration of
        % the transition matrix will change
        shadow_transition_matrix=Q{2};
        Vargs=Q{3};
        [Q0,retcode]=utils.code.evaluate_transition_matrices(shadow_transition_matrix,y,Vargs{:});
    end
    if retcode
        return
    end
    PAI=Q0'*PAI;
    csp=[0;cumsum(PAI)];
    st=find(csp>rand,1,'first')-1;
end
end
