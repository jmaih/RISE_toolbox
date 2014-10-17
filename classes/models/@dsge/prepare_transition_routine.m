function Qfunc=prepare_transition_routine(obj)
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


if obj.is_endogenous_switching_model
    M=obj.parameter_values;
    % transition matrix should be invariant. And so, hopefully, the first
    % argument to Q{2} could be the updated data in any state.
    % order of the input arguments is y,x,ss,param,sparam,def,s0,s1 
    % remaining arguments of Q{2} after the first one
    %-------------------------------------------------
    defs=mean([obj.solution.definitions{:}],2);
    Vargs={[],obj.solution.ss{1},mean(M,2),[],defs,[],[]};%Q{3}={[],obj.solution.ss{1},mean(M,2),[],[],[],[]};
    Qfunc=@do_evaluation;
else
    Qfunc=prepare_transition_routine@rise_generic(obj);
end

    function [Q,retcode]=do_evaluation(y)
        [Qall,retcode]=utils.code.evaluate_transition_matrices(obj.routines.transition_matrix,y,Vargs{:});
        Q=Qall.Q;
    end
end
