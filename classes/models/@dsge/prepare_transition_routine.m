function Qfunc=prepare_transition_routine(obj)
% INTERNAL FUNCTION
%

if obj.is_endogenous_switching_model
    
    M=obj.parameter_values;
    
    % transition matrix should be invariant. And so, hopefully, the first
    % argument to Q{2} could be the updated data in any state.
    % order of the input arguments is y,x,ss,param,def,s0,s1 
    % remaining arguments of Q{2} after the first one
    %-------------------------------------------------
    defs=cell2mat(obj.solution.definitions);
    
    Vargs={[],cell2mat(obj.solution.ss),M,defs,[],[]};
    
    Qfunc=memoizer(obj.routines.transition_matrix,Vargs{:});
    
else
    
    Qfunc=prepare_transition_routine@generic(obj);
    
end

end

function Qfunc=memoizer(routine,varargin)

Qfunc=@engine;

    function [Q,retcode]=engine(y)
        % exponentiation of logvars already done in load solution, which
        % calls this function
        %------------------------------------------------------------------
        [Qall,retcode]=utils.code.evaluate_transition_matrices(routine,y,...
            varargin{:});
        
        Q=Qall.Q;
        
    end

end
