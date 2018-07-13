function Qfunc=prepare_transition_routine(obj)
% INTERNAL FUNCTION
%

Qfunc=memoizer(obj.solution.transition_matrices.Q);

end

function Qfunc=memoizer(Q0)
Qfunc=@engine;
    function [Q,retcode]=engine(~)
        retcode=0;
        Q=Q0;
    end
end
