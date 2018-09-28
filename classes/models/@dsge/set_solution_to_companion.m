function [A,B,steady_state]=set_solution_to_companion(obj)
% INTERNAL FUNCTION
%

if isempty(obj)

    if nargout>1

        error([mfilename,':: number of output arguments cannot exceed 1 when the object is empty'])

    end

    A=cell(0,4);

    return

end

if isempty(obj.solution)

    error('model has not been solved')

end

[A,B,Qfunc,steady_state]=load_solution(obj,'iov'); %#ok<ASGLU>

end
