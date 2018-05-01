function [A,B,steady_state]=set_solution_to_companion(obj)
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

[A,B]=load_solution(obj,'iov');

steady_state=obj.solution.ss;

end
