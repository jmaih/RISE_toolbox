function [objective,vargs]=user_function_to_rise_function(objective)
% extract_objective_and_varargin - extracts objective function and
% additional input arguments for user-defined functions
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

vargs={};

if ~isa(objective,'function_handle')
    
    if ischar(objective)
        
        objective=str2func(objective);
        
    else
        
        if ~isa(objective,'cell')
            
            error('objective function must be a function handle, a cell or a char')
            
        end
        
        if ischar(objective{1})
            
            objective{1}=str2func(objective{1});
            
        end
        
        if numel(objective)<2
            
            error('when objective function is a cell, it must have have at least two elements')
            
        end
        
        vargs=objective(2:end);
        
        objective=objective{1};
        
    end
    
end

end