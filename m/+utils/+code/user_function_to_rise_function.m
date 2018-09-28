function [objective,vargs,is_required]=user_function_to_rise_function(objective)
% INTERNAL FUNCTION: extracts objective function and
% additional input arguments for user-defined functions
%

% signals that information is needed from the model object. This happens
% when the objective is a char or a cell array with a char in its first
% entry and the char starts with a "*" (star). In that case the filtering
% procedure will be called first with the model object and a structure in
% which the user can add some information relevant for his filter.

is_required=false;

vargs={};

if ~isa(objective,'function_handle')
    
    if ~ischar(objective)
        
        if ~isa(objective,'cell')
            
            error('objective function must be a function handle, a cell or a char')
            
        end
        
        if numel(objective)<2
            
            error('when objective function is a cell, it must have have at least two elements')
            
        end
        
        vargs=objective(2:end);
        
        objective=objective{1};
        
    end
    
    if ischar(objective)
        
        objective=objective(~isspace(objective));
        
        is_required=strcmp(objective(1),'*');
        
        if is_required
            
            objective=objective(2:end);
            
        end
        
        objective=str2func(objective);
        
    end
    
end

end