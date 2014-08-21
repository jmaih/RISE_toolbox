function obj=apply_linear_restrictions_to_parameters(obj)

if ~isempty(obj.routines.linear_restrictions)
    linrest=obj.routines.linear_restrictions;
    M=obj.parameter_values;
    for irest=1:size(linrest,1)
        row=linrest{irest,1};
        cols=linrest{irest,2};
        M(row,cols)=linrest{irest,3}(M);
    end
    obj.parameter_values=M;
end