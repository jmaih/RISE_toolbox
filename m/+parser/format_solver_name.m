function solver=format_solver_name(solver)

if isnumeric(solver)
    
    solver=int2str(solver);
    
elseif iscell(solver)
    
    other_args=solver(2:end);
    
    for ii=1:numel(other_args)
        
        if isnumeric(other_args{ii})||islogical(other_args{ii})
            
            if ~isequal(size(other_args{ii}),[1,1])
                
                other_args{ii}='?';
                
            else
                
                other_args{ii}=num2str(other_args{ii});
                
            end
            
        elseif isempty(other_args{ii})
            
            other_args{ii}='empty';
            
        end
        
    end
    
    solver=solver{1};
    
    if ~isempty(other_args)
        
        other_args=cell2mat(strcat(other_args,','));
        
        solver=[solver,'(',other_args(1:end-1),')'];
        
    end
    
end

end