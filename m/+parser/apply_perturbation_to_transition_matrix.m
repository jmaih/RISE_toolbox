function TM=apply_perturbation_to_transition_matrix(TM,coef)

h=size(TM,1);

if h==1||isempty(coef)
    
    return
    
end

for s0=1:h
    
    the_row=TM(s0,:);
    
    for s1=1:h
        
        if s1==s0
            
            continue
            
        end
        
        TM{s0,s1}=[coef,'*(',TM{s0,s1},')'];
        
    end
    
    TM{s0,s0}=set_diagonal(the_row,s0);
    
end

    function out=set_diagonal(in,s0)
        
        in(s0)=[];
        
        out=cell2mat(strcat(in,'+'));
        
        out=['1-',coef,'*(',out(1:end-1),')'];
        
    end

end