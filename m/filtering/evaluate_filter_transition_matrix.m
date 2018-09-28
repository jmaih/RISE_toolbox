function [Q,retcode]=evaluate_filter_transition_matrix(Qfunc,att,PAItt,ergodic)

if nargin<4
    
    ergodic=false;
    
end

celltype=iscell(att);

if ergodic
    
    if celltype
        
        att_all=att{1}*PAItt(1);
        
        for rt=2:numel(PAItt)
            
            att_all=att_all+att{rt}*PAItt(rt);
            
        end
        
    else
        
        att_all=att(:,1)*PAItt(1);
        
        for rt=2:numel(PAItt)
            
            att_all=att_all+att(:,rt)*PAItt(rt);
            
        end
        
    end
    
elseif celltype
    
    att_all=cell2mat(att);
    
end

[Q,retcode]=Qfunc(att_all);

end