function [A,B]=set_solution_to_companion(obj)
A=obj.solution.m_x;
B=obj.solution.m_e;
if isa(obj,'svar')
    endo_nbr=obj.endogenous.number(end);
    exo_nbr=sum(obj.exogenous.number);
    cc=(obj.nlags-1)*endo_nbr;
    for ireg=1:numel(A)
        A{ireg}=[A{ireg}(:,1:obj.nlags*endo_nbr);
            eye(cc),zeros(cc,endo_nbr)];
        B{ireg}=[B{ireg};
            zeros(cc,exo_nbr)];
    end
end
end
