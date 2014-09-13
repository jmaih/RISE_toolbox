function [A,B]=set_solution_to_companion(obj)
if isempty(obj)
    if nargout>1
        error([mfilename,':: number of output arguments cannot exceed 1 when the object is empty'])
    end
    A=struct();
    return
end

if isempty(obj.solution)
    error('model has not been solved')
end

endo_nbr=obj.endogenous.number(end);
exo_nbr=sum(obj.exogenous.number);
reg_nbr=obj.markov_chains.regimes_number;

[A,~,B]=vartools.resolve(obj.solution,obj.nlags,reg_nbr);
cc=(obj.nlags-1)*endo_nbr;
for ireg=1:reg_nbr
    A{ireg}=[A{ireg}(:,1:obj.nlags*endo_nbr);
        eye(cc),zeros(cc,endo_nbr)];
    B{ireg}=[B{ireg};
        zeros(cc,exo_nbr)];
end

end
