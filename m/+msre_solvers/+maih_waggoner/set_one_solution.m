function [S,retcode]=set_one_solution(X,Sproto,doall,refine_data)
% INTERNAL FUNCTION
%

S=Sproto;

[n,~,h]=size(S);

retcode=0;

for ireg=1:h
    
    row_range=(ireg-1)*n+(1:n);
    
    S(:,:,ireg)=X(row_range,row_range);
    
end

if ~isempty(refine_data)
    
    [S,retcode]=msre_solvers.maih_waggoner.refine(S,refine_data{:});
    
end

if ~retcode
    
    [S,retcode]=msre_solvers.maih_waggoner.realize_solution(S,doall);
    
end

end
