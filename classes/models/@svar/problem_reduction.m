function [T,Q,n,h]=problem_reduction(obj)

T=load_solution(obj,'ov',false);

Q=obj.solution.transition_matrices.Q;

h=numel(T);

n=size(T{1},1);

for ireg=1:h
    T{ireg}=T{ireg}(:,1:n);
end

end