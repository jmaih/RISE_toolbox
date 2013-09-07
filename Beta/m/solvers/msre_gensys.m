function [T,F,retcode]=msre_gensys(~,Aplus,A0,Aminus,~,~,hh)

if hh>1
    A0=expand_array(A0,1);
    Aminus=expand_array(Aminus,1);
    Aplus=expand_array(Aplus,1);
%     if ~isempty(B)
%         B=expand_array(B,2);
%     end
%     C=expand_array(C,3);
    % Recompute the number of equations
%     nn=size(A0,1);
end

[T,~,~,retcode]=dsge_solve_gensys(Aplus,A0,Aminus,[],[],true);

if retcode==0
    T=collapse_array(T,1,hh);
	F=0;
else
    F=nan;
end
