function [T,F,retcode]=msre_aim(~,Aplus,A0,Aminus,~,~,hh)

if hh>1
    A0=expand_array(A0,1);
    Aminus=expand_array(Aminus,1);
    Aplus=expand_array(Aplus,1);
end

[T,~,~,retcode]=dsge_solve_aim(Aplus,A0,Aminus,[],[],true);

if retcode==0
    T=collapse_array(T,1,hh);
	F=0;
else
    F=nan;
end
