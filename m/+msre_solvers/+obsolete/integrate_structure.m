function [Gplus,A0,Aminus,B,C]=integrate_structure(structural_matrices,h)
Aminus=cell(h,1);
A0=cell(h,1);
B=cell(h,1);
Gplus=structural_matrices.Gp;
C=structural_matrices.Gt;
for s0=1:h
    a0=0;am=0;b=0;
    for s1=1:h
        a0=a0+structural_matrices.Gc{s0,s1};
        am=am+structural_matrices.Gm{s0,s1};
        b=b+structural_matrices.Ge{s0,s1};
    end
    A0{s0}=a0;
    Aminus{s0}=am;
    B{s0}=b;
end
end
