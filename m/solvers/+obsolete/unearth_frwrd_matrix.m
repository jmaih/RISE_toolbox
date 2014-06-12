function Aplus=unearth_frwrd_matrix(Gplus,Q)
% this function extracts Aplus from Gplus and is used for solving models
% using the Waggoner-Zha approach, which states the problem to solve as
% Aplus(st)*X_{t+1}+A0(st)*X_{t}+Aminus(st)*X_{t-1}+B(st)*Et=0 whereas in
% my approach, the problem solved would be
% Aplus(st+1)*X_{t+1}+A0(st)*X_{t}+Aminus(st)*X_{t-1}+B(st)*Et=0   
h=size(Q,1);
Aplus=cell(h,1);
for reg=1:h
    Aplus{reg}=Gplus{reg,reg}/Q(reg,reg);
end
end
