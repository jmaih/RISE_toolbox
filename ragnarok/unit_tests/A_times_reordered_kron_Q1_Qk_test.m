nr=7;

nc=10;

nf=20;

A=rand(nr,nc);

f=rand(nf,nr);

reorder_rows=randperm(nr);

reorder_mat(reorder_rows)=1:nr;

M=speye(nr);M=M(:,reorder_mat);

fA=f*A(reorder_rows,:);

fMA=f*M*A;

max(max(abs(fA-fMA)))

%%
clc
n=struct();
n.v=10;
n.z=7;
n.f=22;

Q1=rand(n.v^2,n.z);

Q2=rand(n.v,n.z^2);

fvvv_=rand(n.f,n.v^3);


matsizes=[
    n.v,n.z
    n.v,n.z
    n.v,n.z
    ];

orders={[2,1,3],[2,3,1]};%[1,2,3], First added automatically

options=[];

fQ1=utils.kronecker.A_times_reordered_kron_Q1_Qk(fvvv_,matsizes,orders,options,Q1,Q2);

% obj.fv_kron_Q1_Qk(W,matsizes,orders,options,varargin);
% kron_Q1_Qk_times_A(obj.debug,a1z_,a1z_,EU2);
% A_times_kron_Q1_Qk(obj.debug,varargin)

%%

Q12=utils.kronecker.sum_permutations(kron(Q1,Q2),matsizes,options,orders{:});

fQ2=fvvv_*Q12;

%%

max(max(abs(fQ1-fQ2)))