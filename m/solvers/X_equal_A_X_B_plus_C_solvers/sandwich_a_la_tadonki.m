function [V,retcode]=sandwich_a_la_tadonki(A,B,C,solve_options)
% sandwich_a_la_tadonki attempts to solve the equation V=A*V*B+C
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

%
if nargin<4
    solve_options=[];
end
if isempty(solve_options)
    solve_options=struct('TolFun',1e-9,'MaxIter',7000,'verbose',false);
end

solver=6;
V0=C(:); % <-- V=zeros(size(A));
switch solver
    case 1
        [V,retcode]=transpose_free_quasi_minimum_residual(@I_kron_A_B_V,...
            V0,C(:),solve_options,A,B);
    case 2
        [V,retcode] = gmres(@(v)I_kron_A_B_V(v,A,B),C(:),[],solve_options.TolFun,...
            min(numel(C),solve_options.MaxIter),[],[],V0);
    case 3
        [V,retcode] = tfqmr(@(v)I_kron_A_B_V(v,A,B),C(:),solve_options.TolFun,...
            min(numel(C),solve_options.MaxIter),[],[],V0);
    case 4
        [V,retcode] = bicgstab(@(v)I_kron_A_B_V(v,A,B),C(:),solve_options.TolFun,...
            min(numel(C),solve_options.MaxIter),[],[],V0);
    case 5
        [V,retcode] = bicgstabl(@(v)I_kron_A_B_V(v,A,B),C(:),10*solve_options.TolFun,...
            min(numel(C),solve_options.MaxIter),[],[],V0);
    case 6
        [V,retcode] = cgs(@(v)I_kron_A_B_V(v,A,B),C(:),solve_options.TolFun,...
            min(numel(C),solve_options.MaxIter),[],[],V0);
    otherwise
end
V=reshape(V,size(A));

% check whether the result should be symmetrized
if isequal(A,B') && isequal(C,C')
    V=utils.cov.symmetrize(V);
end

function out=I_kron_A_B_V(v,A,B)
n=size(A,1);

out=v-utils.kronecker.kron_times_vector(B',A,reshape(v,n,n));

