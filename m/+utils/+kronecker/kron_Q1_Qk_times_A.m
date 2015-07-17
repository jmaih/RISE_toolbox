function C=kron_Q1_Qk_times_A(A,varargin)
% kron_Q1_Qk_times_A -- computes (Q1*Q2*...*Qk)A where * denotes the
% kronecker product
%
% Syntax
% -------
% ::
%
% - C=kron_Q1_Qk_times_A(A,Q1,Q2,...,Qk)
%
% Inputs
% -------
%
% - **A** [matrix]: matrix on the left-hand side
%
% - **Qi** [matrix]: matrix on the kronecker block
%
% Outputs
% --------
%
% - **C** [matrix]: result
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: A_times_kron_Q1_Qk

% transpose all inputs
A=A.';
for iarg=1:length(varargin)
    varargin{iarg}=varargin{iarg}.';
end

C=utils.kronecker.A_times_kron_Q1_Qk(A,varargin{:});

% transpose final result
C=C.';

end