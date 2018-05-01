function C=kron_I_A_times_B(A,B,n)
% H1 line
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

% computes C=kron(In,A)*B

[qn,m]=size(B);
[p,q]=size(A);
if qn/n~=q
    error('matrices sizes inconsistent')
end

C=reshape(A*reshape(B,q,n*m),p*n,m);

end
