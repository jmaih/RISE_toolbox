% SCHUR  Schur decomposition.
%    [U,T] = SCHUR(X) produces a quasitriangular Schur matrix T and
%    a unitary matrix U so that X = U*T*U' and U'*U = EYE(SIZE(U)).
%    X must be square.
% 
%    T = SCHUR(X) returns just the Schur matrix T.
% 
%    If X is real, two different decompositions are available.
%    SCHUR(X,'real') has the real eigenvalues on the diagonal and the
%    complex eigenvalues in 2-by-2 blocks on the diagonal.
%    SCHUR(X,'complex') is triangular and is complex if X has complex
%    eigenvalues.  SCHUR(X,'real') is the default.
% 
%    If X is complex, the complex Schur form is returned in matrix T.
%    The complex Schur form is upper triangular with the eigenvalues
%    of X on the diagonal. The second input is ignored in this case.
% 
%    See RSF2CSF to convert from Real to Complex Schur form.
% 
%    See also ORDSCHUR, QZ, RSF2CSF.
%
%    Documentation for schur
%       doc schur
%
%