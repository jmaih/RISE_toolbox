function [P0,retcode]=diffuse_lyapunov_equation(T,RQR,diffuse_factor,...
check_unit_roots,diffuse_all)
% diffuse_lyapunov_equation attempts to solve V=T*V*T'+RQR under
% nonstationarity
%
% ::
%
%   [P0,retcode]=diffuse_lyapunov_equation(T,RQR)
%   [P0,retcode]=diffuse_lyapunov_equation(T,RQR,diffuse_factor)
%   [P0,retcode]=diffuse_lyapunov_equation(T,RQR,diffuse_factor,check_unit_roots)
%   [P0,retcode]=diffuse_lyapunov_equation(T,RQR,diffuse_factor,check_unit_roots,diffuse_all)
%
% Args:
%    - T :
%    - RQR :
%    - diffuse_factor : [ positive scalar | {10} ]
%    - check_unit_roots : [ {true} | false ]
%    - diffuse_all : [ true | {false} ]
%
% Returns:
%    :
%    - P0 :
%    - Retcode :
%
% Note:
%
% Example:
%
%    See also:

if nargin<5
    diffuse_all=false;
    if nargin<4
        check_unit_roots=true;
        if nargin<3
            diffuse_factor=10;
        end
    end
end
endo_nbr=size(T,1);

P0=diag(diffuse_factor*ones(endo_nbr,1));
retcode=0;
if ~diffuse_all
    if check_unit_roots
        % check whether there are diffuse_all elements as they require special treatment
        [Ts,Us,stable]=schur_form(T);
        RQRs=Us'*RQR*Us;
    else
        Ts=T;
        Us=1;
        RQRs=RQR;
        stable=true(1,endo_nbr);
    end
    if any(stable)
        [tmp,retcode]=lyapunov_equation(Ts(stable,stable),RQRs(stable,stable));
        if ~retcode
            % this does not harm
            P0(stable,stable)=symmetrize(tmp);
            % now transform back to y units
            % y=Us*x
            P0=Us*P0*Us';
        end
    else
        retcode=1;
    end
end
% tmp=Sandwich_solve(Ts(stable,stable),Ts(stable,stable)',RQRs(stable,stable))

function [Ts,Us,stable]=schur_form(T)
[UU,TT]=schur(T);
E=ordeig(TT);
tol=1e-9;
unit=abs(abs(E)-1)<tol;
stable=abs(E)<1-tol;
if any(unit)
    clusters=zeros(size(E));
    clusters(unit)=1;
    [Us,Ts]=ordschur(UU,TT,clusters);
    stable=[false(1,sum(unit)),true(1,sum(stable))];
else
    Us=1;
    Ts=T;
end