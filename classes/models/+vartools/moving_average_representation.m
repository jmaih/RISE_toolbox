function [MA,max_eigen]=moving_average_representation(B,L0,nlags,endo_nbr,nsteps)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% MA is of size n x n x (nsteps+1) the first element is the immediate
% impact at time t+1 and the last one is the impact at time t+nsteps

nz=size(B,2)-nlags*endo_nbr;
T=vartools.companion_form(B,nlags,endo_nbr,nz);
C0=[L0
    zeros(endo_nbr*(nlags-1),endo_nbr)];
MA=nan(endo_nbr,endo_nbr,nsteps);
Ti=C0;
for isteps=1:nsteps
    MA(:,:,isteps)=Ti(1:endo_nbr,:);
    Ti=T*Ti;
end
max_eigen=max(abs(eig(T)));
end
