function R=build_grid(R,vp)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 


% see also mygrid

[rg,cg]=size(R);
rg_star=max(rg,1);

Ip=repmat(transpose(1:vp),1,rg_star);
G0=nan(rg_star*vp,cg);
for ii=1:rg
    G0((ii-1)*vp+1:ii*vp,:)=R(ii*ones(vp,1),:);
end
R=[G0,Ip(:)];
end