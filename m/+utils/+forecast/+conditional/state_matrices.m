function [Tt,Rt,bt,Qt,Record]=state_matrices(T,R,mut,OMG,DPHI,DT,Record,ExpandedFlag)
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

if nargin<8
    ExpandedFlag=0;
end
ShockSpan=size(DPHI,2);

% indicate that an element is a covariance matrix by putting it into a cell
[mut,DPHI,DT,OMG]=utils.forecast.conditional.remove_holes(mut(:),DPHI,DT,{OMG});

if ~isempty(Record) && isequal(Record{1},DPHI)
	M1=Record{2};
else
	[M1,M2,RM2i]=utils.forecast.conditional.null_and_column_spaces(DPHI);
	Record={DPHI,M1,M2*RM2i};
end

MDPHIM=Record{3};
[endo_nbr,exo_nbr,NumberOfAnticipatedPeriods]=size(R);
ex_nbrstar=exo_nbr*NumberOfAnticipatedPeriods;
R=reshape(R,endo_nbr,ex_nbrstar);
% selector of the shocks
S=zeros(ex_nbrstar,ShockSpan);
S(:,1:ex_nbrstar)=eye(ex_nbrstar);

RSMDPHIM=R*S*MDPHIM;
% adjusted state matrices
if isempty(OMG)
    OMG=DPHI*transpose(DPHI);
end
Qt=MDPHIM*OMG*transpose(MDPHIM)+M1*transpose(M1);
wt=transpose(chol(Qt));
if ExpandedFlag
    Tt=[T-RSMDPHIM*DT,zeros(endo_nbr,ShockSpan)
        -MDPHIM*DT,zeros(ShockSpan)];
    bt=[RSMDPHIM
        MDPHIM]*mut;
    Rt=[R*S*wt
        wt];
else
    Tt=T-RSMDPHIM*DT;
    bt=RSMDPHIM*mut;
    Rt=R*S*wt;
end
