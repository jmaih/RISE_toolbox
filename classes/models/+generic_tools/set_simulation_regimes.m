function State=set_simulation_regimes(simul_regime,simul_periods,simul_burn)
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

if nargin<3
    simul_burn=[];
end
if isempty(simul_burn)
    simul_burn=0;
end

State=nan(simul_burn+simul_periods+1,1);
if ~isempty(simul_regime)
    simul_regime=simul_regime(:);
    if isscalar(simul_regime)
        State(1:end)=simul_regime;
    else
        nregs=numel(simul_regime);
        tmp=min(nregs,simul_periods);
        State(1:simul_burn)=simul_regime(1);
        State(simul_burn+(1:tmp))=simul_regime(1:tmp);
        State(simul_burn+tmp+1:end)=State(simul_burn+tmp);
    end
end
end
