function Regimes=mygrid(v)
% v is a vector of the number of states in each dimension
% Regimes is an prod(v) x NumberOfDimensions matrix describing the different
% composite regimes
% example
% Regimes=chain_grid(10); : one chain with 10 regimes
% Regimes=chain_grid([10,5]) : 2 chains with 5 regimes for the second one
% Regimes=chain_grid([10,5,3]): 3 chains with 3 regimes on the third one
n=numel(v);
% first construct the indexes of the intervals
Regimes=transpose(1:v(1)); % first to v(1)th interval

for p=2:n
    [rg,cg]=size(Regimes);
    vp=v(p);
    Ip=repmat(transpose(1:vp),1,rg);
    G0=nan(rg*vp,cg);
    for ii=1:rg
        G0((ii-1)*vp+1:ii*vp,:)=Regimes(ii*ones(vp,1),:);
    end
    Regimes=[G0,Ip(:)];
end
end