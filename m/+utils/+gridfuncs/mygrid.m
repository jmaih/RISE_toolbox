function Regimes=mygrid(v)
% v is a vector of the number of states in each dimension
% Regimes is an prod(v) x NumberOfDimensions matrix describing the different
% composite regimes
% example
% Regimes=mygrid(10); : one chain with 10 regimes
% Regimes=mygrid([10,5]) : 2 chains with 5 regimes for the second one
% Regimes=mygrid([10,5,3]): 3 chains with 3 regimes on the third one
n=numel(v);
% first construct the indexes of the intervals
Regimes=[]; % first to v(1)th interval

for p=1:n
    Regimes=utils.gridfuncs.build_grid(Regimes,v(p));
end
end
