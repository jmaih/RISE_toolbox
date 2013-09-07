function flag=price_puzzle(obj)

retcode=0;
if isempty(obj.T)
    [obj,retcode]=solve(obj);
end
if retcode
    flag=false;
else
    % inflation
    vloc=strcmp('PIE',{obj.varendo.name});
    % monetary policy shock
    sloc=strcmp('EI',{obj.varexo.name});
    response=obj.R(vloc,sloc);
    flag=response>0;
end

