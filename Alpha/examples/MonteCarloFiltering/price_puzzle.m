function flag=price_puzzle(obj)

[obj,retcode]=solve(obj);
if retcode
    flag=false;
else
    % inflation
    vloc=strcmp('PIE',obj.endogenous.name);
    % monetary policy shock
    sloc=strcmp('EI',obj.exogenous.name);
    response=obj.solution.m_e{1}(vloc,sloc);
    flag=response>0;
end

