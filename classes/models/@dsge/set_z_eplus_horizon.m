function obj=set_z_eplus_horizon(obj)
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

shock_horizon=obj.exogenous.shock_horizon;
kmax=max(obj.exogenous.shock_horizon);
% deterministic variables can also be anticipated...
if kmax>0
    solve_order=obj.options.solve_order;
    zzzz=repmat('z',1,solve_order);
    % locate the eplus variables that are 0
    %--------------------------------------
    e_0=obj.locations.after_solve.z.e_0(end);
    nz=size(obj.solution.Tz{1},2);
    zproto=false(1,nz);
    exo_nbr=sum(obj.exogenous.number);
    offset=e_0;
    for iplus=1:kmax
        bad=shock_horizon<iplus;
        bad_locs=offset+find(bad);
        zproto(bad_locs)=true;
        offset=offset+exo_nbr;
    end
    zkz=zproto;
    number_of_regimes=obj.markov_chains.regimes_number;
    for io=1:solve_order
        tz=['T',zzzz(1:io)];
        for ireg_=1:number_of_regimes
            obj.solution.(tz){ireg_}(:,zkz)=0;
        end
        if io<solve_order
            zkz=logical(kron(zkz,zproto));
        end
    end
end
end
