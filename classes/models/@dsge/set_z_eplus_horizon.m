function obj=set_z_eplus_horizon(obj)
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

% get the shock horizon for all the shocks
%------------------------------------------
shock_horizon=obj.exogenous.shock_horizon;
kmax=max(obj.exogenous.shock_horizon(:));

% deterministic variables can also be anticipated...
if kmax>0
    solve_order=obj.options.solve_order;
    zzzz=repmat('z',1,solve_order);
    % locate the eplus variables that are 0
    %--------------------------------------
    e_0=obj.locations.after_solve.z.e_0(end);
    nz=size(obj.solution.Tz{1},2);
    exo_nbr=sum(obj.exogenous.number);
    number_of_regimes=obj.markov_chains.regimes_number;
    for ireg_=1:number_of_regimes
        do_one_regime()
    end
end

    function do_one_regime()
        this_regime_shock_horizon=shock_horizon(ireg_,:);
        % the current shocks are unrestricted
        %------------------------------------
        offset=e_0;
        % future shocks are restricted to their respective horizons
        %----------------------------------------------------------
        zproto=false(1,nz);
        for iplus=1:kmax
            bad=this_regime_shock_horizon<iplus;
            bad_locs=offset+find(bad);
            zproto(bad_locs)=true;
            offset=offset+exo_nbr;
        end
        zkz=zproto;
        for io=1:solve_order
            tz=['T',zzzz(1:io)];
            obj.solution.(tz){ireg_}(:,zkz)=0;
            if io<solve_order
                zkz=logical(kron(zkz,zproto));
            end
        end
    end
end
