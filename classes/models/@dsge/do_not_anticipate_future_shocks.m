function obj=do_not_anticipate_future_shocks(obj)
% do_not_anticipate_future_shocks - kills the effects of future shocks
%
% Syntax
% -------
% ::
%
%   obj=do_not_anticipate_future_shocks(obj)
%
% Inputs
% -------
%
% - **obj** [dsge|rise]: model object
%
% Outputs
% --------
%
% - **obj** [dsge|rise]: model object
%
% More About
% ------------
%
% - this function sets to 0 all the impacts of future shocks irrespective
%   of whether some have a longer horizon than others or not. All this is
%   triggered in the presence of future shocks, if option
%   **irf_anticipate** is set to false. 
%
% Examples
% ---------
%
% See also: 

kmax=max(obj.exogenous.shock_horizon(:));
% deterministic variables can also be anticipated...
if kmax>0 && ~obj.options.irf_anticipate
    solve_order=obj.options.solve_order;
    zzzz=repmat('z',1,solve_order);
    % locate the eplus variables that are 0
    %--------------------------------------
    e_0=obj.locations.after_solve.z.e_0(end);
    nz=size(obj.solution.Tz{1},2);
    zproto=false(1,nz);
    exo_nbr=sum(obj.exogenous.number);
    offset=e_0;
    bad=1:exo_nbr;
    for iplus=1:kmax
        bad_locs=offset+bad;
        zproto(bad_locs)=true;
        offset=offset+exo_nbr;
    end
    zkz=zproto;
    number_of_regimes=obj.markov_chains.regimes_number;
    for io=1:solve_order
        tz=['T',zzzz(1:io)];
        for ireg_=1:number_of_regimes
            obj.solution.(tz){ireg_}(:,zkz)=0;
            obj.solution.(tz){ireg_}=sparse(obj.solution.(tz){ireg_});
        end
        if io<solve_order
            zkz=logical(kron(zkz,zproto));
        end
    end
end

end
