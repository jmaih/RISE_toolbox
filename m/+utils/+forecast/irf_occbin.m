function [irfs,retcode]=irf_occbin(y0,T,ss,state_vars_location,...
    which_shocks,det_vars,Initcond,options)
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

% find the first-order approximation of the system
%--------------------------------------------------

A0=options.occbin.A0;
Aplus=A0;
h=size(A0,3);
for ireg=1:h
    Aplus(:,:,ireg)=options.occbin.Gplus01(:,:,ireg,ireg);
end
Aminus=options.occbin.Aminus;
B=options.occbin.B;
c=options.occbin.c;

ref_state=options.solve_occbin;
other_state=setdiff(1:h,ref_state);
if numel(other_state)>1
    error('occbin can only handle one constraint')
end
nstate_vars=numel(state_vars_location);

[endo_nbr,nlags]=size(y0.y);

nsteps=Initcond.nsteps;

compl=@(x)isempty(Initcond.complementarity)||Initcond.complementarity(x);

% find the solution for the reference regime
%--------------------------------------------
[H,G,k]=reference_solution();

exo_nbr=numel(which_shocks);

nshocks=sum(which_shocks);

if any(which_shocks & det_vars)
    error('trying to compute irfs for deterministic shocks')
end

irfs=zeros(endo_nbr,nlags+Initcond.nsteps,nshocks,Initcond.nsimul);

iter=0;

retcode=0;
orig_shocks=y0.econd.data(:,:,1);
for ishock=1:exo_nbr
    if ~which_shocks(ishock)
        continue
    end
    iter=iter+1;
    shock_id=ishock;
    for isimul=1:Initcond.nsimul
        if ~retcode
            shocks=utils.forecast.create_shocks(exo_nbr,shock_id,det_vars,Initcond);
            % make sure that impulses that are inherited stay on in the
            % reference simulation. This may imply over-riding the impulse
            % itself if it happens to be on the path of an inherited shock
            %--------------------------------------------------------------
            inherited_shocks=orig_shocks~=0;
            shocks(inherited_shocks)=orig_shocks(inherited_shocks);
            y0.econd.data=shocks(:,:,ones(3,1));
            if ~retcode
                sim1=do_one_occbin_path();
                if ~retcode
                    path1=[y0.y,sim1];
                    irfs(:,:,iter,isimul)=path1;
                end
            end
        end
    end
end

    function [H,G,k]=reference_solution()
        H=zeros(endo_nbr);
        H(:,state_vars_location)=T{ref_state}(:,1:nstate_vars);
        k=T{ref_state}(:,nstate_vars+1);
        G=T{ref_state}(:,nstate_vars+2:end);
    end

    function sim1=do_one_occbin_path()
        % try and forecast from the ref
        %------------------------------
        Tfunc=@(t)T{ref_state};
        [sim1,is_viol]=forecaster(Tfunc);
        if ~is_viol
            return
        end
        
        % make assumption about the duration of the violation and retry
        %---------------------------------------------------------------
        bigt=1;
        while is_viol
            [Ht,Gt,kt]=resolve_occbin(bigt);
            Tfunc=@(t)[Ht(:,state_vars_location,t),kt(:,t),Gt(:,:,t)];
            [sim1,is_viol]=forecaster(Tfunc);
            bigt=bigt+1;
        end        
    end

    function [sim1,is_viol]=forecaster(Tfunc)
        sim1=zeros(endo_nbr,Initcond.nsteps);
        y00=y0.y;
        for istep=1:nsteps
            y00_s=y00-ss{1}; % steady state is the same???
            state=[y00_s(state_vars_location);1;shocks(:,istep)];
            sim1(:,istep)=ss{1}+Tfunc(istep)*state;
            y00=sim1(:,istep);
            is_viol=~compl(y00);
            if is_viol
                break
            end
        end
    end

    function [Ht,Gt,kt]=resolve_occbin(bigt)
        Ht=zeros(endo_nbr,endo_nbr,nsteps);
        Gt=zeros(endo_nbr,exo_nbr,nsteps);
        kt=zeros(endo_nbr,nsteps);
        for t=bigt+1:nsteps
            Ht(:,:,t)=H;
            Gt(:,:,t)=G;
            kt(:,t)=k;
        end
        maxed_out=bigt==nsteps;
        if maxed_out
            error('maxed out: maybe increasing the horizon may help...')
        end
        for t=bigt:-1:1
            AHAi=-(Aplus(:,:,other_state)*Ht(:,:,t+1)+A0(:,:,other_state))\eye(endo_nbr);
            Ht(:,:,t)=AHAi*Aminus(:,:,other_state);
            Gt(:,:,t)=AHAi*B(:,:,other_state);
            kt(:,t)=AHAi*(c(:,other_state)+Aplus(:,:,other_state)*kt(:,t+1));
        end
    end
end