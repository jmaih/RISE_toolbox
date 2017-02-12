function [sim1,regimes,retcode]=simul_occbin(y0,T,ss,state_vars_location,...
    options,shocks)

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

use_pinv=true;
is_accelerated=true;

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

endo_nbr=size(y0.y,1);
exo_nbr=size(B{1},2);

span=options.nsteps+options.burn;

compl=@(x)isempty(options.complementarity)||options.complementarity(x);

% find the solution for the reference regime
%--------------------------------------------
[H,G,k]=reference_solution();

retcode=0;
y0.econd.data=shocks(:,:,ones(3,1));

last_step=0;
bigt=0;
Ht0=[];
Gt0=[];
kt0=[];

% initialize with nans in order to see if/when simulation breaks down
sim1=nan(endo_nbr,span);
regimes=ref_state*ones(span,1);
do_one_occbin_path()

    function [H,G,k]=reference_solution()
        H=zeros(endo_nbr);
        H(:,state_vars_location)=T{ref_state}(:,1:nstate_vars);
        k=full(T{ref_state}(:,nstate_vars+1));
        G=full(T{ref_state}(:,nstate_vars+2:end));
    end

    function do_one_occbin_path()
        % try and forecast from the ref
        %------------------------------
        Tfunc=@(t)T{ref_state};
        is_viol=forecaster(y0.y,Tfunc);
        if ~is_viol
            return
        end
        bigt=1;
        last_step=1;
        % accelerate things by preparing the longest possible stretch
        %-------------------------------------------------------------
        if is_accelerated
            [Ht0,Gt0,kt0]=resolve_occbin(span-1);
        end
        
        % make assumption about the duration of the violation and retry
        %---------------------------------------------------------------
        if last_step<=1
            y00=y0.y;
        else
            y00=sim1(:,last_step-1);
        end
        while is_viol && last_step-1+bigt<span
            [Ht,Gt,kt]=load_solutions(bigt);
            if retcode
                % we exit as quickly as possible, it is not possible to
                % tame the system.
                bigt=inf;
            else
                Tfunc=@(t)[Ht(:,state_vars_location,t),kt(:,t),Gt(:,:,t)];
                is_viol=forecaster(y00,Tfunc,last_step);
                if is_viol
                    bigt=bigt+1;
                    if last_step<=1
                        y00=y0.y;
                    else
                        y00=sim1(:,last_step-1);
                    end
                end
            end
        end
        if is_viol
            retcode=701;
        end
    end

    function [Ht,Gt,kt]=load_solutions(bigt)
        if is_accelerated
            last_step_=max(1,last_step);
            wing=last_step_:span;
            nwing=numel(wing);
            stretch=span-bigt:span-1;
            Ht=H(:,:,ones(nwing,1));Ht(:,:,1:bigt)=Ht0(:,:,stretch);
            Gt=G(:,:,ones(nwing,1));Gt(:,:,1:bigt)=Gt0(:,:,stretch);
            kt=k(:,ones(nwing,1));kt(:,1:bigt)=kt0(:,stretch);
            regimes(last_step_:last_step_+bigt-1)=other_state;
            regimes(last_step_+bigt:end)=ref_state;
        else
            [Ht,Gt,kt]=resolve_occbin(bigt);
        end
    end

    function is_viol=forecaster(y00,Tfunc,istart)
        if nargin<3
            istart=1;
        end
        iter=0;
        for istep=istart:span
            iter=iter+1;
            y00_s=y00-ss{1}; % steady state is the same???
            state=[y00_s(state_vars_location);1;shocks(:,istep)];
            sim1(:,istep)=ss{1}+Tfunc(istep-istart+1)*state;
            y00=sim1(:,istep);
            if any(~isfinite(y00))||~isreal(y00)
                retcode=701;
            end
            is_viol=~compl(y00);
            if is_viol||retcode
                break
            end
        end
        % if the period after the guess does not violate, then the guess
        % was correct. Update the last step
        %----------------------------------------------------------------
        if ~is_viol||(iter>bigt && compl(sim1(:,istart+bigt)))
            last_step=istep;
            bigt=0;
        end
    end

    function [Ht,Gt,kt]=resolve_occbin(bigt)
        new_span=span-last_step+1;
        Ht=zeros(endo_nbr,endo_nbr,new_span);
        Gt=zeros(endo_nbr,exo_nbr,new_span);
        kt=zeros(endo_nbr,new_span);
        
        t=bigt+1:new_span;
        Ht(:,:,t)=H(:,:,ones(1,new_span-bigt));
        Gt(:,:,t)=G(:,:,ones(1,new_span-bigt));
        kt(:,t)=k(:,ones(1,new_span-bigt));
        regimes(t+last_step-1)=ref_state;
        
        maxed_out=bigt>=new_span;
        if maxed_out
            retcode=701;
        else
            for t=bigt:-1:1
                if use_pinv
                    AHAi=-pinv(Aplus(:,:,other_state)*Ht(:,:,t+1)+...
                        A0(:,:,other_state));
                else
                    AHAi=-(Aplus(:,:,other_state)*Ht(:,:,t+1)+...
                        A0(:,:,other_state))\eye(endo_nbr);
                end
                Ht(:,:,t)=AHAi*Aminus(:,:,other_state);
                Gt(:,:,t)=AHAi*B{other_state};
                kt(:,t)=AHAi*(c(:,other_state)+Aplus(:,:,other_state)*kt(:,t+1));
                regimes(t+last_step-1)=other_state;
            end
        end
    end
end

