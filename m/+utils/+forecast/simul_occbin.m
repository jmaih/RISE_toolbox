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
% exo_nbr=size(B{1},2);

span=options.nsteps+options.burn;

if isempty(options.sep_compl)
    
    sep_compl=@(x)ones(2,1);
    
else
    
    sep_compl=@(x)options.sep_compl(x);
    
end

% find the solution for the reference regime
%--------------------------------------------
[H,G,k]=reference_solution();

retcode=0;

y0.econd.data=shocks(:,:,ones(3,1));

bigt=0; Ht0=[]; Gt0=[]; kt0=[];

% initialize with nans in order to see if/when simulation breaks down
sim1=nan(endo_nbr,span);

regimes=ref_state*ones(span,1);

reference=true;

Tfunc_ref=@(t)T{ref_state};

do_one_occbin_path()

    function [H,G,k]=reference_solution()
        
        H=zeros(endo_nbr);
        
        H(:,state_vars_location)=T{ref_state}(:,1:nstate_vars);
        
        k=full(T{ref_state}(:,nstate_vars+1));
        
        G=full(T{ref_state}(:,nstate_vars+2:end));
        
    end

    function do_one_occbin_path()
        
        y00=y0.y;
        
        % accelerate things by preparing the longest possible stretch
        %-------------------------------------------------------------
        [Ht0,Gt0,kt0]=resolve_occbin(span-1);
        
        t=0;
        
        while t<span
            
            if retcode
                
                %retcode=701;
                
                return
                
            end
            
            t=t+1;
            
            if reference
                
                regimes(t)=ref_state;
                
                % 1-step forecast
                [sim1(:,t),viol_last]=forecaster(y00,Tfunc_ref,1,t);
                
                if viol_last
                    
                    t=t-1;
                    
                end
                
            else
                
                t0=t;
                
                t=t-1;
                
                keep_going=true;
                
                place_holder=[];
                
                while keep_going
                    
                    t=t+1;
                    
                    regimes(t)=other_state;
                    
                    bigt=t-t0+1;
                    
                    [Ht,Gt,kt]=load_solutions(bigt);
                    
                    Tfunc=@(t)[Ht(:,state_vars_location,t),kt(:,t),Gt(:,:,t)];
                    % bigt-step forecast
                    [sim1(:,t0:t),viol_last]=forecaster(y00,Tfunc,bigt,t);
                    
                    if isempty(place_holder)
                        
                        place_holder={sim1(:,t0:t),[]};
                        
                    elseif isempty(place_holder{2})
                        
                        place_holder{2}=sim1(:,t0:t);
                        
                    else
                        
                        place_holder{1}=place_holder{2};
                        
                        place_holder{2}=sim1(:,t0:t);
                        
                    end
                    
                    keep_going=t+bigt<span && ~viol_last;
                    
                end
                
                if viol_last
                    
                    t=t-1;
                    
                    if t>=t0
                        
                        sim1(:,t0:t)=place_holder{1};
                        
                    end
                    
                end
                
            end
            
            y00=sim1(:,max(1,t));
            
            if viol_last
                
                reference=~reference;
                
            end
            
        end
        
    end

    function [sim1,is_viol]=forecaster(y00,Tfunc,nsteps,t)
        
        sim1=y00(:,ones(1,nsteps));
        
        is_viol=false;
        
        for istep=1:nsteps
            
            y00_s=y00-ss{1}; % steady state is the same???
            
            state=[y00_s(state_vars_location);1;shocks(:,t+istep-1)];
            
            sim1(:,istep)=ss{1}+Tfunc(istep)*state;
            
            y00=sim1(:,istep);
            
            if any(~isfinite(y00))||~isreal(y00)
                
                error('non-real or infinite values in simulation')
                
            end
            
            if check_viol(y00)
                
                is_viol=true;
                
                break
                
            end
            
        end
        
        function v=check_viol(x)
            
            test=sep_compl(x);
            
            if numel(test)==1
                
                test=[test,-test];
                
            end
            
            if reference
                
                test=test(1);
                
            else
                
                test=test(2);
                
            end
            
            v=test<0;
            
        end
        
    end

    function [Ht,Gt,kt]=load_solutions(bigt)
                                
        stretch=span-bigt:span-1;
        
        Ht=Ht0(:,:,stretch);
        
        Gt=Gt0(:,:,stretch);
        
        kt=kt0(:,stretch);
        
    end

    function [Ht,Gt,kt]=resolve_occbin(bigt)
        
        wingspan=bigt+1;
        
        Ht=H(:,:,ones(1,wingspan));
        
        Gt=G(:,:,ones(1,wingspan));
        
        kt=k(:,ones(1,wingspan));
        
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
            
        end
            
    end

end