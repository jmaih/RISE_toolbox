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

map_regime=double.empty;

if iscell(ref_state)
    
    map_regime=ref_state{2};
    
    ref_state=ref_state{1};
    
end

nstate_vars=numel(state_vars_location);

endo_nbr=size(y0.y,1);

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

% initialize with nans in order to see if/when simulation breaks down
sim1=nan(endo_nbr,span);

regimes=ref_state*ones(span,1);

next_state=ref_state;

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
        
        t=0;
        
        while t<span
            
            if retcode
                
                %retcode=701;
                
                return
                
            end
            
            t=t+1;
            
            if next_state==ref_state
                
                regimes(t)=ref_state;
                
                % 1-step forecast
                [sim1(:,t),next_state]=forecaster(y00,Tfunc_ref,ref_state,t);
                
                if next_state~=ref_state
                    % last step failed
                    t=t-1;
                    
                end
                
            else
                
                t0=t;
                
                t=t-1;
                
                keep_going=true;
                
                place_holder=[];
                
                the_path=zeros(1,0);
                                
                while keep_going
                    
                    t=t+1;
                    
                    the_path(end+1)=next_state; %#ok<AGROW>
                    
                    bigt=t-t0+1;
                    
                    regimes(t:t+bigt-1)=the_path;
                    
                    [Ht,Gt,kt]=load_solutions(the_path);
                    
                    Tfunc=@(t)[Ht(:,state_vars_location,t),kt(:,t),Gt(:,:,t)];
                    % bigt-step forecast
                    [sim1(:,t0:t),next_state]=forecaster(y00,Tfunc,the_path,t);
                    
                    if isempty(place_holder)
                        
                        place_holder={sim1(:,t0:t),[]};
                        
                    elseif isempty(place_holder{2})
                        
                        place_holder{2}=sim1(:,t0:t);
                        
                    else
                        
                        place_holder{1}=place_holder{2};
                        
                        place_holder{2}=sim1(:,t0:t);
                        
                    end
                    
                    keep_going=t+bigt<span && (next_state~=ref_state);
                    
                    if keep_going
                        
                        if next_state~=the_path(end)
                            % last step failed
                            t=t-1;
                            
                            % shorten the path since the step is to be redone
                            the_path(end)=[];
                            
                        end
                        
                    elseif next_state==ref_state
                        % last step failed
                        t=t-1;
                        
                        if t>=t0
                            
                            sim1(:,t0:t)=place_holder{1}(:,1:t-t0+1);
                            
                        end
                        
                    end
                    
                end
                
            end
            
            y00=sim1(:,max(1,t));
            
        end
        
    end

    function [sim1,r,is_viol]=forecaster(y00,Tfunc,the_path,t)
        
        nsteps=numel(the_path);
        
        sim1=y00(:,ones(1,nsteps));
                
        for istep=1:nsteps
            
            y00_s=y00-ss{1}; % steady state is the same???
            
            state=[y00_s(state_vars_location);1;shocks(:,t+istep-1)];
            
            sim1(:,istep)=ss{1}+Tfunc(istep)*state;
            
            y00=sim1(:,istep);
            
            if any(~isfinite(y00))||~isreal(y00)
                
                error('non-real or infinite values in simulation')
                
            end
            
            r0=the_path(istep);
            
            [is_viol,r]=check_viol(y00,r0);
            
%             if is_viol && istep<nsteps
%                 
%                 error('violation before end or road')
%                 
%             end
            
        end
        
        function [v,r]=check_viol(x,r0)
            
            test=sep_compl(x);
            
            disp(test)
            
            if numel(test)==1
                
                test=[test,-test];
                
            end
            
            if isempty(map_regime)
                
                if numel(test)~=2
                    
                    error('A function generalizing occbin is needed')
                    
                end
                
                v=test(r0)<0;
                
                if v
                    
                    r=setdiff(1:h,r0); 
                
                else
                    
                    r=r0; 
                
                end
                
            else
                
                [v,r]=map_regime(test,r0);
                                
            end
            
        end
        
    end

    function [Ht,Gt,kt]=load_solutions(the_path)
        
        bigt=numel(the_path);
        
        wingspan=bigt+1;
        
        Ht=H(:,:,ones(1,wingspan));
        
        Gt=G(:,:,ones(1,wingspan));
        
        kt=k(:,ones(1,wingspan));
        
        for t=bigt:-1:1
            
            st=the_path(t);
            
            if use_pinv
                
                AHAi=-pinv(Aplus(:,:,st)*Ht(:,:,t+1)+...
                    A0(:,:,st));
                
            else
                
                AHAi=-(Aplus(:,:,st)*Ht(:,:,t+1)+...
                    A0(:,:,st))\eye(endo_nbr);
                
            end
            
            Ht(:,:,t)=AHAi*Aminus(:,:,st);
            
            Gt(:,:,t)=AHAi*B{st};
            
            kt(:,t)=AHAi*(c(:,st)+Aplus(:,:,st)*kt(:,t+1));
            
        end
            
    end

end