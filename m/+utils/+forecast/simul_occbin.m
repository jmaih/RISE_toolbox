function [sim1,regimes,retcode]=simul_occbin(y0,T,ss,state_vars_location,...
options,shocks)

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

% evaluate stretch t:t+n
% evaluate t+n+1 with reference regime
% if no violations continue to t+n+2 from t+n+1
% if violations: extend the batch with the new conditions i.e. the
% violating regimes

use_pinv=false;

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

if iscell(ref_state)
    
    restr_map=ref_state{2};
    
    my_regimes=ref_state{3}(:,2:end);
    
    ref_state=ref_state{1};
    
    fields=fieldnames(restr_map);
    
    mycols=locate_variables(fields,my_regimes(1,:));
    
    my_regimes=cell2mat(my_regimes(2:end,mycols));
    
    ref_regime=my_regimes(ref_state,:);
    
else
    
    ref_regime=ref_state;
    
    my_regimes=(1:2).';
    
end

nstate_vars=numel(state_vars_location);

endo_nbr=size(y0.y,1);

span=options.nsteps+options.burn;

if isempty(options.sep_compl)
    
    sep_compl=@(x)1;
    
else
    
    sep_compl=@(x)options.sep_compl(x);
    
end

% find the solution for the reference regime
%--------------------------------------------
[H,G,k]=reference_solution();

% reference solution
Tref=@(varagin)[H(:,state_vars_location),k,G];

retcode=0;

y0.econd.data=shocks(:,:,ones(3,1));

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
        
        y00=y0.y;
                
        next_state=ref_state;
        
        t=0;
        
        while t<span
            
            if retcode
                
                return
                
            end
            
            t=t+1;
            
            [y1,regs1]=simulate_forecast(y00,shocks(:,t));
            
            rest_t=span-t+1;
            
            ny1=min(size(y1,2),rest_t);
            
            y1=y1(:,1:ny1);
            
            regs1=regs1(1:ny1);
            
            future_shocks=shocks(:,t+1:end);
            
            if isempty(future_shocks)||all(vec(future_shocks(:))==0)
                
                sim1(:,t:t+ny1-1)=y1;
                
                regimes(t:t+ny1-1)=regs1;
                
                t=t+ny1-1;
                
            else
                
                sim1(:,t)=y1(:,1);
                
                regimes(t)=regs1(1);
                
            end
            
            % initial conditions for next step
            y00=sim1(:,t);
            
        end
                
        function [y1,the_path]=simulate_forecast(y0,shkt)
            
            % steady state is the same irrespective of the regime followed
            y00_s=y0-ss{1};
            
            state=[y00_s(state_vars_location);1;shkt];
            
            y1=ss{1}+Tref()*state;
            
            next_state=map_regime(y1);
            
            if next_state==ref_state
                
                the_path=ref_state;
                
                return
                
            end
            
            good=false;
            
            nsteps=0;
            
            the_path=zeros(1,0);
            
            while ~good
                
                nsteps=nsteps+1;
                
                the_path(nsteps)=next_state;
                                
                [Ht,Gt,kt]=load_solutions(the_path);
                
                Tfunc=@(t)[Ht(:,state_vars_location,t),kt(:,t),Gt(:,:,t)];
                
                y01=y0(:,ones(1,nsteps));
                
                new_y0=y0;
                
                shkti=shkt;
                
                for istep=1:nsteps
                    
                    y00_s=new_y0-ss{1};
                    
                    state=[y00_s(state_vars_location);1;shkti];
                    
                    y01(:,istep)=ss{1}+Tfunc(istep)*state;
                    
                    if istep==1
                        % change the shocks going forward
                        shkti=0*shkt;
                        
                    end
                    
                    % next step
                    new_y0=y01(:,istep);
                    
                end
                
                % do one clean step to check we are back to normal times
                y00_s=new_y0-ss{1};
                
                state=[y00_s(state_vars_location);1;shkti];
                
                ytest=ss{1}+Tref()*state;
                
                next_state=map_regime(ytest);
                
                good=next_state==ref_state;
                
            end
            
            % format new output
            %------------------
            y1=y01;
            
        end
        
    end

    function r=map_regime(x)
        
        test=sep_compl(x);
        
        test=test(:,1);
        
        this_regime=ref_regime;
        
        r=ref_state;
        
        for ii=1:size(my_regimes,2)
            
            if test(ii)<0
                
                this_regime(ii)=setdiff([1,2],this_regime(ii));
                
            end
            
        end
        
        for ii=1:size(my_regimes,1)
            
            if all(my_regimes(ii,:)-this_regime==0)
                
                r=ii;
                
                break
                
            end
            
        end
        
    end

    function [Ht,Gt,kt]=load_solutions(the_path)
        
        bigt=numel(the_path);
        
        wingspan=bigt+1;
        
        Ht=H(:,:,ones(1,wingspan));
        
        Gt=G(:,:,ones(1,wingspan));
        
        kt=k(:,ones(1,wingspan));
        
        if all(the_path==ref_state)
            % quick exit if we are in the reference state
            return
            
        end
        
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