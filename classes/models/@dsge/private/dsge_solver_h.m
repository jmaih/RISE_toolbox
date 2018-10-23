function [T,eigval,retcode,obj,nsols]=dsge_solver_h(obj,sm)
% INTERNAL FUNCTION
%

% the obj going out probably contains the changed options

if isempty(obj)
    
    if nargout>1
        
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
        
    end
    
    T=struct();
    
    return
    
end

%% begin
T=struct();

% options.solve_order 1
%--------
if obj.options.solve_order<1
    
    return
    
end

[pos,siz,shock_horizon]=dsge_tools.rehash_topology(obj,sm);

% Structure of elements that will move across different orders
%-------------------------------------------------------------

[T.Tz,eigval,retcode,obj.options]=solve_first_order(sm,...
    siz,pos,obj.options,shock_horizon);

if all(retcode)
    
    nsols=0;
    
    return
    
end

nsols=size(T.Tz,3);

ss=sm.y_ordered;

g=sm.bgp_ordered;

order=obj.options.solve_order;

if obj.options.solve_order>1
    
    k=max(obj.exogenous.shock_horizon(:));
    
    % higher orders
    %--------------
    if isempty(obj.hops)||obj.hops.k~=k
        
        ns=sum(obj.endogenous.is_static);
        
        np=sum(obj.endogenous.is_predetermined);
        
        nb=sum(obj.endogenous.is_pred_frwrd_looking);
        
        nf=sum(obj.endogenous.is_frwrd_looking);
        
        ne0=obj.exogenous.number(1);
        
        h=obj.markov_chains.regimes_number;
        
        sig_in_v=siz.nsig==1;
        
        obj.hops=hops(ns,np,nb,nf,ne0,k,h,sig_in_v);
        
    end
    
    this0=obj.hops;
    
    % transfer derivatives
    %-----------------------
    
    f=struct();
    
    vorder=repmat('v',1,order);
    
    zorder=strrep(vorder,'v','z');
    
    for ii=1:order
        
        v=vorder(1:ii);
        
        f.(v)=sm.(['d',v]);
        
    end
    
    clear sm
    
    relev_fields=fieldnames(this0.options);
    
    opts=obj.options;
    
    for ii=1:numel(relev_fields)
        
        ff=relev_fields{ii};
        
        if isfield(opts,ff)
            
            this0.options.(ff)=opts.(ff);
            
        end
        
    end
    
    % compute moments for the shocks
    %-------------------------------
    M=moments_generator();
    
    % solve for higher orders
    %--------------------------    
    for isol=1:nsols
        
        if retcode(isol)
            
            continue
            
        end
        
        [this,retcode(isol)]=solve(this0,f,T.Tz(:,:,isol),ss,g,order,M);
        
        for io=1:order
            
            Tz_io=['T',zorder(1:io)];
            
            try
                
                T.(Tz_io)(:,:,isol)=this.(Tz_io);
                
            catch
                
                disp(['Order #',int2str(io),' failed because "',decipher(retcode(isol)),'"'])
                
            end
            
        end
        
    end
    
end

% solve for growth constant
%--------------------------
T=growth_component_solver(T,ss,g,order,pos,retcode);

    function M=moments_generator()
        
        M=[];
        
        homf=obj.options.solve_user_defined_shocks;
        
        if isempty(homf)
            
            return
            
        end
        
        engine=homf(obj,true,order);
        
        template_moments=cell2mat(hops.set_up_normal_moments(order,1));
        
        template_moments=template_moments(:).';
        
        nshocks=sum(obj.exogenous.number);
        
        template_moments=template_moments(ones(1,nshocks),:);
        
        xlist_user=fieldnames(engine);
        
        xlist=obj.exogenous.name;
        
        locs=locate_variables(xlist_user,xlist);
        
        M=template_moments;
        
        % push the user's defined moments
        %---------------------------------
        for ivar=1:numel(xlist_user)
            
            M(locs(ivar),:)=engine.(xlist_user{ivar})(:).';
            
        end
        
        if ~all(M(:,1)==0)
            
            error('all first-order moments should be zero')
            
        end
        
        if ~all(abs(M(:,2)-1)<1e-10)
            
            error('all second-order moments should be 1')
            
        end
        
        % Now construct the cell array for each order
        %--------------------------------------------
        M=hops.set_up_moments(order,nshocks,M);
        
        % we will need to find a way to deal with observed exogenous as
        % some point and in particular, finding a way to partial them out
        % so that the user doesn't have to be bothered by them.
        % Alternatively the user could just consider that the observed
        % exogenous have 0 higher-order moments...
        
        
        %         for o=1:order
        %
        %             ping=repmat({locs},1,o);
        %
        %             M{o}=M{o}(ping{:});
        %
        %         end
        
    end

end


function T=growth_component_solver(T,ss,g,solve_order,pos,retcode)
% INTERNAL FUNCTION
%

zzz=repmat('z',1,solve_order);

isstate=pos.z.pb;

nstate_vars=numel(isstate);

[~,h,nsols]=size(T.Tz);

nshocks=size(T.Tz{1,1,1},2)-(nstate_vars+1);

zero_shocks=zeros(nshocks,1);

sig_0=0;

for isol=1:nsols
    
    if retcode(isol)
        
        continue
        
    end
    
    for ireg=1:h
        
        do_one_regime();
        
    end
    
end

% now push final result back into T: write is as imaginary so as to
% 1- make it readily available if it is to be looked at
% 2- avoid nonlinear computations trying to mix it with Tsig
% 3- give the ability to turn it off separately
% The costs are as follows
% a- hiding it into Tz_sig implies readapting the routines calling
% one_step_engine
% b- how does it square with pruning? A pruned model will not prevent this
% from exploding...
%----------------------------------------------------------------------

    function do_one_regime()
        
        if ~any(g(:))
            
            return
            
        end
        
        x0=ss(:,ireg);
        
        x1=x0+g(:,ireg);
        
        x1_k=compute_forecast(x0);
        
        k=x1-x1_k;
        
        k(abs(k)<1e-10)=0;
        
        T.Tz{1,ireg,isol}(:,nstate_vars+1)=T.Tz{1,ireg,isol}(:,nstate_vars+1)+k;
        
        function y=compute_forecast(y0)
            
            y0t=y0-ss(:,ireg);
            
            zt=[y0t(isstate);sig_0;zero_shocks];
            
            y=ss(:,ireg);
            
            zkron=zt;
            
            for io=1:solve_order
                
                y=y+T.(['T',zzz(1:io)]){1,ireg,isol}*zkron;
                
                if io<solve_order
                    
                    zkron=kron(zkron,zt);
                    
                end
                
            end
            
        end
        
    end

end