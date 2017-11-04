function [Histdec,obj]=hd(obj,varargin)
% historical_decomposition Computes historical decompositions of a DSGE model
%
% Syntax
% -------
% ::
%
%   [Histdec,obj]=hd(obj)
%   [Histdec,obj]=hd(obj,varargin)
%
% Inputs
% -------
%
% - obj : [rise|dsge|rfvar|svar] model(s) for which to compute the
%   decomposition. obj could be a vector of models
%
% - varargin : standard optional inputs **coming in pairs**. Among which:
%   - **histdec_start_date** : [char|numeric|{''}] : date at which the
%     decomposition starts. If empty, the decomposition starts at he
%     beginning of the history of the dataset
%
% Outputs
% --------
%
% - Histdec : [struct|cell array] structure or cell array of structures
%   with the decompositions in each model. The decompositions are given in
%   terms of:
%   - the exogenous variables
%   - **InitialConditions** : the effect of initial conditions
%   - **sig** : measure of the effect of non-certainty equivalence
%   - **switch** : the effect of switching (which is also a shock!!!)
%   - **steady_state** : the contribution of the steady state
%
% Remarks
% --------
%
% - the elements that do not contribute to any of the variables are
%   automatically discarded.
%
% - **N.B** : a switching model is inherently nonlinear and so, strictly
%   speaking, the type of decomposition we do for linear/linearized
%   constant-parameter models is not feasible. RISE takes an approximation
%   in which the variables, shocks and states matrices across states are
%   averaged. The averaging weights are the smoothed probabilities.
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        Histdec=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

nobj=numel(obj);

if nobj>1
    
    Histdec=cell(1,nobj);
    
    for iobj=1:nobj
        
        [Histdec{iobj},obj(iobj)]=hd(obj(iobj),varargin{:});
        
    end
    
    return
    
end

obj=set(obj,varargin{:});

[obj,~,~,retcode]=filter(obj);

if retcode
    
    error(decipher(retcode))
    
end

init=filter_initialization(obj);

ov=obj.order_var;

order_var_endo_names=obj.endogenous.name(ov);

shock_names=obj.exogenous.name;

[Ty,Te,Tsig]=rebuild_state_space();

ss=init.steady_state;

histdec_start_date=obj.options.histdec_start_date;

first=[];

[smoothed_vars,dn]=load_values(obj.filtering.smoothed_variables,...
    order_var_endo_names);

[smoothed_shocks]=load_values(obj.filtering.smoothed_shocks,shock_names);

smoothed_probabilities=double(...
    ts.collect(obj.filtering.smoothed_regime_probabilities));

smoothed_probabilities=permute(smoothed_probabilities,[2,1]);

select_range()

smoothed_vars(init.is_log_var_new_order,:,:)=log(smoothed_vars(...
    init.is_log_var_new_order,:,:));

% enforce 1
%-----------
nrounds=10;
for iround=1:nrounds
    
    smoothed_probabilities=bsxfun(@rdivide,smoothed_probabilities,...
        sum(smoothed_probabilities,1));
    
%     max(abs(sum(smoothed_probabilities,1)-1))
    
end

[D,~,switch_id,sig_id,trend_id]=decomposition_engine(smoothed_vars,...
    smoothed_shocks,Ty,Te,Tsig,ss);

contrib_names=[{'y0','switch','sig','trend'},shock_names];

discard=false(1,size(D,3));

for ii=[switch_id,sig_id,trend_id]
    
    di=D(:,:,ii,:);
    
    discard(ii)=all(abs(di(:))<1e-9);
    
end

contrib_names=contrib_names(~discard);

D=D(:,:,~discard,:);

D=sum(bsxfun(@times,D,permute(smoothed_probabilities,[3,2,4,1])),4);

D=permute(D,[2,3,1]);

Histdec=struct();

for v=1:numel(order_var_endo_names)
    
    if v==1
        
        proto=ts(histdec_start_date,D(:,:,v),contrib_names);
        
    else
        
        proto=reset_data(proto,D(:,:,v),contrib_names);
        
    end
    
    Histdec.(order_var_endo_names{v})=proto;
    
end

    function select_range()
        
        hist_start_date=dn(1);
        
        hist_end_date=dn(end);
        
        if isempty(histdec_start_date)
            
            histdec_start_date=hist_start_date;
            
        else
            
            histdec_start_date=date2serial(histdec_start_date);
            
        end
        
        if histdec_start_date<hist_start_date || ...
                histdec_start_date>hist_end_date
            
            error([mfilename,':: the decomposition start date must lie between ',...
                serial2date(hist_start_date),' and ',serial2date(hist_end_date)])
            
        end
        
        first=find(histdec_start_date==(hist_start_date:hist_end_date));
        
        
        smoothed_vars=smoothed_vars(:,first:end,:);
        
        smoothed_shocks=smoothed_shocks(:,first:end,:);
        
        smoothed_probabilities=smoothed_probabilities(:,first:end);
        
    end


    function [Ty,Te,Tsig]=rebuild_state_space()
        
        Ty=init.Tx;
        
        ny=size(Ty{1},1);
        
        tmp=zeros(ny);
        
        Te=init.Te;
        
        for ireg=1:numel(Ty)
            
            tmp(:,init.state_vars_location)=Ty{ireg};
            
            Ty{ireg}=tmp;
            
            % expand for anticipation terms
            %------------------------------
            Te{ireg}=Te{ireg}(:,:);
            
        end
        
        Tsig=init.Tsig;
        
    end


    function [v,dn,nt,np]=load_values(smooth,vnames)
        
        nt=smooth.(vnames{1}).NumberOfObservations;
        
        np=smooth.(vnames{1}).NumberOfVariables;
        
        nv=numel(vnames);
        
        v0=nan(nv,nt,np);
        
        for iv=1:nv
            
            if iv==1
                
                dn=smooth.(vnames{iv}).date_numbers;
                
            end
            
            v0(iv,:,:)=permute(double(smooth.(vnames{iv})),[3,1,2]);
            
        end
        
        v=v0;
        
%         v=v0(:,:,1);
%         
%         for ireg=2:np
%             
%             v=[v;v0(:,:,ireg)]; %#ok<AGROW>
%             
%         end
        
    end

end

function [D,y0_id,switch_id,sig_id,trend_id,shocks_id]=decomposition_engine(...
    smoothvars,shocks,Ty,Te,Tsig,ss)

ny=size(Ty{1},1);

[nx,nt,np]=size(shocks);

% shock columns and rows may have been added in the state matrices during
% filtering. Get rid of them
retrim()

nz=1+1+1+1+nx;

y0_id=1;

switch_id=y0_id+1;

sig_id=switch_id+1;

trend_id=sig_id+1;

shocks_id=trend_id+(1:nx);

nonswitch=setdiff(1:nz,switch_id);

D=zeros(ny,nt,nz,np);

for st=1:np
    
    Yst=bsxfun(@minus,smoothvars(:,:,st),ss{st});
    
    y0=Yst(:,1);
    
    for t=1:nt
        
        if t==1
            
            D(:,1,y0_id,st)=y0;
            
            continue
            
        end
        
        % Initial conditions
        %-------------------
        y0=D(:,t-1,y0_id,st);
        
        D(:,t,y0_id,st)=Ty{st}*y0;
        
        % risk/constant term
        %--------------------
        Csig0=D(:,t-1,sig_id,st);
        
        D(:,t,sig_id,st)=Ty{st}*Csig0+get_Csig();
        
        % trend
        %-------
        Ctrend0=D(:,t-1,trend_id,st);
        
        D(:,t,trend_id,st)=Ty{st}*Ctrend0+get_Ctrend();
        
        % shocks
        %-------
        Sh0=D(:,t-1,shocks_id,st);
        
        D(:,t,shocks_id,st)=Ty{st}*permute(Sh0,[1,3,2])+get_shocks();
        
        % discrepancy/switch
        %--------------------
        Css0=sum(D(:,t,nonswitch,st),3);
        
        D(:,t,switch_id,st)=Yst(:,t)-Css0;
        
    end
    
end

    function retrim()
        
        my=size(smoothvars,1);
        
        if my==ny
            
            return
            
        end
        
        % remove the shocks
        ny=my;
        
        xrange=1:my;
        
        for ireg=1:np
             Ty{ireg}=Ty{ireg}(xrange,xrange);
             
             Te{ireg}=Te{ireg}(xrange,:);
             
             Tsig{ireg}=Tsig{ireg}(xrange,:);
             
             ss{ireg}=ss{ireg}(xrange,:);
             
        end
        
    end

    function sh=get_shocks()
        
        sh=bsxfun(@times,Te{st},shocks(:,t,st).');
        
    end

    function Ctrend=get_Ctrend()
        
        Ctrend=imag(Tsig{st});
        
    end

    function Csig=get_Csig()
        
        Csig=real(Tsig{st});
        
    end

end


function d=the_defaults()

d={
    'histdec_start_date','',@(x)is_date(x)||is_serial(x),...
    'histdec_start_date must be a valid date'
    };

end
