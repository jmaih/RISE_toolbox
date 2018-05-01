function [welf,retcode,V,d]=loss(obj,simuls,varargin)

% LOSS calculates welfare
%
% ::
%
%   [welf,retcode,V,d]=loss(obj)
%   [welf,retcode,V,d]=loss(obj,simuls)
%   [welf,retcode,V,d]=loss(obj,simuls,varargin)
%
% Args:
%    - **obj** [rise|dsge]: scalar or vector of model objects.
%
%    - **simuls** [ts|struct|empty]: if empty, the unconditional welfare is
%      returned. if not empty the conditional welfare is calculated instead
%
%    - **varargin** []: optional arguments coming in pairs
%
% Returns:
%    :
%
%    - **welf** [scalar|vector]: conditional loss (scalar) or unconditional
%    loss (scalar if number of regimes = 1)
%
%    - **retcode** [scalar]: return code
%
%    - **V** [n x n x h]: array of equilibrium matrices
%
%    - **d** [h x 1 vector]: constant in loss
%
%
% Note:
%
%    The loss is such that v(x_t,r_t)=x_t.'*Vrt*x_t+drt
%
% Example:
%
%    See also:

% Note W,Tx are stored in alphabetical order for the convenience of the user
%------------------------------------------------------------------------

fslv=~true; % just for debugging purposes one can use fsolve instead of tfqmr

if nargin <2
    
    simuls=[];
    
end

nobj=numel(obj);

if nobj>1
    
    welf=cell(1,nobj);
    
    V=cell(1,nobj);
    
    d=cell(1,nobj);
    
    retcode=cell(1,nobj);
    
    for iobj=1:nobj
        
        [welf{iobj},retcode{iobj},V{iobj},d{iobj}]=loss(obj(iobj),simuls,varargin{:});
        
    end
    
    return
    
end

welf=[]; V=[]; d=[];

if isempty(obj)
    
    welf=cell(0,4);
    
    return
    
end

if ~(obj.is_optimal_policy_model||obj.is_optimal_simple_rule_model)
    
    error('model should have a loss function')
    
end

callers=dbstack;

callers={callers.name};

if any(strcmp(callers,'estimate'))
    
    retcode=0;
    
else
    
    [obj,retcode]=solve(obj,varargin{:});
    
end

if retcode
    
    return % error(decipher(retcode))
    
end

[Tx,Te]=set_solution_to_companion(obj);

h=numel(Tx);

W0=obj.solution.planner.weights;

Q=obj.solution.transition_matrices.Q;

beta=obj.solution.planner.discount;

[n,n]=size(Tx{1});

h=numel(Tx);

if numel(W0)==1 && h>1
    
    W0=W0(1,ones(1,h));
    
end

W=reshape(full(cell2mat(W0(:).')),[n,n,h]);

V0 = W; % V0=zeros(n,n,h);

if fslv
    
    options=struct('Display','iter',...
        'TolFun',obj.options.fix_point_TolFun,...
        'MaxIter',obj.options.fix_point_maxiter);
    
    [V,fval,exitflag]=fsolve(@fslvObjective,V0,options);
    
    exitflag=utils.optim.exitflag(exitflag,V,fval,obj.options.fix_point_TolFun);
    
    if exitflag~=1
        
        retcode=401;
        
    end
    
else
    
    B=W(:);
    
    [V,retcode]=transpose_free_quasi_minimum_residual(@objective,... % coefficient matrix
        B(:),... % right hand side
        V0(:),... % initial guess
        obj.options.fix_point_TolFun,... % tolerance level
        obj.options.fix_point_maxiter,... % maximum number of iterations
        obj.options.fix_point_verbose); % flag for printing progress or not
    
end

if retcode
    
    return % error(decipher(retcode))
    
end

V=reshape(full(V),[n,n,h]);

G=zeros(h,1);

for iii=1:h
    
    G(iii)=trace(V(:,:,iii)*(Te{iii}*Te{iii}.'));
    
end

b=cell2mat(beta);

bQ=bsxfun(@times,b(:),Q);

d=(eye(h)-bQ)\(bQ*G);

if isempty(simuls)
    
    [welf,retcode]=unconditional_welfare();
    
else
    
    simuls=pages2struct(simuls);
    
    [welf,retcode]=conditional_welfare();
    
end

    function [welf,retcode]=conditional_welfare()
        
        retcode=0;
        
        endo_names=obj.endogenous.name;
        
        span=get(simuls.(endo_names{1}),'NumberOfObservations');
        
        regimes=double(simuls.regime);
        
        x=nan(n,span);
        
        for iv=1:n
            
            x(iv,:)=double(simuls.(endo_names{iv}));
            
        end
        
        ss=obj.solution.ss;
        
        welf=0;
        
        for t=1:span
            
            rt=regimes(t);
            
            xt=x(:,t)-ss{rt};
            
            welf=welf+xt.'*V(:,:,rt)*xt+d(rt);
            
        end
        
        welf=welf/span;
        
    end

    function [welf,retcode]=unconditional_welfare()
        
        [ACV,retcode]=theoretical_autocovariances(obj,'autocov_ar',0,...
            'autocov_model_resolve',false);
        
        if retcode
            
            welf=[];
            
            return
            
        end
        
        % keeping only the first, i.e. the variance. But the calculation of the
        % covariance is not entirely satisfactory from a regime switching
        % perspective
        %----------------------------------------------------------------------
        ACV=ACV(:,:,ones(1,h));
        
        welf=nan(1,h);
        
        for ii=1:h
            
            welf(ii)=d(ii)+trace(V(:,:,ii)*ACV(:,:,ii));
            
        end
        
        %    EVV=sum(diag(W0{1}).*diag(ACV))/(1-beta{1});
        
        % EV == EVV?
    end

    function ax=fslvObjective(V)
        
        ax=objective(V);
        
        ax=reshape(ax,[n,n,h])-W;
        
    end

    function ax=objective(V)
        
        V=reshape(full(V),[n,n,h]);
        
        F=zeros(n,n,h);
        
        for ii=1:h
            
            for jj=1:h
                
                F(:,:,ii)=F(:,:,ii)+Q(ii,jj)*Tx{ii}.'*V(:,:,jj)*Tx{ii};
                
            end
            
            F(:,:,ii)=V(:,:,ii)-beta{ii}*F(:,:,ii);
            
        end
        
        ax=F(:);
        
    end

end
