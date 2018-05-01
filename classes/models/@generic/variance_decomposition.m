function [Vardec,obj]=variance_decomposition(obj,varargin)
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


% PURPOSE: Computes variance decompositions of a MSRE model
%-------------------------------------------------------------
% USAGE:
% where:
%
%
%
%
%-------------------------------------------------------------
% EXAMPLE:
%-------------------------------------------------------------
% LOG:
% 1.
%
%
%-------------------------------------------------------------
% SEE ALSO:
%-------------------------------------------------------------

% written [December 28, 2010] by:
% Junior Maih, Dept of Economics
% The Central Bank of Norway
% junior.maih@norges-bank.no
% this update [July 18, 2014]

if isempty(obj)   
    
    mydefaults=the_defaults();
    
    if nargout
        
        Vardec=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

nobj=numel(obj);

if nobj>1
    
    Vardec=cell(1,nobj);
    
    for iobj=1:nobj
        
        [Vardec{iobj},obj(iobj)]=variance_decomposition(obj(iobj),varargin{:});
        
    end
    
    return
    
end

obj=set(obj,varargin{:});

[obj,retcode]=solve(obj);	

% extract the state and transition matrices
if retcode
    
    error([mfilename,':: model could not be solved'])
    
end

vardec_ergodic=obj.options.vardec_ergodic;

vardec_theoretical  =obj.options.vardec_theoretical  ;

vardec_periods  =obj.options.vardec_periods  ;

k=max(100,max(vardec_periods));

% collect the state matrices
Q=obj.solution.transition_matrices.Q;

endo_nbr=obj.endogenous.number;

exo_nbr=sum(obj.exogenous.number);

horizon=1;

if isa(obj,'dsge')
    
    horizon=max(obj.exogenous.shock_horizon(:))+1;
    
end

nregs=obj.markov_chains.regimes_number;

[T,R]=set_solution_to_companion(obj);

if vardec_ergodic
    
    [T,R]=expected_state_matrices(T,R);
    
end

endo_names=obj.endogenous.name;

exo_names=obj.exogenous.name;

if isa(obj,'rfvar') && ~isempty(obj.identification)
    % the names of the exogenous has to change, especially if the model has
    % been identified
    %----------------------------------------------------------------------
    exo_names=obj.structural_shocks.name;
    
end

is_observed=obj.exogenous.is_observed;

for ishock=1:numel(exo_names)
    
    if is_observed(ishock)
        
        for istate=1:nregs
            
            R{istate}(:,ishock,:)=0;
            
        end
        
    end
    
end

Vardec=struct();

if vardec_theoretical
    
    if vardec_ergodic || nregs==1
        
        [Vinfi,Vi]=theoretical_vardec_engine(T{1},R{1});
        
        Vardec.conditional=vardec2rise_time_series(Vi);
        
        Vardec.infinity=vardec2rise_time_series(Vinfi);
        
    else
        
        for istate=1:nregs
            
            [Vinfi,Vi]=theoretical_vardec_engine(T{istate},R{istate});
            
            Vardec(istate).conditional=vardec2rise_time_series(Vi);
            
            Vardec(istate).infinity=vardec2rise_time_series(Vinfi);
            
        end
        
    end
    
else
    % compute from simulation
    error('Simulated variance decomposition not yet implemented')
    
end

    function db=vardec2rise_time_series(V)
        
        db=struct();
        
        this_k=size(V,3);
        
        for iendo=1:endo_nbr
            
            thisdec=transpose(reshape(V(iendo,:,:),exo_nbr,this_k));
            
            db.(endo_names{iendo})=ts(1,thisdec,exo_names);
            
        end
        
    end

    function [Vinfi,Vi]=theoretical_vardec_engine(T,R)
        
        grand_endo_nbr=size(R,1);
        
        R=reshape(R,grand_endo_nbr,exo_nbr*horizon);
        
        V=zeros(grand_endo_nbr);
        
        Vi=zeros(grand_endo_nbr,exo_nbr,k);
        
        Vii=zeros(grand_endo_nbr,grand_endo_nbr,exo_nbr);
        
        RR=R*R';
        
        for i=1:k
            
            V=T*V*T'+RR;
            
            [Vi(:,:,i),Vii]=decompose_variance(Vi(:,:,i),V,Vii);
            
        end
        
        Vinfi=zeros(grand_endo_nbr,exo_nbr);
        
        [Vinf,retcode0]=lyapunov_equation(T,RR,obj.options);
        
        if retcode0 || any(~isfinite(Vinf(:)))
            
            error('Variance could not be solved')
            
        end
        
        % deal with zero variances
        Vinfi=decompose_variance(Vinfi,Vinf);
        
        function [V,Vkk]=decompose_variance(V,total_variance,Vkk)
            
            total_variance=diag(total_variance);
            
            total_variance(total_variance<1e-12)=1;
            
            for iexo=1:exo_nbr
                
                Ri=zeros(size(R));
                
                locs=iexo:exo_nbr:exo_nbr*horizon;
                
                Ri(:,locs)=R(:,locs);
                
                RRi=Ri*Ri';
                
                if nargin<3
                    
                    [V00,retcode1]=lyapunov_equation(T,RRi,obj.options);
                    
                    if retcode1 && any(~isfinite(V00(:)))
                        
                        error('Variance could not be solved')
                        
                    end
                    
                    Vkk=V00;
                    
                else
                    
                    Vkk(:,:,iexo)=T*Vkk(:,:,iexo)*T'+RRi;
                    
                    V00=Vkk(:,:,iexo);
                    
                end
                
                V(:,iexo)=diag(V00)./total_variance;
                
            end
            
        end
        
    end

    function [A,B]=expected_state_matrices(T,R)
        
        probs_t=[eye(nregs)-Q';ones(1,nregs)]\[zeros(nregs,1);1];
        
        A=0;
        
        B=0;
        
        for st=1:nregs
            
            A=A+probs_t(st)*T{st};
            
            B=B+probs_t(st)*R{st};
            
        end
        
        A={A};B={B};
        
    end

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

d={
    'vardec_shocks','',@(x)ischar(x)||iscellstr(x),'vardec_shocks must be char or cellstr'

    'vardec_periods',[1 4 8 16 40 100 200 400],...
    @(x)isvector(x) && all(num_fin_int(x)),...
    'all vardec_periods must be a finite and positive integer'

    'vardec_theoretical',true,@(x)islogical(x),'vardec_theoretical must be a logical'
    
    'vardec_ergodic',false,@(x)islogical(x),'vardec_ergodic must be a logical'
    };

end
