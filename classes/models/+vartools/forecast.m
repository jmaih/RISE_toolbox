function [fkst,info]=forecast(y0,xdet,params,shocks,identification)

if nargin < 5
    
    identification=[];
    
    if nargin < 4
        
        shocks=[];
        
    end
    
end

uncertainty=false;

nperiods=nan;

nvars=size(params(1).S,1);

ny=size(params(1).B,1);

% [nyy,nlags]=size(y0);
        
% nlags=floor(kpx/ny);

check_inputs_for_error();

par_replic=numel(params);

shocks_to_shocks()

if isempty(identification)
    
    identification=vartools.choleski(nvars);
    
end

shocks_replic=size(shocks,3);

kreps=max(par_replic,shocks_replic);

fkst=zeros(nvars,nperiods,kreps);

Aj=[];

Rj=[];

for jj=1:kreps
    
    jpar=min(jj,par_replic);
    
    load_params();
    
    vartools.check_factorization(Rj,params(jpar).S)
    
    jshk=min(jj,shocks_replic);
    
    fkst(:,:,jj)=vartools.simulate(y0,xdet,Aj,Rj,shocks(:,:,jshk));
    
end

info={'nvars','length','nrepetitions'};

    function load_params()
        % identification can be expensive so save some steps if possible,
        % for instance in the case
        if jpar==jj
            
            Aj=params(jpar).B;
            
            Rj=identification(params(jpar));
            
        end
    
    end

    function check_inputs_for_error()
                
        if nvars~=ny
            
            error('number of variables does not match in y0 and in S')
            
        end
        
    end

    function shocks_to_shocks()
        
        if isempty(shocks)
            
            shocks=12;
            
        end
        
        [nys,nper,sreplic]=size(shocks);
        
        if nper==1 && nys==1 && sreplic==1
            % length of forecasts only
            uncertainty=true;
            
            nperiods=shocks;
            
        elseif nys==nvars && nper>1
            % shocks given
            nperiods=nper;
            
        end
        
        if uncertainty
            
            shocks=randn(nvars,nperiods,par_replic);
            
        end
        
    end

end
