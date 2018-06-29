function [myirfs,info]=irf(kdata,irf_periods,params,identification)

% if nargin < 4
% 
%     identification=[];
% 
% end

if isempty(identification)
    
    identification=vartools.choleski(kdata.nvars);
    
end

kreps=numel(params);

nshocks=kdata.nvars*kdata.ng;

nvars=kdata.nvars*kdata.ng;

nx=kdata.nx*kdata.ng;

h=size(params(1).B,3);

myirfs=zeros(nvars,irf_periods,nshocks,kreps,h);

shocks0=zeros(nshocks,irf_periods);

y0=zeros(nvars,kdata.nlags);

failed=false(1,kreps);

% feed in zeros for deterministic variables
%-------------------------------------------
xdet=zeros(nx,irf_periods);

% Only regime-specific IRFs are computed
%---------------------------------------
Qfunc=transition_function();

for jj=1:kreps
    
    Aj=params(jj).B;
        
    [Rj,retcode]=identification(params(jj));
    
    if retcode
        
        failed(jj)=true;
        
        continue
        
    end
    
    vartools.check_factorization(Rj,params(jj).S)
        
    for ishock=1:nshocks
        
        for ireg=1:h
            
            shocks=shocks0;
            
            shocks(ishock,1)=1;
            
            myirfs(:,:,ishock,jj,ireg)=vartools.simulate(y0,xdet,...
                Aj(:,:,ireg),Rj(:,:,ireg),shocks,Qfunc);
            
        end
        
    end
    
end

myirfs=myirfs(:,:,:,~failed,:);

info={'nvars','length','nshocks','nrepetitions','nregimes'};

    function out=transition_function()
        
        out=@engine;
        
        function [Q,retcode]=engine(~)
            
            Q=1;
            
            retcode=0;
            
        end
        
    end

end
