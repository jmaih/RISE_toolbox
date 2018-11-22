function [myirfs,info]=irf(self,irf_periods,params,identification,girf_setup)
% INTERNAL FUNCTION
%

if isempty(identification)
    
    identification=vartools.choleski(self.nvars);
    
end

is_girf=~isempty(girf_setup);

if is_girf
    
    d=girf_default_setup(girf_setup);
    
    nsims=d.nsims;
    
end

kreps=numel(params);

nshocks=self.nvars*self.ng;

nvars=self.nvars*self.ng;

nx=self.nx*self.ng;

h=size(params(1).B,3);

fifth_dimension=h*(~is_girf)+is_girf;

myirfs=zeros(nvars,irf_periods,nshocks,kreps,fifth_dimension);

y0=zeros(nvars,self.nlags);

failed=false(1,kreps);

% feed in zeros for deterministic variables
%-------------------------------------------
xdet=zeros(nx,irf_periods);

retcode=zeros(1,kreps);

for jj=1:kreps
    
    Aj=params(jj).B;
    
    [Rj,retcode(jj)]=identification(params(jj));
    
    if retcode(jj)
        
        failed(jj)=true;
        
        continue
        
    end
    
    % Only regime-specific IRFs are computed
    %---------------------------------------
    Qfunc=transition_function(params(jj).Q.Q);
    
    vartools.check_factorization(Rj,params(jj).S)
    
    for ishock=1:nshocks
        
        if is_girf
            
            myirfs(:,:,ishock,jj)=girfs();
            
        else
            
            myirfs(:,:,ishock,jj,:)=simple_irfs();
            
        end
        
    end
    
end

rcode=unique(retcode);

if all(rcode)
    
    for ii=1:numel(rcode)
        
        decipher(rcode(ii))
        
    end
    
    error('Irfs could not be computed due to the errors above')
    
end

myirfs=myirfs(:,:,:,~failed,:);

info={'nvars','length','nshocks','nrepetitions','nregimes'};

info=info(1:4+~is_girf);


    function rf=girfs()
        
        rf=zeros(nvars,irf_periods);
        
        c=1/nsims;
        
        for isim=1:nsims
            
            shocks=randn(nshocks,irf_periods);
            
            shocks(ishock,1)=1;
            
            [rf1,regs1]=vartools.simulate(y0,xdet,Aj,Rj,shocks,Qfunc);
            
            shocks(ishock,1)=0;
            
            rf2=vartools.simulate(y0,xdet,Aj,Rj,shocks,Qfunc,regs1);
            
            rf12=rf1-rf2;
            
            rf=rf+c*rf12;
            
        end
        
    end


    function rf=simple_irfs()
        
        rf=zeros(nvars,irf_periods,1,1,h);
        
        for ireg=1:h
            
            shocks=zeros(nshocks,irf_periods);
            
            shocks(ishock,1)=1;
            
            rf(:,:,1,1,ireg)=vartools.simulate(y0,xdet,...
                Aj(:,:,ireg),Rj(:,:,ireg),shocks,Qfunc);
            
        end
        
    end


    function out=transition_function(Q0)
        
        out=@engine;
        
        function [Q,retcode]=engine(~)
            
            Q=1;
            
            retcode=0;
            
            if is_girf
                
                Q=Q0;
                
            end
            
        end
        
    end

end


function d=girf_default_setup(girf_setup)

if ~isstruct(girf_setup)
    
    error('girf_setup must be a structure')
    
end

numfint=@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&floor(x)==ceil(x);

d={'nsims',300,@(x)numfint(x),'nsims must be a positive and finite integer'};

d=cell2struct(d,{'name','default','check','errmsg'},2);

f={d.name};

for ii=1:numel(f)
    
    fi=f{ii};
    
    if isfield(girf_setup,fi)
        
        gfi=girf_setup.(fi);
        
        assert(d(ii).check(gfi),d(ii).errmsg)
        
        d.(ii).default=gfi;
        
        girf_setup=rmfield(girf_setup,fi);
        
    end
    
end

f=fieldnames(girf_setup);

if ~isempty(f)
    
    disp(f)
    
    error('the fields above are not valid fields for girf_setup')
    
end

d=rmfield(d,{'check','errmsg'});

d=struct2cell(d);

d=cell2struct(d(2,:),d(1,:),2);

end
