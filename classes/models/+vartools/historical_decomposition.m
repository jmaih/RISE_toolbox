function [Histdec,info,contrib_positions]=historical_decomposition(AA,BB,nz,X,Resids)

% collect the state matrices
n=size(AA,1);
    
% deterministic terms
xx_exogenous=X(1:nz,:);

% Initial conditions
y0=X(nz+1:end,1);

% separate deterministic terms from the rest
C=AA(:,1:nz);
AA=AA(:,nz+1:end);

% get the residuals before companioning BB
smoothed_shocks=BB\Resids;

% Form the companion, explicitly excluding deterministic terms
[AA,BB]=vartools.companion(AA,BB,0);

% add zeros under C to account for companionship
vvv=size(AA,1)-n;

C=[C;zeros(vvv,nz)];

[D,y0_id,z_id,shocks_id]=decomposition_engine(y0,xx_exogenous,smoothed_shocks);

% t0 and retain only the relevant rows
%--------------------------------------
D=D(1:n,2:end,:);

Histdec=permute(D,[2,3,1]);

info={'horizon','contributors','variables'};

contrib_positions=struct('y0',y0_id,'det_vars',z_id,'shocks',shocks_id);

    function [D,y0_id,z_id,shocks_id]=decomposition_engine(y0,xx,shocks)
        
        ny=size(AA,1);
        
        [nx,nt]=size(shocks);
        
        nz=size(C,2);
        
        y0_id=1;
        
        z_id=y0_id+(1:nz);
        
        shocks_id=y0_id+nz+(1:nx);
        
        D=zeros(ny,nt+1,1+nz+nx);
        
        D(:,1,y0_id)=y0;
        
        for t=1:nt
            
            % initial conditions
            %--------------------            
            D(:,t+1,y0_id)=AA*D(:,t,y0_id);
            
            % deterministic terms: There could be many deterministic
            % variables and so, move the second dimension to the end
            %--------------------------------------------------------
            Z0=permute(D(:,t,z_id),[1,3,2]);
            
            D(:,t+1,z_id)=AA*Z0+bsxfun(@times,C,xx(:,t).');
            
            % shocks
            %--------
            Sh0=permute(D(:,t,shocks_id),[1,3,2]);
            
            D(:,t+1,shocks_id)=AA*Sh0+bsxfun(@times,BB,shocks(:,t).');
            
        end
        
    end

end
