function [cfkst]=condition_on_shocks(T,sstate,y0,states,state_cols,k,...
    nsteps,shocks,debug)
% cumul_matrices -- creates impact matrix for contemporaneous and future
% shocks
%
% Syntax
% -------
% ::
%
%   [M,ufkst]=cumul_matrices(T,sstate,y0,states,state_cols,k,nsteps,y_pos,e_pos)
%
% Inputs
% -------
%
% - **T** [{}]:
%
% - **sstate** [{}]:
%
% - **y0** [{}]:
%
% - **states** [{}]:
%
% - **state_cols** [{}]:
%
% - **k** [{}]:
%
% - **nsteps** [{}]:
%
% - **y_pos** [{}]:
%
% - **e_pos** [{}]:
%
% Outputs
% --------
%
% - **M** [{}]:
%
% - **ufkst** [{}]:
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

[ny,nz]=size(T{1});
h=size(T,2);
nx=numel(state_cols);
nshocks=(nz-(nx+1))/(k+1);
[nr,nc,np]=size(shocks);
if nr~=nshocks
    error('inconsistency of number of shocks')
end
[Tsig,Ty,Te]=re_separate_terms();
cfkst2=y0(:,ones(1,nsteps+1),ones(1,np));
cfkst=cfkst2;
for jstep=1:nsteps
    % pick the state
    %----------------
    st=states(jstep);
    
    for ip=1:np % loop through third dimension
        % first strategy : all in one go
        %----------------
        dev=cfkst(:,jstep,ip)-sstate{st};
        chop_max=min(jstep+k,nc);
        z=stateify(dev,shocks(:,jstep:chop_max,ip));
        rz=numel(z);
        cfkst(:,jstep+1,ip)=sstate{st}+T{st}(:,1:rz)*z;
        
        % second strategy: main part then shocks
        %---------------------------------------
        if debug
            cfkst2(:,jstep+1,ip)=sstate{st}+...
                Tsig{st}+...
                Ty{st}*(cfkst2(:,jstep,ip)-sstate{st});
            for icol=0:k
                t=jstep+icol;
                if t<=chop_max
                    which_shock=icol+1;
                    cfkst2(:,jstep+1,ip)=cfkst2(:,jstep+1,ip)+...
                        Te{st}{which_shock}*shocks(:,t,ip);
                end
            end
        end
    end
end

if debug
    disp('checking that both ways of computing forecast give the same answer')
    all(abs(cfkst2(:)-cfkst(:))<1e-9)
end

    function [Tsig,Ty,Te]=re_separate_terms()
        Tsig=cell(1,h);
        Ty=cell(1,h);
        Te=cell(1,h);
        dim1Dist=ny;
        dim2Dist=nshocks*ones(1,k+1);
        for ireg=1:h
            Ti=T{ireg};
            Tx=Ti(:,1:nx);
            Ty{ireg}=add_zeros(Tx);
            Tsig{ireg}=Ti(:,nx+1);
            Te{ireg}=mat2cell(Ti(:,nx+2:end),dim1Dist,dim2Dist);
        end
        function Ty=add_zeros(Tx)
            Ty=zeros(ny);
            Ty(:,state_cols)=Tx;
        end
    end
    function z=stateify(y0,shocks)
        x=y0(state_cols);
        % nullify sig and shocks when id>1
        %---------------------------------
        % form the state vector
        %-----------------------
        sig=1;
        z=[x
            sig
            shocks(:)];
    end
end