function [yf,info]=conditional_forecast(y0,Xh,params,Rfunc,yh,shock_uncertainty)

Aj=[];

Rj=[];

info={'variables','horizon','replications'};

np=numel(params);

for jj=1:np
    
    load_params();
    
    if jj==1
        
        yf=do_one(Aj,Rj,y0,Xh,yh,shock_uncertainty);
        
        yf=yf(:,:,ones(1,np));
        
    else
        
        yf(:,:,jj)=do_one(Aj,Rj,y0,Xh,yh,shock_uncertainty);
        
    end

end

    function load_params()
        % identification can be expensive so save some steps if possible,
        % for instance in the case
        Aj=params(jj).B;
        
        Rj=Rfunc(params(jj));
        
    end

end

function yf=do_one(A,B,y0,Xh,yh,shock_uncertainty)

nx=size(Xh,1);

if nx && any(isnan(Xh(:)))
    
    error('deterministic variables cannot be nan')
    
end

n=size(A,1);

ns=size(B,2);

[~,h]=size(yh);

nnul=[];

[M1,M2,RM2i,mu]=null_colum_spaces(A,B,y0);

g2=RM2i*mu;

if shock_uncertainty
    
    g1=randn(nnul,1);
    
else
    
    g1=zeros(nnul,1);
    
end

e=reshape(M1*g1+M2*g2,ns,[]);

yf=vartools.simulate(y0,Xh,A,B,e);

    function [M1,M2,RM2i,mu]=null_colum_spaces(A,B,y0)
        % separate deterministic terms from the rest
        C=A(:,1:nx);
        A=A(:,nx+1:end);
        
        % Form the companion, explicitly excluding deterministic terms
        [A,B]=vartools.companion(A,B,0);
        
        na=size(A,1);
        
        % add zeros under C to account for companionship        
        C=[C;zeros(na-n,nx)];
        
        AB0=cell(1,h); AC0=cell(1,h);
        
        y0=fliplr(y0);
                
        Ay=cell(h,1); ay=A*y0(:);
        
        ab=B; ac=C;
        
        for t=h:-1:1
            %
            AB0{t}=ab(1:n,:);
            ab=A*ab;
            %
            AC0{t}=ac(1:n,:);
            ac=A*ac;
            % Note that the loop is starting from the end and working
            % backwards!!!
            Ay{h-t+1}=ay(1:n);
            ay=A*ay;
        end
        
        AB0=cell2mat(AB0); AC0=cell2mat(AC0); Ay=cell2mat(Ay);
        
        R=zeros(h*n,h*ns); AC=zeros(h*n,h*nx);
        
        for t=1:h
            rr=(t-1)*n+1:t*n;
            %
            R(rr,1:t*ns)=AB0(:,(h-t)*ns+1:end);
            %
            AC(rr,1:t*nx)=AC0(:,(h-t)*nx+1:end);
        end
        
        vyh=yh(:);
        
        good=~isnan(vyh);
        
        vyh=vyh(good);
        
        R=R(good,:);
        
        AC=AC(good,:);
        
        Ay=Ay(good);
        
        mu=vyh-Ay-AC*Xh(:);
        
        nconstr=size(R,1);
        
        nnul=size(R,2)-nconstr;
        
        [M1,M2]=null_col(R);
        
        RM2i=(R*M2)\eye(nconstr);
        
    end

end

function [M1,M2]=null_col(R)
        
[~,S,V]=svd(R);

[q,cols]=size(S);

% another basis which is different but gives the same end results is
% M22=null(M1')
M2=V(:,1:q); 

M1=V(:,q+1:cols); % == null(R)

%         M1=null(R); M2=null(M1.');
end