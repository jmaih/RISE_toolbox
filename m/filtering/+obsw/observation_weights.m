function w=observation_weights(Z,T,H,Q,R,f,select,y)

if nargin<8
    
    y=[];
    
    if nargin<7
        
        select=[];
        
    end
    
end

debug=~isempty(y);

narginchk(6,8)

finalList={'a','att','v','r','alpha','epsilon','eta'};

if isempty(select)
    
    select=finalList;
    
else
    
    if ischar(select),select=cellstr(select); end
    
    bad=~ismember(select,finalList);
    
    if any(bad)
        
        disp(['Valid variables to decompose are :',cell2mat(strcat(finalList,'|'))])
        
        disp(select(bad))
        
        error('The items above are not valid elements to decompose')
        
    end
    
end

[~,p,n]=size(f.iF);

m=size(f.a,1);

r=size(R,2);

Im=eye(m);

[sp,Z,Z_is_selector]=obsw.time_series_length(Z,T,H,Q,R);

if Z_is_selector
    
    Zmat=Im(Z,:);
    
end

L=computeL();

w=struct();

do_a()

if any(strcmp(select,'v')), do_v(), end

if any(strcmp(select,'att')), do_att(), end

if any(strcmp(select,'r')), do_r(), end

if any(strcmp(select,'alpha')), do_alpha(), end

if any(strcmp(select,'epsilon')), do_epsilon(), end

if any(strcmp(select,'eta')), do_eta(), end

% discard unwanted
%------------------
fields=fieldnames(w);

bad=~ismember(fields,select);

badFields=fields(bad);

w=rmfield(w,badFields);

    function do_a()
        
        w.a=zeros((n+1)*m,n*p);
        
        for t=1:n+1
            
            rows=(t-1)*m+1:t*m;
            
            Ba=Im;
            
            for j=t-1:-1:1
                
                cols=(j-1)*p+1:j*p;
                
                w.a(rows,cols)=Ba*f.K(:,:,j);
                
                Ba=Ba*L(:,:,j);
                
            end
            
        end
        
        if debug
            
            disp('wa*y-a')
            
            disp(max(max(abs(reshape(w.a*y(:),m,[])-f.a))))
            
        end
    end

    function do_v()
        
        if ~isfield(w,'a')
            
            do_a()
            
        end
        
        w.v=zeros(n*p,n*p);
        
        Iy=speye(n*p);
        
        for t=1:n
            
            rows=(t-1)*p+1:t*p;
            
            at_stretch=(t-1)*m+1:t*m;
            
            if Z_is_selector
                
                Za=w.a(at_stretch(Z),:);
                
            else
                
                Zt=Z(:,:,sp.Z(t));
                
                Za=Zt*w.a(at_stretch,:);
                
            end
            
            w.v(rows,:)=Iy(rows,:)-Za;
            
        end
        
        if debug
            
            disp('wv*y-v')
            
            disp(max(max(abs(reshape(w.v*y(:),p,[])-f.v))))
            
        end
        
    end

    function do_att()
        
        if ~isfield(w,'v')
            
            do_v()
            
        end
        
        w.att=zeros(m*n,p*n);
        
        for t=1:n
            
            if Z_is_selector
                
                PZiF=f.P(:,Z,t)*f.iF(:,:,t);
                
            else
                
                PZiF=f.P(:,:,t)*Z(:,:,sp.Z(t)).'*f.iF(:,:,t);
                
            end
            
            rows=(t-1)*m+1:t*m;
            
            vrows=(t-1)*p+1:t*p;
            
            w.att(rows,:)=w.a(rows,:)+PZiF*w.v(vrows,:);
            
        end
        
        if debug
            
            disp('watt*y-att')
            
            disp(max(max(abs(reshape(w.att*y(:),m,[])-f.att))))
            
        end
        
    end

    function do_r()
        
        if ~isfield(w,'v')
            
            do_v()
            
        end
        
        w.r=zeros((n+1)*m,n*p);
        
        for t=n:-1:1
            
            rows=(t-1)*m+1:t*m;
            
            prev_rows=t*m+1:(t+1)*m;
            
            t_stretch=(t-1)*p+1:t*p;
            
            if Z_is_selector
                
                Zt=Zmat;
                
            else
                
                Zt=Z(:,:,sp.Z(t));
                
            end
            
            w.r(rows,:)=Zt.'*f.iF(:,:,t)*w.v(t_stretch,:)+...
                L(:,:,t).'*w.r(prev_rows,:);
            
        end
        
        if debug
            
            disp('wr*y-r')
            
            disp(max(max(abs(reshape(w.r*y(:),m,[])-f.r))))
            
        end
        
    end

    function do_alpha()
        
        if ~isfield(w,'r')
            
            do_r()
            
        end
        
        w.alpha=zeros(m*n,p*n);
        
        for t=1:n
            
            rows=(t-1)*m+1:t*m;
            
            w.alpha(rows,:)=w.a(rows,:)+f.P(:,:,t)*w.r(rows,:);
            
        end
        
        if debug
            
            disp('walpha*y-alpha')
            
            disp(max(max(abs(reshape(w.alpha*y(:),m,[])-f.alpha))))
            
        end
    end

    function do_epsilon()
        
        if ~isfield(w,'r')
            
            do_r()
            
        end
        
        w.epsilon=zeros(p*n,p*n);
        
        for t=1:n
            
            rows=(t-1)*p+1:t*p;
            
            mrows=t*m+1:(t+1)*m;
            
            w.epsilon(rows,:)=H(:,:,sp.H(t))*(f.iF(:,:,t)*w.v(rows,:)-f.K(:,:,t).'*w.r(mrows,:));
            
        end
        
        if debug
            
            disp('wepsilon*y-epsilon')
            
            disp(max(max(abs(reshape(w.epsilon*y(:),p,[])-f.epsilon))))
            
        end
        
    end

    function do_eta()
        
        if ~isfield(w,'r')
            
            do_r()
            
        end
        
        w.eta=zeros(r*n,p*n);
        
        for t=1:n
            
            rows=(t-1)*r+1:t*r;
            
            mrows=t*m+1:(t+1)*m;
            
            w.eta(rows,:)=Q(:,:,sp.Q(t))*R(:,:,sp.R(t)).'*w.r(mrows,:);
            
        end
        
        % check
        % reshape(w.eta*y(:),m,[])-f.eta
        
    end

    function L=computeL()
        % Ref. page 87
        L=zeros(m,m,n);
        
        for it=1:n
            
            if Z_is_selector
                
                Zt=Zmat;
                
            else
                
                Zt=Z(:,:,sp.Z(it));
                
            end
            
            L(:,:,it)=T(:,:,sp.T(it))-f.K(:,:,it)*Zt;
            
        end
        
    end

end