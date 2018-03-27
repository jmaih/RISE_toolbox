function linres=linear_restrictions(A,b)

linres=struct();

[na1,nparam]=size(A);

na2=nparam-na1;

linres.d=sparse(nparam,1);

linres.K=speye(nparam);

if na2==0
    
    linres.a2tilde_to_a=@(x,~)x;
    
    linres.a_to_a2tilde=@(x,~)x;
    
    return
    
end

do_checks()

[Q,R,E]=qr(A);

% [Q2,R2,evec]=qr(A,'vector'); ievec(evec)=1:nparam; % this does not seem
% to work well

iE=sparse(E\eye(nparam));

r=Q.'*b;

% partitioning
%-------------

R1=R(:,1:na1);

R2=R(:,na1+1:end);

if any(r)
    
    R1i_r_0=[R1\r;zeros(na2,1)];
    
else
    
    R1i_r_0=sparse(nparam,1);
    
end

R1i_R2_I2=[-R1\R2;eye(na2)];

linres.d=E*R1i_r_0;

linres.K=E*R1i_R2_I2;

linres.a2tilde_to_a=memoize_alpha(linres.d,linres.K);

linres.a_to_a2tilde=memoize_alpha2_tilde(iE,na1);

    function do_checks()
        
        if rank(full(A))~=na1
            
            error('Linear restriction matrix R (R*a=r) not of full rank. Probably some redundant restrictions')
        
        end
        % remove the rows with no restrictions
        bad_rows=~any(A,2);
        
        if any(bad_rows)
            
            disp(find(bad_rows))
            
            if any(b(bad_rows))
                
                error('no-restrictions rows with non-zero rhs')
                
            end
            
            warning('The no-restriction rows above are removed')
            
            A=A(~bad_rows,:);
            
            b=b(~bad_rows);
            
        end
        
    end

end

function out=memoize_alpha2_tilde(iE,na1)

out=@get_alpha2_tilde;

    function a2tilde=get_alpha2_tilde(a,covflag)
        
        if nargin<2,covflag=false; end
        
        % get atilde first then extract the relevant part
        %------------------------------------------------
        if covflag
            
            atilde=iE*a*iE.';
            
            a2tilde=atilde(na1+1:end,na1+1:end);
            
        else
            
            atilde=iE*a;
            
            a2tilde=atilde(na1+1:end);
            
        end
        
    end

end

function out=memoize_alpha(d,K)

out=@get_alpha;

    function a=get_alpha(a2tilde,covflag)
        
        if nargin<2
            
            covflag=false;
            
        end
        
        % get atilde first then re-order it
        %----------------------------------
        if covflag
            % the covariance may have several pages...
            npages=size(a2tilde,3);
            
            nn=size(K,1);
            
            a=nan(nn,nn,npages);
            
            for ii=1:npages
                
                a(:,:,ii)=K*a2tilde(:,:,ii)*K';
                
            end
            
        else
            
            a=d+K*a2tilde;
            
        end
        
    end

end
