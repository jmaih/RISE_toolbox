function AllSols=grobner_solve(Aplus,A0,Aminus,aplus,aminus)
% solves for all possible solutions of a polynomial
% Aplus,A0,Aminus are cell arrays Aplus is h x h while the others are 1 x h
% inside each cell, A0 are all square
% aplus,aminus are logicals indicating which columns of the square matrix
% are filled.

if nargin<5
    
    aminus=[];
    
    if nargin<4
        
        aplus=[];
        
    end
    
end

h=size(Aplus,1);

n0=size(A0{1},1);

if isempty(aplus),aplus=true(1,n0); end

if isempty(aminus),aminus=true(1,n0); end

% np=sum(aplus);

n_=sum(aminus);

if n_==0
    
    AllSols={zeros(1,n_)};
    
    return
    
end

[T,Ti]=create_symbols(n0,n_,h);

if license('checkout','symbolic_toolbox')&& ~verLessThan('matlab','9.4.0')
    
    AllSols=use_symbolic_toolbox(T,Ti,A0,Aminus,Aplus,aplus,aminus);
    
else
    
    AllSols=use_petschel(T,Ti,A0,Aminus,Aplus,aplus,aminus);
    
end

end

function AllSols=use_petschel(T,Ti,A0,Aminus,Aplus,aplus,aminus)

cellmult=@utils.miscellaneous.cellmult;

T=utils.miscellaneous.cell_str2cell(T);

h=numel(A0);

vnames=Ti(:).';

n=numel(vnames);

proto=zeros(1,n+1);

pols=cell(1,h);

% polynomials
%------------
for s0=1:h
    
    pols{s0}=set_polynomial(A0{s0},T{s0});
    
    pols{s0}=set_constants(pols{s0},Aminus{s0});
    
    for s1=1:h
        
        ATTs10=cellmult(T{s1}(aplus,:),T{s0}(aminus,:));
        
        pols{s0}=set_polynomial(Aplus{s0,s1},ATTs10,pols{s0});
        
    end
    
end

pols=[pols{:}];

sol=petschel.polynsolve(pols(:).','lex',vnames,sqrt(eps));

nsol=size(sol,1);

AllSols=repmat({zeros(size(Ti))},1,nsol);

for ii=1:nsol
    
    soli=sol(ii,:);
    
    for jj=1:numel(Ti)
        
        AllSols{ii}(jj)=soli(jj);
        
    end
    
end

    function out=set_constants(out,C)
        
        [rc,cc]=size(C);
        
        tmp=proto;
        
        for iii=1:rc
            
            for jjj=1:cc
                
                cij=C(iii,jjj);
                
                if cij==0
                    
                    continue
                    
                end
                
                tmp(1)=cij;
                
                out{iii,jjj}=[out{iii,jjj};tmp];
                
            end
            
        end
        
    end

    function out=set_polynomial(A,T0,out)
        
        if nargin<3
            
            out=[];
            
        end
        
        T1=rip_apart(T0);
        
        [ra,ca]=size(A);
        
        [~,ct]=size(T1);
        
        if isempty(out)
            
            out=cell(ra,ct);
            
        end
        
        for icolt=1:ct
            
            for irowa=1:ra
                
                for icola=1:ca
                    
                    tmp=proto;
                    
                    tmp(1)=A(irowa,icola);
                    
                    if tmp(1)~=0
                        
                        V=T1{icola,icolt};
                        
                        for kk=1:size(V,1)
                            
                            tmp(V(kk,1)+1)=V(kk,2);
                            
                        end
                        
                        out{irowa,icolt}=[out{irowa,icolt};tmp];
                        
                    end
                    
                end
                
                %%%% out{irowa,icolt}=sparse(out{irowa,icolt});
                
            end
            
        end
        
    end

    function T=rip_apart(T)
        
        for iii=1:numel(T)
            
            T{iii}=regexp(T{iii},...
                '(?<var>\<[a-zA-Z]+\w*\>)(\^)?(?<power>\d+)?','names');
            
            nv=numel(T{iii});
            
            X=nan(nv,2);
            
            for jjj=1:nv
                
                st=T{iii}(jjj);
                
                if isempty(st.power)
                    
                    st.power='1';
                    
                end
                
                pos=find(strcmp(st.var,vnames));
                
                X(jjj,:)=[pos,str2double(st.power)];
                
            end
            
            T{iii}=X;
            
        end
        
    end

end


function AllSols=use_symbolic_toolbox(T,Ti,A0,Aminus,Aplus,aplus,aminus)

syms(Ti{:});

vars=cell2mat(strcat(Ti(:).',','));

vars=eval(['[',vars(1:end-1),']']);

T=cellfun(@eval,T,'uniformOutput',false);

h=numel(A0);

% polynomials
%------------
for s0=1:h
    
    f0=A0{s0}*T{s0}+Aminus{s0};
    
    for s1=1:h
        
        f0=f0+Aplus{s0,s1}*T{s1}(aplus,:)*T{s0}(aminus,:);
        
    end
    
    if s0==1
        
        pols=f0(:).';
        
    else
        
        pols=[pols,f0(:).']; %#ok<AGROW>
        
    end
    
end

% Groebner basis: pols must be a row vector
%------------------------------------------
B=gbasis(pols,vars);

Sol0=solve(B);

Sol = structfun(@double, Sol0,'UniformOutput',false);

nsol=numel(Sol.T111);

AllSols=repmat({zeros(size(Ti))},1,nsol);

for ii=1:nsol
    
    for jj=1:numel(Ti)
        
        AllSols{ii}(jj)=Sol.(Ti{jj})(ii);
        
    end
    
end

end


function [T,Ti]=create_symbols(n0,n_,h)

T=cell(1,h);

Ti=cell(n0,n_,h);

for ireg=1:h
    
    T{ireg}='[';
    
    for jrow=1:n0
        
        for kcol=1:n_
            
            xx=sprintf('T%0.0f%0.0f%0.0f',ireg,jrow,kcol);
            
            T{ireg}=[T{ireg},xx];
            
            if kcol==n_
                
                if jrow<n0
                    
                    T{ireg}=[T{ireg},';'];
                    
                end
                
            else
                
                T{ireg}=[T{ireg},','];
                
            end
            
            Ti{jrow,kcol,ireg}=xx;
            
        end
        
    end
    
    T{ireg}=[T{ireg},']'];
    
end

end