function [T,R,retcode,optim_opt,itercode,algo_name]=loose_commitment_solver(...
    Aminus,A0,Aplus,B,W,gam,betta,...
    solve_expect_order,reordering_index,T0,optim_opt)

% Reference: Debortoli, Maih, Nunes (2010): "Loose Commitment in Medium
% and Large-Scale Models"
% This version, June 13, 2011

optimization_options=struct('lc_reconvexify',false,...
    'lc_algorithm','short',...
    'lc_MaxIter',3000,...
    'lc_TolFun',sqrt(eps),...
    'lc_qz_criterium',1.000001,...
    'lc_verbose',false,...
    'lc_explosion_limit',1e+6);
if nargin==0
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    T=optimization_options;
    return
end
if nargin<11
    optim_opt=[];
    if nargin<10
        T0=[];
        if nargin<9
            reordering_index=[];
            if nargin<8
                solve_expect_order=1;
                if nargin <7
                    betta=1;
                    if nargin<6
                        gam=1;
                        if nargin<5
                            error([mfilename,':: at least arguments Aminus,A0,Aplus,B,C,Q should be passed'])
                        end
                    end
                end
            end
        end
    end
elseif nargin>11
    error([mfilename,':: number of arguments cannot exceed 10'])
end
algo_name=mfilename;
retcode=0;

if isempty(optim_opt)
    optim_opt=struct();
end

myfields=fieldnames(optimization_options);
for ii=1:numel(myfields)
    vi=myfields{ii};
    if ~isfield(optim_opt,vi)
        optim_opt.(vi)=optimization_options.(vi);
    end
    tmp=optim_opt.(vi); %#ok<NASGU>
    eval([vi,'=tmp;'])
end

if gam==1
    lc_algorithm='long';
end

[mult_nbr,endo_nbr]=size(A0);
m=endo_nbr+mult_nbr;

if isempty(T0)||~isequal(size(T0),[m,m])
    T0=zeros(m);
end
% solving for T
iter=0;
conv_T=0.1*lc_explosion_limit;
while conv_T>lc_TolFun && iter<lc_MaxIter && conv_T<lc_explosion_limit && ~any(any(isnan(T0)))
    iter=iter+1;
    %===== Iterate =====%
    switch lc_algorithm
        case 'short'
            if iter==1
                [GAM0,GAMlag,GAMv]=loose_commitment_matrices(T0,Aminus,A0,Aplus,B,W,betta,gam,lc_reconvexify);
                if det(GAM0)==0
                    [GAM0,GAMlag,GAMv]=loose_commitment_matrices(rand(size(T0)),Aminus,A0,Aplus,B,W,betta,gam,lc_reconvexify);
                end
            else
                GAM0=loose_commitment_matrices(T0,Aminus,A0,Aplus,B,W,betta,gam,lc_reconvexify);
            end
            T=-GAM0\GAMlag;
        case 'long'
            if iter==1
                [OMG0,OMGlag,OMGlead,GAMv]=loose_commitment_matrices_long(T0,Aminus,A0,Aplus,B,W,betta,gam,lc_reconvexify);
            else
                OMG0=loose_commitment_matrices_long(T0,Aminus,A0,Aplus,B,W,betta,gam,lc_reconvexify);
            end
            % 		[T,control]=msre_aim(0,OMGlead,OMG0,OMGlag,0,m,1);
%             [T,control]=msre_gensys(0,{OMGlead},{OMG0},{OMGlag},0,m,1);
            [T,control]=msre_klein(0,{OMGlead},{OMG0},{OMGlag},0,m,1);

            if isnan(control)
                T=nan;
            else
                T=T{1};
            end
            if gam==1
                T0=T;
            end
    end
    %====================%
    conv_T=max(max(max(abs(T-T0))));
    if lc_verbose
        fprintf(1,'%8.0f %12.4f\n',iter,conv_T);
    end
    T0=T;
end
itercode=iter;
if iter>=lc_MaxIter
    retcode=21;
    if lc_verbose
        disp([mfilename,':: Maximum Number of Iterations reached'])
    end
elseif any(any(isnan(T0)))
    retcode=22;
    if lc_verbose
        disp([mfilename,':: NAN elements in the solution'])
    end
elseif conv_T>=lc_explosion_limit
    retcode=23;
    if lc_verbose
        disp([mfilename,':: explosive solution'])
    end
elseif max(abs(eig(T)))>lc_qz_criterium
    retcode=24;
    if lc_verbose
        disp([mfilename,':: Some eigenvalues greater than lc_qz_criterium, Model potentially unstable'])
    end
end


if retcode
    T=[];
    R=[];
else
    if strcmp(lc_algorithm,'long')
        clear OMGlead OMG0 OMGlag
        GAM0=loose_commitment_matrices(T0,Aminus,A0,Aplus,B,W,betta,gam,lc_reconvexify);
    end
    
    iGAM0=GAM0\eye(m);
    [~,exo_nbr,nregs]=size(GAMv);
    R=nan(m,exo_nbr,solve_expect_order,nregs);
    for ireg=1:nregs
        R(:,:,1,ireg)=-iGAM0*GAMv(:,:,ireg);
        if solve_expect_order>1
            AGAiA=iGAM0*[zeros(endo_nbr),betta*gam*transpose(Aminus)
                Aplus,zeros(mult_nbr)];
            for h=2:solve_expect_order
                R(:,:,h,ireg)=-AGAiA*R(:,:,h-1,ireg);
            end
        end
    end
    
    if ~isempty(reordering_index)
        T=T(reordering_index,reordering_index);
        R=R(reordering_index,:,:,:);
    end
    
end

end


function [GAM0,GAMlag,GAMv]=loose_commitment_matrices(T,Aminus,A0,Aplus,B,W,betta,gam,lc_reconvexify)

[mult_nbr,endo_nbr]=size(A0);
m=mult_nbr+endo_nbr;
y=1:endo_nbr;
l=endo_nbr+1:m;

GAM0=zeros(m);

Tyy=T(y,y);
Tyl=T(y,l);
Tly=T(l,y);
Tll=T(l,l);

AplusT=Aplus*Tyy;

GAM0(y,y)=2*W+betta*Aminus'*Tly;
GAM0(y,l)=transpose(gam*betta*Tll'*Aminus+A0+(1-gam)*AplusT);

GAM0(l,y)=A0+AplusT;
GAM0(l,l)=gam*Aplus*Tyl;

if nargout>1
    GAMlag=zeros(m);
    if lc_reconvexify
        GAMlag(y,l)=gam/betta*Aplus';
    else
        GAMlag(y,l)=(gam>0)/betta*Aplus';
    end
    GAMlag(l,y)=Aminus;
    if nargout==3
        GAMv=zeros(m,size(B,2),size(B,3));
        GAMv(l,:,:)=B;
    end
end
end

function [OMG0,OMGlag,OMGlead,GAMv]=loose_commitment_matrices_long(T,Aminus,A0,Aplus,B,W,betta,gam,lc_reconvexify)
[mult_nbr,endo_nbr]=size(A0);
m=endo_nbr+mult_nbr;
y=1:endo_nbr;
l=endo_nbr+1:m;
Tyy=T(y,y);
Tly=T(l,y);

OMG0=zeros(m);
OMG0(y,y)=2*W+betta*(1-gam)*transpose(Aminus)*Tly;
OMG0(y,l)=transpose(A0+(1-gam)*Aplus*Tyy);
OMG0(l,y)=A0+(1-gam)*Aplus*Tyy;

if nargout>1
    OMGlead=zeros(m);
    OMGlead(l,y)=gam*Aplus;
    OMGlead(y,l)=betta*gam*transpose(Aminus);
    
    OMGlag=zeros(m);
    OMGlag(l,y)=Aminus;
    if lc_reconvexify
        OMGlag(y,l)=gam/betta*Aplus';
    else
        OMGlag(y,l)=(gam>0)/betta*Aplus';
    end
    if nargout>3
        GAMv=zeros(m,size(B,2),size(B,3));
        GAMv(l,:,:)=B;
    end
end
end
