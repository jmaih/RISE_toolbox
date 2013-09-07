function [crit,FT,FR]=markov_switching_dsge_objective(Q,Aplus,A0,Aminus,B,T,R)
if nargin<7
    R=[];
end
if isempty(R)
    if ~isempty(B)
        error([mfilename,':: B cannot be given without R'])
    end
end
if ~isempty(B)
    [n,x,~,h]=size(B);
else
    x=0;
    [~,n,h]=size(T);
end
FT=zeros(n*h,n);
FR=zeros(n*h,x);
Penalty=0;
for st=1:h
    QT=0;
    for stlead=1:h
        if st~=stlead
            QT=QT+Q(st,stlead)*T(:,:,stlead);
        end
    end
	App=Q(st,st)*Aplus(:,:,st);
	A00=Aplus(:,:,st)*QT+A0(:,:,st);
    Amm=Aminus(:,:,st);
    AXB=App*T(:,:,st)+A00;
    FT((st-1)*n+1:st*n,:)=AXB*T(:,:,st)+Amm;
    if ~isempty(R)
        FR((st-1)*n+1:st*n,:)=AXB*R(:,:,1,st)+B(:,:,st);
    end
    Penalty=Penalty+max(max(abs(T(:,:,st)+AXB\Amm)));
end
crit=max(max(abs(FT)))+Penalty;
if ~isempty(R)
    crit=max(crit,max(max(abs(FR))));
end

