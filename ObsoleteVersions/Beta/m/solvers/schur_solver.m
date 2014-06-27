function [tt,rr,ss,retcode,genegvals]=schur_solver(Aminus,A0,Aplus,B,C,order)

if nargin<6
    order=1;
	if nargin<5
		C=[];
		if nargin<4
			B=[];
			if nargin<3
				Aminus=[];
				if nargin<2
					error([mfilename,':: number of arguments cannot be lower than 2'])
				end
			end
		end
	end
elseif nargin>6
	error([mfilename,':: number of arguments cannot exceed 6'])
end

tt=[];
rr=[];
ss=[];

n=size(Aplus,1);
if isempty(C)
	C=zeros(n,1);
end
if isempty(Aminus)
	Aminus=zeros(n);
end

F=[zeros(n),eye(n)
    -Aminus,-A0];
G=[eye(n),zeros(n)
    zeros(n),Aplus];

[T,S,Q,Z] = qz(F,G,'real');%complex

genegvals=nan(2*n,1);
ds=diag(S);
dt=diag(T);
genegvals(ds==0)=inf*sign(dt(ds==0));
genegvals(ds~=0)=dt(ds~=0)./ds(ds~=0);
[abs_eig,id]=sort(abs(genegvals));
genegvals=genegvals(id);
if sum(abs_eig<1)==n
	retcode=0;
%    disp([mfilename,':: unique stable solution'])
elseif sum(abs_eig<1)<n
	retcode=1;
%    disp([mfilename,':: no stable solution: solution is explosive'])
elseif sum(abs_eig<1)>n
	retcode=2;
%    disp([mfilename,':: multiple solutions: picking the most stable'])
end

if ~retcode
	[Tx,Sx,Qx,Zx] = ordqz(T,S,Q,Z,'udi'); %#ok<ASGLU>
	
	tt = Zx(n+1:end,1:n)/Zx(1:n,1:n);
	
	if ~isempty(B)
		ATA=-(Aplus*tt+A0)\eye(n);
		rr=ATA*B;
		for ii=2:order
			rr(:,:,ii)=ATA*Aplus*rr(:,:,ii-1);
		end
	end

	if any(C~=0)
		ss=-(Aplus+A0+Aminus)\C;
	else
		ss=C;
	end
end


