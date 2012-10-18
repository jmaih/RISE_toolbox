function flag=is_eigenvalue_solver_candidate(Aplus,A0,Aminus,Q)
flag=true;
[junk,junk,h]=size(A0);
if ~(h==1 || all(diag(Q)==1))
	A0test=A0(:,:,1);
	Aplustest=Aplus(:,:,1);
	Aminustest=Aminus(:,:,1);
	for st=2:h
		t0=max(max(abs(A0test-A0(:,:,st))));
		tplus=max(max(abs(Aplustest-Aplus(:,:,st))));
		tminus=max(max(abs(Aminustest-Aminus(:,:,st))));
		tmax=max([t0,tplus,tminus]);
		if tmax>1e-9
			flag=false;
			break
		end
	end
end