function flag=is_eigenvalue_solver_candidate(Aplus,A0,Aminus,Q)
flag=true;
h=numel(A0);
if ~(h==1 || all(diag(Q)==1))
	A0test=A0{1};
    Aplus=reconfigure(Aplus);
	Aplustest=Aplus{1,1};
	Aminustest=Aminus{1};
	for st=2:h
		t0=get_max(A0test-A0{st});
		tplus=get_max(Aplustest-Aplus{st});
		tminus=get_max(Aminustest-Aminus{st});
		tmax=max([t0,tplus,tminus]);
		if tmax>1e-9
			flag=false;
			break
		end
	end
end

    function Apl=reconfigure(Aplus)
        Apl=cell(h,1);
        for ii=1:h
            Apl{ii}=Aplus{ii,ii}/Q(ii,ii);
        end
    end
end

function m=get_max(x)
m=max(abs(x(:)));
end