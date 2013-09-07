function ysim=simulate_second_order(solution,Qfunc,pai0,y0,sig,nsteps)
% Q=solution.Q;
ss=solution.ss;
bgp=solution.bgp;
mx=solution.mx;
me=solution.me;
msig=solution.msig;
mxx=solution.mxx;
mee=solution.mee;
mxe=solution.mxe;
mxsig=solution.mxsig;
mesig=solution.mesig;
msigsig=solution.msigsig;
clear solution

[ny,ne]=size(me{1});
ysim=zeros(ny,nsteps+1);
ysim(:,1)=y0;
pai=update_probability(y0,x,param,ss,def,pai0);
for t=2:nsteps+1
    % draw state
    cpai=cumsum(pai);
    cpai=[0,cpai(:)'];
    cpai(end)=1;
    st=find(rand<=cpai,1,'first')-1;
    % draw shock
    et=randn(ne,1);
    ysim(:,t)=iterate_simulation(ysim(:,t-1));
pai=update_probability(y0,x,param,ss,def,pai);
end

    function y1=iterate_simulation(y0)
        xhat=y0-ss(:,st);
        y1=ss(:,st)+(0.5*sig*mxsig{st}+mx{st})*xhat+...
            (0.5*sig^2*msigsig{st}+sig*msig{st})+...
            0.5*mxx{st}*kron(xhat,xhat);
        if ~isempty(et)
            y1=y1+...
            (0.5*sig*mesig{st}+me{st})*et+...
            0.5*mee{st}*kron(et,et)+...
            0.5*mxe{st}*kron(xhat,et);
        end
    end

    function pai=update_probability(y,x,param,ss,def,pai0)
        Q=Qfunc(y,param(:,1),x,ss,def);
        pai=Q'*pai0;
    end
end
