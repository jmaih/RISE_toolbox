function objective=klic_lagrange_objective(gam,pai,gy_gbar)
tgam=gam(:).';
objective=0;
n=numel(pai);
for ii=1:n
    objective=objective+pai(ii)*exp(tgam*gy_gbar(:,ii));
end
end
