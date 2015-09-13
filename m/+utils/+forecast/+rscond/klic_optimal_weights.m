function paistar=klic_optimal_weights(gam,pai,gy)

tgam=gam(:).';
paistar=pai;
n=numel(pai);
for ii=1:n
    paistar(ii)=pai(ii)*exp(tgam*gy(:,ii));
end

paistar=paistar/sum(paistar);
end