function hh=plot_correlation_patterns(corrmat,parnames,cutoff,type)

if nargin<4
    type='undeclared';
    if nargin<3
        cutoff=.5;
    end
end

if ischar(parnames),parnames=cellstr(parnames);end

nparam=numel(parnames);
[rr,cc]=size(corrmat);
if rr~=cc
    if rr>cc && cc==nparam
        corrmat=corr(corrmat);
    elseif cc>rr && rr==nparam
        corrmat=corr(corrmat');
    else
        error([mfilename,':: maxtrix inconsistent with # of parameters'])
    end
end
pax=nan(nparam,nparam);
corrmat=tril(corrmat,-1);
hotties=abs(corrmat)>abs(cutoff);
pax(hotties)=corrmat(hotties);

titel=['Correlation patterns in the ',type,' sample '];
hh=figure('name',titel);

imagesc(pax,[-1 1]);

for ip=1:nparam
    text(ip,(0.5),parnames{ip},'HorizontalAlignment','left',...
        'rotation',90,'interpreter','none')
    text(0.5,ip,parnames{ip},'HorizontalAlignment','right',...
        'interpreter','none')
end

colorbar;
ax=colormap;
ax(1,:)=[0.9 0.9 0.9];
colormap(ax);

if nparam>10
    set(gca,'xtick',(5:5:nparam))
    set(gca,'ytick',(5:5:nparam))
end

set(gca,'xgrid','on')
set(gca,'ygrid','on')

xlabel(titel)
