function plot_smirnov(xparam,behave,lb,ub,parameter_names,r0,c0)
%% density plots
titel='CDF behavior vs non-behavior';
if ischar(parameter_names),parameter_names=cellstr(parameter_names);end

npar=numel(parameter_names);

nstar=r0*c0;

nfig=ceil(npar/nstar);
% SaveUnderName='BehaviorVsNonBehavior';
for fig=1:nfig
    if nfig>1
        titelfig=[titel,' ',int2str(fig)];
        %         SaveUnderName0=[SaveUnderName,'_',int2str(fig)];
    else
        %         SaveUnderName0=SaveUnderName;
        titelfig=titel;
    end
    hfig=figure('name',titelfig);
    [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,npar,r0,c0);
    for ii=1:min(nstar,Remains)
        subplot(r,c,ii)
        var_id=(fig-1)*nstar+ii;
        
        x=transpose(xparam(var_id,behave));
        x_=transpose(xparam(var_id,~behave));
        [hh]=smirnov_subplot(x,x_,lb(var_id),ub(var_id),parameter_names{var_id});
        if ii==1
            leg=legend('behavior','non-behavior');
            set(leg,'Orientation','horizontal',...
                'Position',[0.4499 0.0442 0.11962 0.0219])
        end
    end
    %     saveas(hfig,[SaveUnderName0,'.pdf'])
    %     saveas(hfig,[SaveUnderName0,'.fig'])
    %     saveas(hfig,[SaveUnderName0,'.eps'])
end

function [hh]=smirnov_subplot(x,x_,lb,ub,parname)
[xx,f]=distributions.empirical_cdf(x,lb,ub);
[xx_,f_]=distributions.empirical_cdf(x_,lb,ub);
plot(xx,f,xx_,f_)
axis([lb,ub,0,1])
ff=abs(f-f_);
largest=ff==max(ff);
largest=find(largest);
if numel(largest)>1
    disp([parname,' has more than one location with largest distance'])
    largest=largest(1);
end
Dn=ff(largest);
xn=xx(largest);
hold on
plot([xn,xn],[f(largest),f_(largest)])
%----------------
n1     =  length(x);
n2     =  length(x_);
n      =  n1 * n2 /(n1 + n2);
lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * Dn , 0);
%
%  Use the asymptotic Q-function to approximate the 2-sided P-value.
%
j       =  (1:101)';
pValue  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
pValue  =  min(max(pValue, 0), 1);
%----------------
plottitle=[parname,' (Dn=',num2str(Dn),', P-value=',num2str(pValue),')'];
title(plottitle,'interpreter','none')
hh=gca;
