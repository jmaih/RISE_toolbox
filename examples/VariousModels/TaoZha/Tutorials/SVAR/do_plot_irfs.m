function do_plot_irfs(myirfs,m,shocks)
if nargin<3
    shocks=[];
end

shock_names=fieldnames(myirfs);
vnames=m.endogenous.name;
vnames_tex=m.endogenous.tex_name;
nvars=numel(vnames);
if isempty(shocks)
    shocks=shock_names;
end
for ishock=1:numel(shocks)
    this_shock=shocks{ishock};
    figure('name',['Impulse responses to a ',this_shock,' shock']);
    for ivar=1:nvars
        thisvar=vnames{ivar};
        subplot(nvars,1,ivar)
        thedata=myirfs.(this_shock).(thisvar);
        plot(thedata,'linewidth',2);
        title(vnames_tex{ivar});
        if strcmp(this_shock,'EXO_FFR') && strcmp(thisvar,'pi')
            disp('********** Impulse responses of inflation to an interest rate shock ***********');
            myirfs.EXO_FFR.pi
        end
        if ivar==1
            thelegend=get(thedata,'varnames');
            if ~isempty(thelegend{1})
                legend(thelegend)
            end
        end
    end
end

end