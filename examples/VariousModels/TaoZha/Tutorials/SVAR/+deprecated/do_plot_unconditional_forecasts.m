function do_plot_unconditional_forecasts(mycast,m)

vnames=m.endogenous.name;
vnames_tex=m.endogenous.tex_name;
nvars=numel(vnames);

figure('name','unconditional forecasts');
for ivar=1:nvars
    subplot(nvars,1,ivar);
    plot(mycast.(vnames{ivar}),'linewidth',2);
    ylabel(vnames_tex{ivar});
end

disp('********** Unconditional forecasts ***********');
display([mycast.ygap mycast.pi mycast.FFR])

end