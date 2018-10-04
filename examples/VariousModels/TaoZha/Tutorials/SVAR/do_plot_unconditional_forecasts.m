function do_plot_unconditional_forecasts(mycast,m)

vnames=m.endogenous;

nvars=numel(vnames);

figure('name','unconditional forecasts');

for ivar=1:nvars
    
    subplot(nvars,1,ivar);
    
    plot(mycast.(vnames{ivar}),'linewidth',2);
    
    ylabel(vnames{ivar});
    
end

tmp=set([mycast.ygap mycast.pi mycast.FFR],'varnames',{'ygap','pi','FFR'});

disp('********** Unconditional forecasts ***********');
display(tmp)

end