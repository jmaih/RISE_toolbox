function fig_handles=plot_separate(this,varargin)

vnames=this.varnames;
nvars=this.NumberOfVariables;
if ~(nvars==numel(vnames))
    error([mfilename,':: each variable should have a name'])
end
r0=3;
c0=3;
nstar=c0*r0;
nfig=ceil(nvars/nstar);
fig_handles=nan(nfig,1);
for fig=1:nfig
    fig_handles(fig)=figure('name',['multiple time series plots # ',int2str(fig)]);
    [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,nvars,r0,c0);
    for ii=1:min(nstar,Remains)
        id=(fig-1)*nstar+ii;
        S=struct('type','()','subs',{vnames(id)});
        obj=this.subsref(S); % <--- obj=this(vnames{id}); does not work
        subplot(r,c,ii)
        plot(obj,varargin{:});
        title(vnames{id},'interpreter','none')
    end
end

