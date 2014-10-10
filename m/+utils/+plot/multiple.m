function hfig=multiple(plotfunc,vnames,fig_title,r0,c0,varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

% this function plots multiple graphs and figures
% input arguments
%   - plotfunc: a function handle which takes as input a valid variable
%   name and returns (1) the name/description of the variable/parameter
%   plotted, (2) the legend 
%   in case there are multiple plots.
%   - vnames: cellstr of variables to be plotted
%   - fig_title: the main title of the figures
%   - r0: the desired maximum number of rows in each figure
%   - c0: the desired maximum number of columns in each figure
%   - varargin: pairwise elements of entering the title and the legend
% output arguments
%   - hfig: handles to the different figures created

if ischar(vnames)
    vnames=cellstr(vnames);
end
npar=numel(vnames);
rmdoll=@(x)strrep(x,'$','');
nstar=r0*c0;
nfig=ceil(npar/nstar);
hfig=nan(nfig,1);
for fig=1:nfig
    [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,npar,r0,c0);
    if nfig>1
        titelfig=[fig_title,' ',int2str(fig)];
    else
        titelfig=fig_title;
    end
    hfig(fig)=figure('name',titelfig);
    for plt=1:min(nstar,Remains)
        par_id=(fig-1)*nstar+plt;
        figure(hfig(fig))
        subplot(r,c,plt)
        [tex_name,legend_]=plotfunc(vnames{par_id});
        if ~isempty(tex_name)
            title(rmdoll(tex_name),varargin{:})
        end
        if plt==1 && ~isempty(legend_)
            hleg=legend(legend_);
            set(hleg,varargin{:})
        end
    end
end

end
