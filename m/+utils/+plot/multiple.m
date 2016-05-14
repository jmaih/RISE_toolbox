function hfig=multiple(plotfunc,vnames,fig_title,r0,c0,varargin)
% multiple -- multiple plots in graphs
%
% Syntax
% -------
% ::
%
%   hfig=multiple(plotfunc,vnames,fig_title,r0,c0)
%
%   hfig=multiple(plotfunc,vnames,fig_title,r0,c0,varargin)
%
% Inputs
% -------
%
% - **plotfunc** [function handle]: which takes as input a valid variable
%   name and returns (1) the name/description of the variable/parameter
%   plotted, (2) the legend
%
% - **vnames** [cellstr]: names of variables to be plotted
%
% - **fig_title** [char]: main title of the figures to plot
%
% - **r0** [integer]:: the desired maximum number of rows in each figure
%
% - **c0** [integer]:: the desired maximum number of columns in each figure
%
% - **varargin** [|{}]:: pairwise elements of entering the title and the
% legend output arguments
%
% Outputs
% --------
%
% - **hfig** [handles]: handles to the different figures created
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

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
            
            if ~isempty(varargin)
                
                set(hleg,varargin{:});
                
            end
            
        end
        
    end
    
end

end
