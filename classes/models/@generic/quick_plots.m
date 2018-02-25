function hdl=quick_plots(m,batch,var_list,fig_title,r0c0,xrange)
% H1 line
%
% Syntax
% -------
% ::
%
%   hdl=quick_plots(m,batch)
%
%   hdl=quick_plots(m,batch,var_list)
%
%   hdl=quick_plots(m,batch,var_list,fig_title)
%
%   hdl=quick_plots(m,batch,var_list,fig_title,r0c0)
%
%   hdl=quick_plots(m,batch,var_list,fig_title,r0c0,xrange)
%
% Inputs
% -------
%
% - **m** [rise|dsge|svar|rfvar]: model object
%
% - **batch** [struct]: structure with variables to plot
%
% - **var_list** [cellstr]: list of variables of interest
%
% - **fig_title** [char|{'no title'}]: title of the figure
%
% - **r0c0** [vector|{[4,4]}]: number of rows and columns in figure
%
% - **xrange** [serial|char|{[]}]: range over which to plot
%
% Outputs
% --------
%
% - **hdl** [handle]: handle to the plotted objects
%
% See also: QUICK_IRFS

if nargin<6
    
    xrange=[];
    
    if nargin<5
        
        r0c0=[];
        
        if nargin<4
            
            fig_title='no title';
            
            if nargin<3
                
                var_list=[];
                
            end
            
        end
        
    end
    
end

if isempty(xrange)
    
    xrange={};
    
end

if ~isempty(xrange) && ~iscell(xrange)
    
    xrange={xrange};
    
end

if isempty(r0c0)
    
    r0c0=[4,4];
    
end

if isempty(var_list)
    
    var_list=fieldnames(batch);
    
end

r0=r0c0(1);

c0=r0c0(2);

description=get(m,'tex(long)');

hfig=utils.plot.multiple(@plotfunc,var_list,fig_title,r0,c0);

if nargout
    
    hdl=hfig;
    
end

    function [texname,theLegend]=plotfunc(var_name)
        
        if ischar(var_name)
            
            var_name=cellstr(var_name);
            
        end
        
        nnames=numel(var_name);
        
        theLegend=cell(1,nnames); texname=cell(1,nnames); theIrf=cell(1,nnames);
        
        for ii=1:nnames
            
            [theLegend{ii},texname{ii},theIrf{ii}]=extract_one(var_name{ii});
            
            if iscellstr(theLegend{ii})
                
                if nnames~=1
                    
                    error('multiple legends per variable and multiple variables')
                    
                end
                
                theLegend=theLegend{ii};
                
            end
            
        end
        
        if nnames>1||any(cellfun(@(x)isempty(x),theLegend,'uniformOutput',true))
            
            theLegend='';
            
        end
               
        plot(xrange{:},[theIrf{:}],'linewidth',2)
        
        function [theLegend,texname,theIrf]=extract_one(vname)
            
            bv=batch.(vname);
            
            if isstruct(bv)
                
                theLegend=bv.regime_1.varnames;
                
                theIrf=bv.regime_1;
                
                warning('Multiple regimes detected. extracting first only')
                
            else
                
                theLegend=bv.varnames;
                
                theIrf=bv;
                
            end
            
            if numel(theLegend)==1||...
                    all(cellfun(@(x)isempty(x),theLegend,'uniformOutput',true))
                
                theLegend='';
                
            end
            
            var_tex=description.(vname);
            
            texname=[var_tex,'(',vname,')'];
            
        end
        
    end

end