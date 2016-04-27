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

if isempty(r0c0)
    
    r0c0=[4,4];
    
end

if isempty(var_list)
    
    var_list=fieldnames(batch);
    
end

r0=r0c0(1);

c0=r0c0(2);

description=get(m,'tex');

hfig=utils.plot.multiple(@plotfunc,var_list,fig_title,r0,c0);

if nargout
    
    hdl=hfig;
    
end

    function [texname,theLegend]=plotfunc(var_name)
        
        theLegend='';
        
        var_tex=description.(var_name);
        
        texname=[var_tex,'(',var_name,')'];
        
        theIrf=batch.(var_name);
        
        if isempty(xrange)
            
            plot(theIrf,'linewidth',2)
            
        else
            
            plot(xrange,theIrf,'linewidth',2)
            
        end
        
    end

end