function hdl=quick_irfs(m,myirfs,var_list,shock_list,r0c0,xrange)
% QUICK_IRFS -- quick plotting of impulse response functions
%
% Syntax
% -------
% ::
%
%   hdl=quick_irfs(m,myirfs)
%
%   hdl=quick_irfs(m,myirfs,var_list)
%
%   hdl=quick_irfs(m,myirfs,var_list,shock_list)
%
%   hdl=quick_irfs(m,myirfs,var_list,shock_list,r0c0)
%
%   hdl=quick_irfs(m,myirfs,var_list,shock_list,r0c0,xrange)
%
% Inputs
% -------
%
% - **m** [rise|dsge|svar|rfvar]: model object
%
% - **myirfs** [struct]: structure with irfs
%
% - **var_list** [cellstr]: list of variables of interest
%
% - **shock_list** [cellstr]: list of shocks of interest
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: QUICK_PLOTS

if nargin<6
    
    xrange=[];
    
    if nargin<5
        
        r0c0=[];
        
        if nargin<4
            
            shock_list=[];
            
            if nargin<3
                
                var_list=[];
                
            end
            
        end
        
    end
    
end

if isempty(shock_list)
    
    shock_list=get(m,'exo_list');
    
end

description=get(m,'tex');

nshocks=numel(shock_list);

hfig=cell(1,nshocks);

for ishock=1:nshocks
    
    shock_name=shock_list{ishock};
    
    shocktex=description.(shock_name);
    
    fig_title=['Impulse responses to a ',shocktex];
    
    hfig{ishock}=quick_plots(m,myirfs.(shock_name),var_list,fig_title,r0c0,xrange);
    
end

if nargout
    
    hdl=hfig;
    
end

end