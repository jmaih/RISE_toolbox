function hdl=quick_irfs(m,myirfs,var_list,shock_list,r0c0,xrange)
% QUICK_IRFS -- quick plotting of impulse response functions
%
% ::
%
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
% Args:
%
%    - **m** [rise|dsge|svar|rfvar]: model object
%
%    - **myirfs** [struct]: structure with irfs
%
%    - **var_list** [cellstr]: list of variables of interest
%
%    - **shock_list** [cellstr]: list of shocks of interest
%
%    - **r0c0** [vector|{[4,4]}]: number of rows and columns in figure
%
%    - **xrange** [serial|char|{[]}]: range over which to plot
%
% Returns:
%    :
%
%    - **hdl** [handle]: handle to the plotted objects
%
% Note:
%
% Example:
%
%    See also: QUICK_PLOTS

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
    
    shock_list=get(m(1),'exo_list');
    
end

if ischar(shock_list)
    
    shock_list=cellstr(shock_list);
    
end

description=get(m(1),'tex(long)');

nshocks=numel(shock_list);

hfig=cell(1,nshocks);

for ishock=1:nshocks
    
    shock_name=shock_list{ishock};
    
    shocktex=description.(shock_name);
    
    irf_batch=myirfs.(shock_name);
    
    regime_names=fieldnames(irf_batch);
    
    test=regexp(regime_names,'regime_\d+','match');
    
    test=[test{:}];
    
    is_multiple_regimes=~isempty(test) && ...
        all(cellfun(@(x)~isempty(x),test,'uniformOutput',true));
    
    if is_multiple_regimes
        
        nregs=numel(regime_names);
        
        tmp=cell(1,nregs);
        
        for ireg=1:nregs
            
            tmp{ireg}=irf_batch.(regime_names{ireg});
            
            regime_names{ireg}=['(',regime_names{ireg},')'];
            
        end
        
        irf_batch=tmp;
        
    else
        
        nregs=1;
        
        regime_names={''};
        
        irf_batch={irf_batch};
        
    end
    
    for ireg=1:nregs
        
        fig_title=['Impulse responses to a ',shocktex,regime_names{ireg}];
        
        hfig{ishock,ireg}=quick_plots(m,irf_batch{ireg},var_list,fig_title,r0c0,xrange);
    end
    
end

if nargout
    
    hdl=hfig;
    
end

end