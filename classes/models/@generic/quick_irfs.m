function hdl=quick_irfs(m,myirfs,var_list,shock_list,r0c0,xrange,suplab)
% Quick plotting of impulse response functions
%
% ::
%
%
%   hdl=quick_irfs(m,myirfs)
%   hdl=quick_irfs(m,myirfs,var_list)
%   hdl=quick_irfs(m,myirfs,var_list,shock_list)
%   hdl=quick_irfs(m,myirfs,var_list,shock_list,r0c0)
%   hdl=quick_irfs(m,myirfs,var_list,shock_list,r0c0,xrange)
%   hdl=quick_irfs(m,myirfs,var_list,shock_list,r0c0,xrange,suplab)
%
% Args:
%
%    m (rise | dsge | svar | rfvar): model object
%
%    myirfs (struct): structure with irfs
%
%    var_list (cellstr): list of variables of interest
%
%    shock_list (cellstr): list of shocks of interest
%
%    r0c0 (vector | {[4,4]}): number of rows and columns in figure
%
%    xrange (serial | char | {[]}): range over which to plot
%
%    suplab (false | {true}): add a sup-label to the figure or not
%
% Returns:
%    :
%
%    - **hdl** [handle]: handle to the plotted objects
%
% See also: quick_plots
%

if nargin<7
    
    suplab=[];
    
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
    
end

if isempty(suplab)
    
    suplab=true;
    
end

if isempty(shock_list)
    
    shock_list=get(m(1),'exo_list');
    
end

if ischar(shock_list)
    
    shock_list=cellstr(shock_list);
    
end

if isempty(var_list)
    
    var_list=get(m(1),'endo_list');
    
end

description=get(m(1),'tex(long)');

nshocks=numel(shock_list);

for ishock=1:nshocks
    
    shock_name=shock_list{ishock};
    
    shocktex=description.(shock_name);
    
    irf_batch=myirfs.(shock_name);
    
    v1=var_list{1};
    
    if iscell(v1)
        
        v1=v1{1};
        
    end
    
    is_multiple_regimes=isstruct(irf_batch.(v1));
    
    if is_multiple_regimes
        
        regime_names=fieldnames(irf_batch.(v1));
        
        nregs=numel(regime_names);
        
        tmp=repmat({struct()},1,nregs);
        
        for ireg=1:nregs
            
            for ii=1:numel(var_list)
                
                v=var_list{ii};
                
                if iscell(v)
                    
                    error('multiple variables in one plot not allowed in multiple regimes')
                    
                end
                
                tmp{ireg}.(v)=irf_batch.(v).(regime_names{ireg});
                
            end
            
            regime_names{ireg}=['(',regime_names{ireg},')'];
            
        end
        
        irf_batch=tmp;
        
    else
        
        nregs=1;
        
        regime_names={''};
        
        irf_batch={irf_batch};
        
    end
    
    if ishock==1
        
        hfig=cell(nshocks,nregs);
        
    end
    
    for ireg=1:nregs
        
        fig_title=['Impulse responses to a ',shocktex,regime_names{ireg}];
        
        hfig{ishock,ireg}=quick_plots(m,irf_batch{ireg},var_list,fig_title,r0c0,xrange);
        
        if suplab
            
            [~,h]=sup_label(fig_title,'t');
            
            set(h,'fontsize',12)
            
        end
        
    end
    
end

hfig=[hfig{:}];

if nargout
    
    hdl=hfig;
    
end

end