function plot_probabilities(sv,specs)

if nargin<2
    
    specs=[3,3];
    
end

myLimits=[-sqrt(eps),1+sqrt(eps)];

r0=specs(1);

c0=specs(2);

[f,the_regimes]=load_filters(sv);

cellform = iscell(f);

if cellform
    
    all_regimes=fieldnames(f{1}.smoothed_regime_probabilities);
    
    all_states=fieldnames(f{1}.smoothed_state_probabilities);
    
    regs=the_regimes{1};
    
    for ii=2:numel(f)
        
        [all_regimes,IA,IB]=intersect(all_regimes,...
            fieldnames(f{ii}.smoothed_regime_probabilities)); %#ok<ASGLU>
        
        regs=regs(IA);
        
        [all_states]=intersect(all_states,...
            fieldnames(f{ii}.smoothed_state_probabilities));
        
    end
    
    the_regimes=regs;
    
else
    
    all_regimes=fieldnames(f.smoothed_regime_probabilities);
    
    all_states=fieldnames(f.smoothed_state_probabilities);
    
end

utils.plot.multiple(@plotfuncr,all_regimes,...
    'smoothed regime probabilities',r0,c0);

utils.plot.multiple(@plotfuncs,all_states,...
    'smoothed state probabilities',r0,c0);

    function [tex,leg]=plotfuncr(vname)
        
        d=load_item('smoothed_regime_probabilities',vname);
        
        plot(d,'linewidth',2)
        
        ylim(myLimits)
        
        tex=vname;
        
        vname(1:numel('regime'))=[];
        
        vname=strrep(vname,'_','');
        
        tex=sprintf('%s(%s)',tex,the_regimes{str2double(vname)});
        
        leg='';
        
    end

    function [tex,leg]=plotfuncs(vname)
        
        d=load_item('smoothed_state_probabilities',vname);
        
        plot(d,'linewidth',2)
        
        ylim(myLimits)
        
        tex=vname;
        
        leg='';
        
    end

    function d=load_item(a,b)
        
        if cellform
            
            d=f{1}.(a).(b);
            
            for jj=2:numel(f)
                
                d=[d,f{jj}.(a).(b)]; %#ok<AGROW>
                
            end
            
        else
            
            d=f.(a).(b);
            
        end
        
    end

end