function plot_probabilities(sv,specs)

if nargin<2
    
    specs=[3,3];
    
end

myLimits=[-sqrt(eps),1+sqrt(eps)];

r0=specs(1);

c0=specs(2);

[f,the_regimes]=load_filters(sv);

all_regimes=fieldnames(f.smoothed_regime_probabilities);

all_states=fieldnames(f.smoothed_state_probabilities);

utils.plot.multiple(@plotfuncr,all_regimes,...
    'smoothed regime probabilities',r0,c0);

utils.plot.multiple(@plotfuncs,all_states,...
    'smoothed state probabilities',r0,c0);


    function [tex,leg]=plotfuncr(vname)
        
        plot(f.smoothed_regime_probabilities.(vname),'linewidth',2)
        
        ylim(myLimits)
        
        tex=vname;
        
        vname(1:numel('regime'))=[];
        
        vname=strrep(vname,'_','');
        
        tex=sprintf('%s(%s)',tex,the_regimes{str2double(vname)});
        
        leg='';
        
    end

    function [tex,leg]=plotfuncs(vname)
        
        plot(f.smoothed_state_probabilities.(vname),'linewidth',2)
        
        ylim(myLimits)
        
        tex=vname;
        
        leg='';
        
    end

end