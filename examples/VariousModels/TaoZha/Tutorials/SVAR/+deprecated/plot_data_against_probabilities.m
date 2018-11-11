function plot_data_against_probabilities(sv,type,specs)

if nargin<3
    
    specs=[3,3];
    
    if nargin<2
        
        type='state';
        
    end
    
end

% myLimits=[-sqrt(eps),1+sqrt(eps)];

[f,the_regimes,endog,data,tex]=load_filters(sv);

switch type
    
    case 'state'
        
        probs=f.smoothed_state_probabilities;
        
    case 'regime'
        
        probs=f.smoothed_regime_probabilities;
        
    otherwise
        
    error('second input must be either "state" or "regime"')
    
end

r0=specs(1);

c0=specs(2);

all_regimes=fieldnames(probs);

all_regimes(strcmp(all_regimes,'const_1'))=[];

for iv=1:numel(endog)
    
    vn=endog{iv};
    
    utils.plot.multiple(@plotfuncr,all_regimes,...
        [vn,' against smoothed ',type,' probabilities'],r0,c0);
    
    if ~isempty(tex)
        
        [~,h]=sup_label(tex.(vn),'t');
        
        set(h,'fontsize',12)
        
    end

end


    function [mytex,leg]=plotfuncr(vname)
        
        dd=data.(vn)(probs.(vname).date_numbers);
        
        [AX,H1,H2]=plotyy(probs.(vname),dd,'linewidth',2);
        
        axis(AX,'tight')
        
        if isempty(tex)
            
            mytex=vname;
            
        else
            
            mytex=tex.(vname);
            
        end
        
        if strcmp(type,'regime')
            
            vname(1:numel('regime'))=[];
            
            vname=strrep(vname,'_','');
            
            mytex=sprintf('%s(%s)',mytex,the_regimes{str2double(vname)});
            
        end
        
        leg='';
        
    end

end