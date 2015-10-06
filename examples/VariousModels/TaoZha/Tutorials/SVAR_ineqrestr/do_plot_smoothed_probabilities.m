function [hfig_st,hfig_reg]=do_plot_smoothed_probabilities(m)

plot_regimes_also = m.markov_chains.chains_number>2;

state_names=m.markov_chains.state_names-'const_1';

smooth_state_probs=m.filtering.smoothed_state_probabilities;

multiple=@utils.plot.multiple;

fig_title='Smoothed State probabilities';
nrows=3;
ncols=3;

hfig=multiple(@plot_state,state_names,fig_title,nrows,ncols);
if nargout
    hfig_st=hfig;
    hfig_reg=[];
end
if plot_regimes_also
    fig_title='Smoothed Regime probabilities';
    smooth_regime_probs=m.filtering.smoothed_regime_probabilities;
    regime_names=fieldnames(smooth_regime_probs);
    hfig=multiple(@plot_regime,regime_names,fig_title,nrows,ncols);
    if nargout
        hfig_reg=hfig;
    end
end

    function [tex_name,legend_]=plot_regime(name)
        tex_name=name;
        legend_=[];
        plot(smooth_regime_probs.(name),'linewidth',2)
    end

    function [tex_name,legend_]=plot_state(name)
        tex_name=name;
        legend_=[];
        plot(smooth_state_probs.(name),'linewidth',2)
    end
end