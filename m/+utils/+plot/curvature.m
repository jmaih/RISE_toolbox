function [h,legend_,tex_name]=curvature(xname,db)
% curvature -- plots the curvature of an estimated parameter at the model
%
% ::
%
%
%   [h,legend_]=curvature(xname,db)
%
% Args:
%
%    - **xname** [char]: name of the parameter to plot
%
%    - **db** [struct]: structure containing the various parameters. Each
%    parameter field is itself a structure with the following
%      - **tex_name** [char]: name of the parameter as it should appear in the
%      title of the plot
%      - **mode** [scalar]: value of the parameter at the mode
%      - **log_post_mode** [scalar]: value of the log posterior mode
%      - **log_lik_mode** [scalar]: value of the log likelihood at the mode
%      - **x** [vector]: x-axis values
%      - **log_post** [vector]: value of the log-posterior for each value of x
%      - **log_lik** [vector]: value of the log-likelihood for each value of x
%
% Returns:
%    :
%
%    - **h** [handle]: handle for the plot
%
%    - **legend_** [cellstr]: names of the lines in the plot
%
%    - **tex_name** [char]: name of the parameter
%
% Note:
%
% Example:
%
%    See also:


legend_={'log post','log lik','mode'};
pp=db.(xname);
tex_name=pp.tex_name;

low_f=min(min([pp.log_post,pp.log_lik]));
high_f=max(max([pp.log_post,pp.log_lik]));
rescale=(high_f-low_f)/10000;

posj=find(abs(pp.x-pp.mode)==min(abs(pp.x-pp.mode)),1,'first');
plot(pp.x,pp.log_post,...
    pp.x,pp.log_lik,...
    [pp.x(posj),pp.x(posj)],[(1-rescale)*low_f,high_f],...
    'linewidth',1.5)
title(pp.tex_name)

if any(isnan(pp.log_post))
    hold on
    locs=find(isnan(pp.log_post));
    plot(pp.x(locs)',min(pp.log_post)*ones(numel(locs),1),'.r','markersize',10)
    hold off
end

% size of the plot
%-----------------
x_min=min(pp.x);

x_max=max(pp.x);

llf=(1-sign(low_f)*rescale)*low_f;

hhf=(1+sign(high_f)*rescale)*high_f;

if abs(hhf-llf)<eps
    
    hhf=inf;
    
end

axis([x_min x_max llf hhf]);
% % axis tight % xlim([min(pp.x),max(pp.x)])
h=gca();

end
