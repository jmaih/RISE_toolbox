function out=fanchart(this,ci)
% Creates data for fanchart
%
% ::
%
%    out=fanchart(this,ci)
%
% Args:
%
%    this (ts): time series with either multiple pages or multiple
%      columns, but not both
%
%    ci (vector): numbers in [0,1] or in [0,100] to be used in the
%      calculation of the width of the fans
%
% Returns:
%    :
%
%    - **out** [struct]: time series object
%
%       - **ci** [vector]: same as above
%       - **median** [vector]: median of the data
%       - **variance** [vector]: variance of the data
%       - **quantiles** [matrix]: quantiles defined by **ci**
%       - **prob_index** [vector]: locations of the cutoffs in the data
%       - **probs** [vector]: probabilities for the quantiles (function of
%         ci)
%       - **date_numbers** [vector]: serial dates for reconstructing time
%         series
%
% Example:
%    ::
%
%           this=ts('1990Q2',rand(100,1000));
%           out=fanchart(this,[30,50,70,90])
%           out=fanchart(this,[30,50,70,90]/100)
%           plot_fanchart(out,'r',10)
%
% Note:
%    If the input time series contains several variables, the output is a
%    structure with the names of the different variables in the first level of
%    the fields
%
% See also:
%    - :func:`plot_fanchart <ts.plot_fanchart>`
%

out=struct();

nvar=this.NumberOfVariables;

if nvar>1

    vnames=this.varnames;

    if ~isempty(vnames{1})

        for ivar=1:nvar

            S=struct('type','()','subs',{vnames(ivar)});

            obj=subsref(this,S);%<--obj=this(vnames{ivar});

            out.(vnames{ivar})=fanchart(obj,ci);

        end

        return

    end

end

if any(ci<0)

    error('confidence bands cannot be negative')

end

ci=sort(ci(:));

large=ci>=1;

ci(large)=ci(large)/100;

datax=squeeze(double(this));

if size(datax,3)>1

    error('dataset must have multiple columns or multiple pages but not both')

end

probs = 0.5*(1+[-ci,ci]);

probs=[flipud(probs(:,1));probs(:,2)];

emp_moms=distributions.empirical_moments(datax,[],[],probs);

out=struct();

out.ci=ci;

out.mean=[emp_moms.mean];

out.median=[emp_moms.median];

out.variance=[emp_moms.variance];

out.quantiles=vertcat(emp_moms.quantiles);

out.prob_index=emp_moms(1).prob_index;

out.probs=emp_moms(1).probs;

out.date_numbers=this.date_numbers;

end