%--- help for ts/fanchart ---
%
%  Creates data for fanchart
% 
%  ::
% 
%     out=fanchart(this,ci)
% 
%  Args:
% 
%     this (ts): time series with either multiple pages or multiple
%       columns, but not both
% 
%     ci (vector): numbers in [0,1] or in [0,100] to be used in the
%       calculation of the width of the fans
% 
%  Returns:
%     :
% 
%     - **out** [struct]: time series object
% 
%        - **ci** [vector]: same as above
%        - **median** [vector]: median of the data
%        - **variance** [vector]: variance of the data
%        - **quantiles** [matrix]: quantiles defined by **ci**
%        - **prob_index** [vector]: locations of the cutoffs in the data
%        - **probs** [vector]: probabilities for the quantiles (function of
%          ci)
%        - **date_numbers** [vector]: serial dates for reconstructing time
%          series
% 
%  Example:
%     ::
% 
%            this=ts('1990Q2',rand(100,1000));
%            out=fanchart(this,[30,50,70,90])
%            out=fanchart(this,[30,50,70,90]/100)
%            plot_fanchart(out,'r',10)
% 
%  Note:
%     If the input time series contains several variables, the output is a
%     structure with the names of the different variables in the first level of
%     the fields
% 
%  See also:
%     - :func:`plot_fanchart <ts.plot_fanchart>`
% 
%