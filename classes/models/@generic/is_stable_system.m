%--- help for generic/is_stable_system ---
%
%  Checks the stability of a linear markov switching system.
% 
%  ::
% 
%     flag = is_stable_system(obj)
%     flag = is_stable_system(obj,varargin)
% 
%  Args:
% 
%     obj (dsge | rise | svar | rfvar): model object
% 
%     varargin (name,value): pairwise valid options for RISE. The most
%       relevant in this case are
% 
%        - **stability_criterion** [numeric\|{1.000001}]: stability criterion.
%          All eigenvalues must be smaller than this criterion for the system to
%          be MSS
%        - **stability_algorithm** ['cfm'\|{'hmg'}]: CFM stands for
%          Costa-Fragoso-Marques while HMG stands for Hassibi-Murray-Gupta.
% 
%  Returns:
%     :
% 
%     - **flag** [false\|true]: result of the investigation on whether the
%       system is stable or not.
% 
%  Note:
% 
%     RISE implements two algorithms from the engineering literature to check
%     for the stability. They are
% 
%        - Costa-Fragoso-Marques :cite:`costa2006discrete`
%        - Hassibi-Murray-Gupta :cite:`gupta2003control`
% 
%     Refer to the references to see the specific algorithms. However, for most
%     applications, one can just use the default options.
% 
%  References:
% 
%     - :cite:`costa2006discrete`
%     - :cite:`gupta2003control`
% 
%