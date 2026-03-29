%--- help for generic/is_stable_system ---
%
%  Checks the stability of a linear markov switching system.
% 
%  Syntax ::
% 
%     flag = is_stable_system(obj)
% 
%     flag = is_stable_system(obj,varargin)
% 
%  Args:
% 
%     - obj (dsge|rise): model object
% 
%     - varargin (name,value): pairwise valid options for RISE. The most
%       relevant in this case are
% 
%        - **stability_criterion** (numeric\|{1.000001}): stability criterion.
%          All eigenvalues must be smaller than this criterion for the system to
%          be MSS
% 
%        - **stability_algorithm** ('cfm'\|{'gmh'}): CFM stands for
%          Costa-Fragoso-Marques while HMG stands for Hassibi-Murray-Gupta.
% 
%  Returns:
% 
%     - **flag** (false\|true|vector|cell array): result of the
%       investigation on whether the system is stable or not.
% 
%     - **retcode** (0|25|vector|cell array): result of the investigation on
%       whether the system is stable or not.
% 
%     - **overall_retcode** (0|25|vector): result of the investigation on whether
%       the system is stable or not.
% 
%     - **cellFlag**  (cell array): stability flag for each parameterization
% 
%  Note:
% 
%     RISE implements two algorithms from the engineering literature to check
%     for the stability. They are
% 
%        - Costa-Fragoso-Marques :cite:`CostaFM2005`
% 
%        - Hassibi-Murray-Gupta :cite:`GuptaMH2003`
% 
%     Refer to the references to see the specific algorithms. However, for most
%     applications, one can just use the default options.
% 
%  References:
% 
%     - :cite:`CostaFM2005`
%     - :cite:`GuptaMH2003`
%