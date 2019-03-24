%--- help for dsge/create_estimation_blocks ---
%
%  INTERNAL FUNCTION
% 
%  ::
% 
%     blocks = create_estimation_blocks(obj,blocks)
% 
%  Args:
% 
%     obj (rise | dsge): model object
%     blocks (cell array): blocking information
% 
%  Returns:
%     :
% 
%     - **blocks**:
% 
%  Note:
% 
%     this function separates the parameters to estimate into blocks controlled
%     by one or several markov chains, depending on the information provided in
%     blocks:
% 
%        #. blocks={[1,3,5,7],[2,9,10],[6,20],...} The numbers represent the
%           order in which the estimated parameters are entered
% 
%        #. blocks={{'alpha','beta'},'gam',{'delta','upsil','omicr'},...}. In
%           this case, the first block includes alpha and beta, which are
%           either parameter names or chain names. If they are parameters,
%           then alpha and beta will be estimated simultaneously. If they are
%           chain names, then all the parameters controlled by alpha and all
%           the parameters controlled by beta will be estimated
%           simultaneously
% 
%     The estimated parameters that are controlled by a markov chain can be
%     entered either as pname(chain,state) or as pname_chain_state
% 
%