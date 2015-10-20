function [sep_cf,cf]=complementarity_memoizer(obj)
% complementarity_memoizer - memoizes linear and nonlinear restrictions
%
% Syntax
% -------
% ::
%
%   [sep_cf,cf]=complementarity_memoizer(obj)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|svar|rfvar]: model object
%
% Outputs
% --------
%
% - **sep_cf** [function handle]: function that takes as input a vector of
% variable values and returns a vector of values that are expected to be
% greater than or equal to 0.
%
% - **cf** [function handle]: function that takes as input a vector of
% variable values and returns a true if all the restrictions are satisfied
% and returns false otherwise
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if ~isfield(obj.routines,'complementarity')||... I do not expect VARs to have this although in theory they could.
        isempty(obj.routines.complementarity)
    routines=[];
else
    routines=obj.routines.complementarity;
    def=obj.solution.definitions{1}; %#ok<NASGU>
    ss=obj.solution.ss{1}; %#ok<NASGU>
    param=obj.parameter_values(:,1); %#ok<NASGU>
    % missing inputs
    %-----------------
    x=[];sparam=[];s0=1;s1=1; %#ok<NASGU>
    rise_inputs=parser.input_list;
end
nrout=numel(routines);
hc=nrout>0 && obj.options.simul_honor_constraints;
simul_honor_constraints_through_switch=obj.markov_chains.regimes_number>1 &&...
    obj.options.simul_honor_constraints_through_switch;
if hc && ~(simul_honor_constraints_through_switch||...
        any(obj.exogenous.shock_horizon(:)>0)||...
        ~isempty(obj.options.solve_occbin))
    error('restrictions detected but no anticipatory or switching behavior to satisfy them')
end
clear obj

sep_cf=[];
cf=[];
if ~hc
    return
end
dd=cell(1,nrout);
for ii=1:nrout
    if ii==1
        % make sure the inputs [EXCEPT THE FIRST] are all here
        %------------------------------------------------------
        for inp=2:numel(rise_inputs)
            eval([rise_inputs{inp},'=',rise_inputs{inp},';'])
        end
    end
    % those functions are expected to return scalar values
    %-----------------------------------------------------
    rout_i=func2str(routines{ii});
    right_par=find(rout_i==')',1,'first');
    if ii==1
        header=rout_i(1:right_par);
    end
    % check that all inputs are correct
    %-----------------------------------
    your_inputs=regexp(header(3:end-1),',','split');
    if ~all(strcmp(your_inputs,rise_inputs))
        error('wrong input names')
    end
    ge_=strfind(rout_i,'>=');
    if isempty(ge_)||numel(ge_)>1
        error('">=" missing or too many ">="')
    end
    dd{ii}=rout_i(right_par+1:ge_-1);
end
dd=cell2mat(strcat(dd,';'));

sep_cf=str2func(['@(',rise_inputs{1},')[',dd(1:end-1),']']);
cf=str2func(['@(',rise_inputs{1},')all([',dd(1:end-1),']>=0)']);

% no need for memoization at this point. all params, defs, ss, x, are
% indexed as if they were matrices and str2func in that case treats them as
% variables to be included in the workspace.

end