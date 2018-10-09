function [sep_cf,cf]=complementarity_memoizer(obj)
% Memoize linear and nonlinear restrictions
%
% ::
%
%   [sep_cf,cf]=complementarity_memoizer(obj)
%
% Args:
%
%    obj (rise | dsge | svar | rfvar): model object
%
% Returns:
%    :
%
%    - **sep_cf** [function handle]: function that takes as input a vector of
%      variable values and returns a vector of values that are expected to be
%      greater than or equal to 0.
%
%    - **cf** [function handle]: function that takes as input a vector of
%      variable values and returns a true if all the restrictions are satisfied
%      and returns false otherwise
%

if ~isfield(obj.routines,'complementarity')||... I do not expect VARs to have this although in theory they could.
        isempty(obj.routines.complementarity)

    routines=[];

else

    routines=obj.routines.complementarity;

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

sep_cf=[];

cf=[];

if ~hc

    return

end

dd=cell(1,nrout);

for ii=1:nrout

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

    dd{ii}=strrep(rout_i(right_par+1:end),';','');

end

dd=cell2mat(strcat(dd,';'));

dd=dd(1:end-1);

is_log_var=obj.log_vars;

additional_inputs=load_inputs(obj,rise_inputs);

[sep_cf, cf]=do_memo(dd, is_log_var,additional_inputs{:});

end

function ainp=load_inputs(obj,rise_inputs)

% inputs and missing inputs
%---------------------------
def=obj.solution.definitions{1}; %#ok<NASGU>

ss=obj.solution.ss{1}; %#ok<NASGU>

param=obj.parameter_values(:,1); %#ok<NASGU>

x=[];s0=1;s1=1; %#ok<NASGU>

% load the additional inputs
%----------------------------
ninp=numel(rise_inputs);

ainp=cell(1,ninp-1);

for ii=2:ninp

    ainp{ii-1}=eval(rise_inputs{ii});

end

end

function [sep_cf, cf]=do_memo(dd, is_log_var,varargin)

args=cell2mat(strcat(parser.input_list,','));

sep_cf0=str2func(['@(',args(1:end-1),')[',dd,']']);

cf0=str2func(['@(',args(1:end-1),')all([',dd,'])']);

sep_cf=@engine_sepcf;

cf=@engine_cf;

    function out=engine_sepcf(y)

        y(is_log_var)=exp(y(is_log_var));

        out=sep_cf0(y,varargin{:});

    end


    function out=engine_cf(y)

        y(is_log_var)=exp(y(is_log_var));

        out=cf0(y,varargin{:});

    end

end