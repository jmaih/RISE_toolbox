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

    rise_inputs=rdico.input_list;

end

nrout=numel(routines);

hc=nrout>0;

sep_cf=[];

cf=[];

if ~hc

    return

end

dd=cell(1,nrout);

tol=obj.options.simul_tol;

for ii=1:nrout

    % those functions are expected to return scalar values
    %-----------------------------------------------------
    rout_i=func2str(routines{ii});

    right_par=find(rout_i==')',1,'first');

    if ii==1

        header=rout_i(1:right_par);

    end

    rest=rout_i(right_par+1:end);

    rest=if_elseif(any(rest=='>'),['tol+',rest],...
        any(rest=='<'),[rest,'+tol']);

    % check that all inputs are correct
    %-----------------------------------
    your_inputs=regexp(header(3:end-1),',','split');

    if ~all(strcmp(your_inputs,rise_inputs))

        error('wrong input names')

    end

    dd{ii}=strrep(rest,';','');

end

dd=cell2mat(strcat(dd,';'));

dd=dd(1:end-1);

additional_inputs=load_inputs(obj,rise_inputs);

additional_inputs=[additional_inputs,{tol}];

[sep_cf, cf]=do_memo(dd,additional_inputs{:});

% N.B: log vars are automatically transformed during shadowization
%------------------------------------------------------------------

end


function ainp=load_inputs(obj,rise_inputs)

% inputs and missing inputs
%---------------------------
c=obj.constants.values; %#ok<NASGU>

s0=1;s1=1; %#ok<NASGU>

% the discretionary-policy inputs: rdico.input_list() ends with D and L, so
% they have to exist before the eval loop below. Both come back empty for a
% model that is not an optimal-policy model.
[D,L]=dsge_tools.optimal_policy.discretion.derivatives_locators(obj); %#ok<ASGLU>

% load the additional inputs
%----------------------------
ninp=numel(rise_inputs);

ainp=cell(1,ninp-1);

for ii=2:ninp

    ainp{ii-1}=eval(rise_inputs{ii});

end

end


function [sep_cf, cf]=do_memo(dd,varargin)

args=cell2mat(strcat(rdico.input_list,','));

args=[args,'tol,'];

sep_cf0=str2func(['@(',args(1:end-1),')[',dd,']']);

cf0=str2func(['@(',args(1:end-1),')all([',dd,'])']);

sep_cf=@engine_sepcf;

cf=@engine_cf;

    function out=engine_sepcf(y)

        out=sep_cf0(y,varargin{:});

    end


    function out=engine_cf(y)

        out=cf0(y,varargin{:});

    end

end
