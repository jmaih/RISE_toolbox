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
    def=obj.solution.definitions{1};
    ss=obj.solution.ss{1};
    param=obj.parameter_values(:,1);
end
nrout=numel(routines);
hc=nrout>0 && obj.options.simul_honor_constraints;
if hc && ~any(obj.exogenous.shock_horizon>0)
    error('restrictions detected but no anticipatory behavior to satisfy them')
end
clear obj
sep_cf=separate_complementarity;
if nargout>1
    cf=@build_complementarity;
end

    function c=separate_complementarity()
        c=[];
        if ~hc
            return
        end
        dd=cell(1,nrout);
        for ii=1:nrout
            if ii==1
                x=[];
                sparam=[];
                s0=1;
                s1=1;
            end
            % those functions are expected to return scalar values
            %-----------------------------------------------------
            rout_i=func2str(routines{ii});
            right_par=find(rout_i==')',1,'first');
            if ii==1
                header=rout_i(1:right_par);
            end
            ge_=strfind(rout_i,'>=');
            if isempty(ge_)||numel(ge_)>1
                error('">=" missing or too many ">="')
            end
            dd{ii}=rout_i(right_par+1:ge_-1);
            right=rout_i(ge_+2:end);
            if ~strcmp(right,'0')
                error('right-hand side of constraints must be 0')
            end
            if ii==nrout
                dd=cell2mat(strcat(dd,';'));
                dd=str2func([header,'[',dd(1:end-1),']']);
                c=utils.code.anonymize_to_one_input(dd,x,ss,param,sparam,def,s0,s1);
            end
        end
    end

    function c=build_complementarity(y)
        c=true;
        if ~hc
            return
        end
        for ii=1:nrout
            if ii==1
                x=[];
                sparam=[];
                s0=1;
                s1=1;
            end
            c= c && routines{ii}(y,x,ss,param,sparam,def,s0,s1);
            if ~c
                break
            end
        end
    end

end