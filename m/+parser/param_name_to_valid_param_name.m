function [pname_out,capture_errors]=param_name_to_valid_param_name(pname)

% change the parameter names from name(chain,state) to name_chain_state

pname_out=parser.valid_names_in_text(pname);

% check that we now have valid variable names
%-------------------------------------------
char_flag=ischar(pname);
capture_flag=nargout>1;
pname_test=pname_out;
if ischar(pname_test)
    pname_test=cellstr(pname_test);
end
npar=numel(pname_test);
if capture_flag
    capture_errors=cell(1,npar);
    discard=false(1,npar);
end
for itest=1:numel(pname_test)
    if isvarname(pname_test{itest})
        if capture_flag
            discard(itest)=true;
        end
    else
        if char_flag
            offender=pname;
        else
            offender=pname{itest};
        end
        if capture_flag
            capture_errors{itest}=offender;
        else
            error(['"',offender,'" is not a valid parameter name'])
        end
    end
end
if capture_flag
    capture_errors(discard)=[];
end
% pname=strrep(pname,'(','_');
% pname=strrep(pname,',','_');
% pname=strrep(pname,')','');

