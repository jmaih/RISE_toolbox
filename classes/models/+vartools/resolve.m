function [B,SIG,SIG05]=resolve(solution,nlags,h)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 


B=cell(1,h);
SIG=cell(1,h);
SIG05=cell(1,h);
is_svar=isfield(solution,'a0');
has_det=isfield(solution,'c');
for ireg=1:h
    for ilag=1:nlags
        B{ireg}=[B{ireg},solution.(sprintf('a%0.0f',ilag)){ireg}];
    end
    if has_det
        B{ireg}=[B{ireg},solution.c{ireg}];
    end
    SIG05{ireg}=solution.omg{ireg}*solution.sig{ireg};
    if is_svar
        B{ireg}=solution.a0{ireg}\B{ireg};
        SIG05{ireg}=solution.a0{ireg}\SIG05{ireg};
    end
    SIG{ireg}=SIG05{ireg}*SIG05{ireg}';
end