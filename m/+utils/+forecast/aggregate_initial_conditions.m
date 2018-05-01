function [y0,T,steady_state]=aggregate_initial_conditions(PAI,y0,T,steady_state)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<4
    steady_state=[];
    if nargin<3
        T=[];
    end
end
h=numel(PAI);
y00_=0;
s_state=0;
solve_order=1;
do_T= ~isempty(T);
do_ss= ~isempty(steady_state);
if do_T
    solve_order=size(T,1); 
    % the last row might be the position of the state variables but this
    % should not be a problem
end
for io=1:solve_order
    tsol=0;
    for ireg=1:h
        if io==1
            y00_=y00_+y0(ireg).y*PAI(ireg);
            if do_ss
                s_state=s_state+steady_state{ireg}*PAI(ireg);
            end
        end
        if do_T
            tsol=tsol+T{io,ireg}*PAI(ireg);
        end
    end
    if do_T
        T(io,:)=repmat({tsol},1,h);
    end
end
if do_ss
    steady_state=repmat({s_state},1,h);
end
y0(1).y=y00_;
y0(1:end)=y0(1);

end
