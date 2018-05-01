function is_changed=parameters_changed_in_ssmodel(shadow_sstate,pstring,n)
% parameters_changed_in_ssmodel - flags the parameters that are modified in
% the steady state model
%
% ::
%
%
%   is_changed=parameters_changed_in_ssmodel(shadow_sstate,pstring,n)
%
% Args:
%
%    - **shadow_sstate** [cellstr]: list of steady state equations
%
%    - **pstring** [char]: generic name of shadow parameters
%
%    - **n** [integer]: number of parameters
%
% Returns:
%    :
%
%    - **is_changed** [logical]: 1 x n vector flagging the changed
%      parameters
%
% Note:
%
% Example:
%
%    See also:

len=size(pstring,2);
bingo=find(strncmp(pstring,shadow_sstate,len));
is_changed=false(1,n);
for irow=1:numel(bingo)
    thisrow=shadow_sstate{bingo(irow)};
    eqs=find(thisrow=='=');
    pos=str2double(thisrow(len+2:eqs-2));
    is_changed(pos)=true;
end

end
