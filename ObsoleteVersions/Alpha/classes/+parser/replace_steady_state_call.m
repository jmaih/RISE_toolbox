function string=replace_steady_state_call(string,flag)
if nargin<2
    flag='';
end
loc_=strfind(string,'steady_state');
span=length('steady_state');
while ~isempty(loc_)
    loc_=loc_(1);
    left=string(1:loc_-1);
    right=string(loc_+span:end);
    closing=strfind(right,')');
    closing=closing(1);
    number=right(2:closing-1);
    right=right(closing+1:end);
    switch flag
        case 'symbolic'
            in_between=['ss_',number];
        otherwise
            in_between=['ss(',number,')'];
    end
    string=[left,in_between,right];
    loc_=strfind(string,'steady_state');
end
end
