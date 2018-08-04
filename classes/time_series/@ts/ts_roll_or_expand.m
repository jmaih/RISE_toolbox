function db=ts_roll_or_expand(db,func,window,varargin)
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


data=db.data;

% create a function handle
%-------------------------
fhandle=ts_create_function_handle();

rolling=~isempty(window);
expanding=~rolling;

tmax=ts_set_sample_length();

for t=1:tmax
    d=ts_extract_data();
    cc=fhandle(d);
    if t==1
        sizcc=size(cc);
        ts_style=sizcc(1)==1;
        if ts_style
            m=nan([tmax,sizcc(2:end)]);
        else
            m=cell(1,tmax);
        end
    end
    if ts_style
        m(t,:,:)=cc;
    else
        m{t}=cc;
    end
end

new_start=1;
if rolling
    new_start=window;
end

db=ts(db.date_numbers(new_start:end),m,db.varnames);

    function tmax=ts_set_sample_length()
        tmax=size(data,1);
        if rolling
            tmax=tmax-window+1;
        end
    end

    function d=ts_extract_data()
        if expanding
            d=data(1:t,:,:);
        elseif rolling
            d=data(t+(0:window-1),:,:);
        end
    end

    function fh=ts_create_function_handle()
        char_type=ischar(func);
        if char_type
            func(isspace(func))=[];
            w=what('utils/stat');
            w=strrep(w.m,'.m','');
            if ~any(strcmp(func,w))
                disp(w(:)')
                error([func,' is not among the functions listed above'])
            end
            fh=@(x)utils.stat.(func)(x,varargin{:});
            % for the time series, the dates start at observation # window
        elseif isa(func,'function_handle')
            fh=@(x)func(x,varargin{:});
        else
            error('func must be a string or a function handle')
        end
    end
end