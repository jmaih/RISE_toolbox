function info=describe(self)
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


nobs=self.NumberOfObservations;

% sort the data
%--------------
data=sort(self.data,1);

% header
%-------
info=set_row('',repmat(self.varnames,[1,1,size(data,3)]));

% mean
%-----
info=[info
    set_row('mean',utils.stat.nanmean(data,1))
    ];

% std
%-----
info=[info
    set_row('std',nanstd(data,0,1))
    ];

% min
%-----
info=[info
    set_row('min',utils.stat.nanmin(data,[],1))
    ];

% 25%
%-----
target=round(0.25*nobs);
info=[info
    set_row('25%',data(target,:,:))
    ];

% 50%
%-----
target=round(0.5*nobs);
info=[info
    set_row('50%',data(target,:,:))
    ];

% 75%
%-----
target=round(0.75*nobs);
info=[info
    set_row('75%',data(target,:,:))
    ];

% max%
%-----
info=[info
    set_row('max',utils.stat.nanmax(data,[],1))
    ];


    function c=set_row(name,data)
        cell_style=iscell(data);
        if cell_style
            c=[{name},data(:,:,1)];
        else
            c=[{name},num2cell(data(:,:,1))];
        end
        for ipage=2:size(data,3)
            if cell_style
                c=cat(3,c,[{name},data(:,:,ipage)]);
            else
                c=cat(3,c,[{name},num2cell(data(:,:,ipage))]);
            end
        end
    end
end