function main_frame=save_objects(filename,obj)

main_classes={'model','sensei','rise',...
    'dyn_ts','rise_time_series','rise_date','tseries',...
    'rise_param','rise_estim_param',...
    'rise_equation',...
    'rise_variable'};

main_frame=obj2struct(obj);
main_frame.class=class(obj);
save(filename,'main_frame')

    function main_frame=obj2struct(obj)
        nobj=numel(obj);
        main_frame=struct();
        fields=fieldnames(obj);
        for ii=1:numel(fields)
            for jj=1:nobj
                if isstruct(obj(jj).(fields{ii}))
                    main_frame(jj).(fields{ii})=struct2struct(obj(jj).(fields{ii}));
                elseif ismember(class(obj(jj).(fields{ii})),main_classes)
                    main_frame(jj).(fields{ii})=obj2struct(obj(jj).(fields{ii}));
                    [main_frame(jj).(fields{ii}).class]=deal(class(obj(jj).(fields{ii})));
                else
                    main_frame(jj).(fields{ii})=obj(jj).(fields{ii});
                end
            end
        end
    end

    function main_frame=struct2struct(strct)
        nstrct=numel(strct);
        main_frame=struct();
        fields=fieldnames(strct);
        for ii=1:numel(fields)
            for jj=1:nstrct
                if isstruct(strct(jj).(fields{ii}))
                    main_frame(jj).(fields{ii})=struct2struct(strct(jj).(fields{ii}));
                elseif ismember(class(strct(jj).(fields{ii})),main_classes)
                    main_frame(jj).(fields{ii})=obj2struct(strct(jj).(fields{ii}));
                    [main_frame(jj).(fields{ii}).class]=deal(class(strct(jj).(fields{ii})));
                else
                    main_frame(jj).(fields{ii})=strct(jj).(fields{ii});
                end
            end
        end
    end
end
