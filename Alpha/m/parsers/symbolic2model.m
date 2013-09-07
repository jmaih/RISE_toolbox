function [structural,shadow]=symbolic2model(symbolic,Lead_lag_incidence,word_list)
delimiters=[char([9:13,32]),'.(){};/*-+=^,[] '];
structural=symbolic;
shadow=symbolic;
for eq=1:size(symbolic,1)
    eqtn=symbolic{eq};
    shad='';
    strct='';
    while ~isempty(eqtn)
        [tok,rest]=strtok(eqtn,delimiters);
        if isempty(tok)
            shad=[shad,eqtn]; %#ok<*AGROW>
            strct=[strct,eqtn];
        else
            position = strfind(eqtn,tok); position=position(1);
            shad=[shad,eqtn(1:position-1)];
            strct=[strct,eqtn(1:position-1)];
            if strncmp(tok,'param_',6)
                shad=[shad,'param(',tok(7:end),')'];
                strct=[strct,word_list.parameters{str2double(tok(7:end))}];
            elseif strncmp(tok,'def_',4)
                shad=[shad,'def(',tok(5:end),')'];
                strct=[strct,word_list.definition_list{str2double(tok(5:end))}];
            elseif strncmp(tok,'y_',2)
                shad=[shad,'y(',tok(3:end),')'];
                nn=str2double(tok(3:end));
                [ii,jj]=find(Lead_lag_incidence==nn);
                strct=[strct,word_list.varendo{ii}];
                if jj-2~=0
                    strct=[strct,'{',int2str(jj-2),'}'];
                end
            elseif strncmp(tok,'x_',2)
                shad=[shad,'x(',tok(3:end),')'];
                strct=[strct,word_list.varexo{str2double(tok(3:end))}];
            elseif strncmp(tok,'ss_',3)
                shad=[shad,'ss(',tok(4:end),')'];
                nn=str2double(tok(4:end));
                strct=[strct,word_list.varendo{nn}];
            else
                shad=[shad,tok];
                strct=[strct,tok];
            end
        end
        eqtn=rest;
    end
    shadow{eq}=shad;
    structural{eq}=strct;
end