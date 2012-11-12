function correction=rise_saveaspdf(fig,filename)%%
theDot=find(filename=='.');
if ~isempty(theDot)
    filename=filename(1:theDot-1);
end

orient(fig,'landscape')

% rise_pdflatex=getappdata(0, 'rise_pdflatex');

success=false;
% if rise_pdflatex ||
% this thing lets me down a times
for itry=1:10
    try
        print(fig,'-depsc',filename)
        success=true;
        break
    end
end
if success
    % mask the filename, with "" in case it contains spaces.
    system(['epstopdf "',filename,'.eps"']);
    correction=-90;
    delete([filename,'.eps'])
end
% end

if ~success %|| ~ rise_pdflatex
    print(fig,'-dpdf',filename)%,sprintf('-r%d',dpi)
    correction=0;
end


% % % % return
% % % % Items={'PaperUnits','inches'
% % % %     'Units','inches'
% % % % %     'PaperOrientation','landscape'
% % % %     'PaperPosition',[]
% % % %     };
% % % % % Collect previous settings
% % % % for pass=1:2
% % % %     for it=1:size(Items,1)
% % % %         if pass==1
% % % %             eval(['pre_',Items{it,1},'=get(fig,Items{it,1});'])
% % % %         elseif pass==2 && ~isempty(Items{it,2})
% % % %             set(fig,Items{it,1},Items{it,2})
% % % %         end
% % % %     end
% % % % end
% % % %
% % % % % rect = [left, bottom, width, height];
% % % % PaperSize=get(fig,'PaperSize');
% % % % left=0.1;
% % % % bottom=left;
% % % % width=PaperSize(1)-2*left;
% % % % height=PaperSize(2)-2*bottom;
% % % % NewPosition=[left, bottom, width, height];
% % % % % Set the page size and position to match the figure's dimensions
% % % % set(fig,'PaperPosition',NewPosition);
% % % %
% % % % % Save the pdf (this is the same method used by "saveas")
% % % % print(fig,'-dpdf',filename)%,sprintf('-r%d',dpi)
% % % %
% % % % % Restore the previous settings
% % % % for it=1:size(Items,1)
% % % %     eval(['set(fig,Items{it,1},pre_',Items{it,1},');'])
% % % % end
% % % %
% % % % %{
% % % % Items={'PaperType','<custom>'
% % % %     'PaperUnits','inches'
% % % %     'Units','inches'
% % % % % % %     'PaperOrientation','landscape'
% % % %     'PaperPosition',[]
% % % %     'PaperSize',[]
% % % %     };
% % % % % Collect previous settings
% % % % for pass=1:2
% % % %     for it=1:size(Items,1)
% % % %         if pass==1
% % % %             eval(['pre_',Items{it,1},'=get(fig,Items{it,1});'])
% % % %         elseif pass==2 && ~isempty(Items{it,2})
% % % %             set(fig,Items{it,1},Items{it,2})
% % % %         end
% % % %     end
% % % % end
% % % %
% % % % position = get(fig,'Position');
% % % % % Set the page size and position to match the figure's dimensions
% % % % set(fig,'PaperPosition',[0,0,position(3:4)]);
% % % % set(fig,'PaperSize',position(3:4));
% % % %
% % % % theDot=find(filename=='.');
% % % % if ~isempty(theDot)
% % % %     filename=filename(1:theDot-1);
% % % % end
% % % %
% % % % % Save the pdf (this is the same method used by "saveas")
% % % % print(fig,'-dpdf',filename)%,sprintf('-r%d',dpi)
% % % %
% % % % % Restore the previous settings
% % % % for it=1:size(Items,1)
% % % %     eval(['set(fig,Items{it,1},pre_',Items{it,1},');'])
% % % % end
% % % %
% % % % print(fig,'-depsc',filename)%,sprintf('-r%d',dpi)
% % % % system(['epstopdf ',filename])
% % % %
% % % %
% % % % %}
