function correction=saveaspdf(fig,filename)
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

filename=parser.remove_file_extension(filename);

orient(fig,'landscape')

correction=0;
success=false;
switch computer
    case 'MAC'
        remove_figure_margins()
%    case {'PCWIN','PCWIN64'}
%        % this thing lets me down at times
%        for itry=1:10
%            try %#ok<TRYNC>
%                print(fig,'-depsc',filename)
%                success=true;
%                break
%            end
%        end
%        if success
%            rise_epstopdf=getappdata(0,'rise_epstopdf');
%            % mask the filename, with "" in case it contains spaces.
%            system([rise_epstopdf,' "',filename,'.eps"']);
%            correction=-90;
%            delete([filename,'.eps'])
%        end
    otherwise
        remove_figure_margins()
end

if ~success 
    print(fig,'-dpdf',filename)%,sprintf('-r%d',dpi)
end


    function remove_figure_margins()
        
        haxes=	findobj(fig, 'Type', 'axes', '-and', 'Tag', '');
        for n=1:length(haxes)
            xl=		get(haxes(n), 'XLim');
            yl=		get(haxes(n), 'YLim');
            lines=	findobj(haxes(n), 'Type', 'line');
            for m=1:length(lines)
                x=				get(lines(m), 'XData');
                y=				get(lines(m), 'YData');
                inx=			(xl(1) <= x) & (x <= xl(2));	% Within the x borders.
                iny=			(yl(1) <= y) & (y <= yl(2));	% Within the y borders.
                keep=			inx & iny;						% Within the box.
                
                if(~strcmp(get(lines(m), 'LineStyle'), 'none'))
                    crossx=		((x(1:end-1) < xl(1)) & (xl(1) < x(2:end))) ...	% Crossing border x1.
                        |	((x(1:end-1) < xl(2)) & (xl(2) < x(2:end))) ...	% Crossing border x2.
                        |	((x(1:end-1) > xl(1)) & (xl(1) > x(2:end))) ...	% Crossing border x1.
                        |	((x(1:end-1) > xl(2)) & (xl(2) > x(2:end)));	% Crossing border x2.
                    crossy=		((y(1:end-1) < yl(1)) & (yl(1) < y(2:end))) ...	% Crossing border y1.
                        |	((y(1:end-1) < yl(2)) & (yl(2) < y(2:end))) ...	% Crossing border y2.
                        |	((y(1:end-1) > yl(1)) & (yl(1) > y(2:end))) ...	% Crossing border y1.
                        |	((y(1:end-1) > yl(2)) & (yl(2) > y(2:end)));	% Crossing border y2.
                    crossp=	[(	(crossx & iny(1:end-1) & iny(2:end)) ...	% Crossing a x border within y limits.
                        |	(crossy & inx(1:end-1) & inx(2:end)) ...	% Crossing a y border within x limits.
                        |	crossx & crossy ...							% Crossing a x and a y border (corner).
                        ),	false ...
                        ];
                    crossp(2:end)=	crossp(2:end) | crossp(1:end-1);		% Add line segment's second end point.
                    keep=			keep | crossp;
                end
                set(lines(m), 'XData', x(keep))
                set(lines(m), 'YData', y(keep))
            end
        end
        
    end

end