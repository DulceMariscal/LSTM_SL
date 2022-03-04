function [xn, yn] = data2norm(ax, x, y)
axPos = get(ax, 'Position');

xRange = get(ax,'XLim');
yRange = get(ax,'YLim');

dxR    = diff(xRange);
dyR    = diff(yRange);

xn     = (x./dxR).*axPos(3) + axPos(1);
yn     = (y./dyR).*axPos(4) + axPos(2);


end