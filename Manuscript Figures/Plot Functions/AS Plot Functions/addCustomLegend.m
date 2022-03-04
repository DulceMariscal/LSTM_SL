function addCustomLegend(strings, fontSize, x, y, dx, dy,...
         lineCol, lineWidth,...
         markSize, markType, markFaceCol, markEdgeCol, markLW)
% xlegs = [120, 300];
% ylegs = [-.07   .07];
% ylegs = [.1   2*.1];
for j = 1:length(y)
%     if j==1
%         cmark = markers.a1;
%     else
%         cmark = markers.a2;
%     end
    plot(x, y(j)*[1 1], 'color', lineCol, 'LineWidth', lineWidth)
    
    scatter(x(1), y(j), markSize, markFaceCol, 'filled', markType{j},...
    'MarkerEdgeColor', markEdgeCol, 'LineWidth', markLW);

    text(x(2) + dx, y(j) + dy, strings{j}, 'FontSize', fontSize,...
        'HorizontalAlignment','left','VerticalAlignment','middle')
end
end
