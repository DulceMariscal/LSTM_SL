function addTexts(sph, xs, ys, names, fs)
    ns = length(xs);
    for i=1:ns
        text(sph, xs(i),ys(i),names{i},'FontSize',fs,'HorizontalAlignment','center');
    end
end
