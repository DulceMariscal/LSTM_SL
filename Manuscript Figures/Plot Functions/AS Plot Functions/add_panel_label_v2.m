function finalPositions = add_panel_label_v2(sph, varargin)
spnames = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W', 'X', 'Y', 'Z'};


if isa(sph,'double')
    [nsp, d4] = size(sph);
    assert(d4==4);
    dt=1;
elseif strcmp(get(sph(1),'type'),'axes')
    nsp = length(sph);
    dt=2;
elseif strcmp(get(sph(1),'type'),'polaraxes')
    nsp = length(sph);
    dt=3;
else
    error('Unrecognized datatype')
end

defOptions = struct('dx', .04, 'dy', .04, 'labels', {spnames(1:nsp)}, 'yl', 1,...
                    'FontSize', 22, 'FontWeight', 'bold', 'FlippedYaxis', false(1,nsp),'tiles',false,...
                    'xAlignment',[],'yAlignment',[],'latex',false);
                                
if ~isempty(varargin)
    inOptions = varargin{1};
    options = myParseArgs(inOptions,  defOptions);
else
    options = defOptions;
end

if strcmp(options.FontWeight,'bold')
    if options.latex
        options.labels = cellfun(@(x) ['\bf{' x '}'],options.labels,'UniformOutput',false);
    end
end
if options.tiles
    dt = 4;
end

if length(options.dx)>1
    cdx = options.dx;
else
    cdx = options.dx*ones(1,nsp);
end
if length(options.dy)>1
    cdy = options.dy;
else
    cdy = options.dy*ones(1,nsp);
end

if dt==4
  
    oldPos = nan(nsp,4);
    
    % Get old positions    
    for sp = 1:nsp
       oldPos(sp,:) = sph(sp).Position;
    end
    corPos = oldPos;
    
    % Correct x
    if ~isempty(options.xAlignment)
        nxalign = length(options.xAlignment);
        for ixa = 1:nxalign
            corPos(options.xAlignment{ixa},1) = min(oldPos(options.xAlignment{ixa},1));
        end
    end
    
    % Correct y
    if ~isempty(options.yAlignment)
        nyalign = length(options.yAlignment);
        for iya = 1:nyalign
            corPos(options.yAlignment{iya},2) = min(oldPos(options.yAlignment{iya},2));
        end
    end
    
end

for sp=1:nsp
    if options.yl==0 % Use axes position as reference
        switch dt
            case 1
                cpos = sph(sp,:);
            case 2
                cpos = plotboxpos(sph(sp));
            case 3
                cpos = sph(sp).Position;
            case 4
                cpos = corPos(sp,:);
        end
        %     cpos = get(sph(sp), 'Position');

        cpos(1) =  cpos(1) - cdx(sp);
        cpos(2) =  cpos(2) + cdy(sp);
        annotation('textbox',cpos,'String',options.labels{sp},'FitBoxToText','on',...
            'linestyle', 'none','FontSize', options.FontSize, 'FontWeight',options.FontWeight);
        finalPositions(sp).x = cpos(1);
        finalPositions(sp).y = cpos(2);
    else %Default
        %         yl.Pos2    = ax.YLabel.Position
        %         titl.Pos2 = ax.Title.Position
        ax = sph(sp);
        yl.Ex    = ax.YLabel.Extent;
        titl.Ex  = ax.Title.Extent;

        xunit = diff(ax.XLim)/100;
        yunit = diff(ax.YLim)/100;

        xt = yl.Ex(1)    + cdx(sp)*xunit;
        if options.FlippedYaxis(sp)
            yt = (titl.Ex(2) - titl.Ex(4)) - cdy(sp)*yunit;
        else
            yt = (titl.Ex(2) + titl.Ex(4)) + cdy(sp)*yunit;
        end
        th = text(ax, xt, yt, options.labels{sp}, 'FontSize', options.FontSize, ...
            'FontWeight',options.FontWeight,   'VerticalAlignment', 'bottom');
        [xtn, ytn] = data2norm(ax, xt, yt);
        
        finalPositions(sp).x = xt;
        finalPositions(sp).y = yt;
        finalPositions(sp).xnorm = xtn;
        finalPositions(sp).ynorm = ytn;
    end
    
    
end
end