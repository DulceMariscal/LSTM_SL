function [defCols, cbCols, grayCols, defColsMat, cbColsMat, grayColsMat] = getColors(varargin)
defColsMat = [get(groot,'defaultAxesColorOrder'); 0 0 0];
cbColsMat  = cbrewer('qual', 'Set1', 9);
dx=.1;
grayColsMat   = repmat([0:dx:.9]',1, 3);

defColNames = {'blue','orange','yellow','purple','green','liBlue','red','black'};
cbColNames  = {'red','blue','green','purple','orange','yellow',...
               'brown', 'pink', 'gray'};
ngr = size(grayColsMat,1);
grayColsNames = cell(1,ngr);
for ig=1:ngr
    grayColsNames{ig} = ['gr' num2str(10*(ig-1)*dx)];
end
defCols  = cell2struct(num2cell(defColsMat,2),   defColNames);
cbCols   = cell2struct(num2cell(cbColsMat,2),    cbColNames);
grayCols = cell2struct(num2cell(grayColsMat,2),  grayColsNames);

if nargin>0
    PLOT = varargin{1};
    if PLOT
        height=1;
        figure,
        
        %SetA------------------------------------------------------------
        subplot(3,1,1)
        n = size(defColsMat,1);
        for i=1:n
            hold on
            scatter(i,height,100,defColsMat(i,:),'filled')
        end
%         legend(defColNames,'NumColumns',n,'Location','north','Box','off');
        xticks([1:n])
        xticklabels(defColNames)
        xtickangle(45)
        title('Matlab''s default colors');
        
        %SetB------------------------------------------------------------
        subplot(3,1,2)
        n = size(cbColsMat,1);
        for i=1:n
            hold on
            scatter(i,height,100,cbColsMat(i,:),'filled')
        end
%         legend(cbColNames,'NumColumns',n,'Location','north','Box','off');
        xticks([1:n])
        xticklabels(cbColNames)
        xtickangle(45)
        title('CBrewer''s qual-Set1 colors');
        
        %SetC------------------------------------------------------------
        subplot(3,1,3)
        n = ngr;
        for i = 1:n
            hold on
            scatter(i,height,100,grayColsMat(i,:),'filled')
        end
%         legend(cbColNames,'NumColumns',n,'Location','north','Box','off');
        xticks([1:n])
        xticklabels(grayColsNames)
        xtickangle(45)
        title('Gray shades');

    end
end
end