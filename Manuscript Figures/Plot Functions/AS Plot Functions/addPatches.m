function ph = addPatches(axh, col, xleft, xright, ybot, ytop, dx, dy, alpha) 
nPatches   = size(col,1);
for ip = 1:nPatches        
        x([1 4]) = xleft(ip)  - dx; % X Start patch
        x([2 3]) = xright(ip) + dx; % X End patch
        y([1 2]) = ybot(ip)   - dy;         % Y Start patch
        y([3 4]) = ytop(ip)   + dy;         % Y End patch
        hold on, ph(ip) = patch(axh,x,y,col(ip,:),'EdgeColor','none',...
                                'FaceAlpha', alpha);
        uistack(ph(ip),'bottom')
end
end