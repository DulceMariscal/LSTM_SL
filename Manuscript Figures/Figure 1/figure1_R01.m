clc
% close all
clear all

%% Load system parameters
p = mfilename('fullpath'); %Full path of current file

load('./Manuscript Figures/Data/Parameters/s82 learningTrace/allv.mat',...
    'systpEst','estimPar','usts');
indValStates = [1 2 5 6 9 10]; %Valid states
sFuns        = struct('syspec',@syspec_ssim_v15,'ffun',@ffun_ssim_v15,'gfun',@gfun_ssim_v13);
[cols.default, cols.cbr, cols.grays] = getColors();

systpEst.p(22)=0.8;

% cols.default = get(groot,'defaultAxesColorOrder');
n.states = length(indValStates);

%% Other initializations
fs.def   = 20;
fs.sma   = 18;
lw.mark  = 2;
lw.lines = 3;
my_set_default(fs.def,lw.lines,10);

% Load ANA paradigm------------------------------------------------------
d              = abs(usts.extra.dest.IS);
% [parad, names] = loadParadigm('ANA',d);
[parad, names] = loadParadigm('Switch',d);

simp = struct('oneStep', false, 'nt', parad.n.ttot,'u', parad.u,...
              't', parad.t,  'udt2', [], 'Ts', []);
          
% Simulate system
cout      = mySystSim(sFuns, systpEst, simp);
parad.y   = cout.y;
parad.z   = parad.u - parad.y;
statesAmp = diag([1 1 1 1 100 100]);
parad.x   = statesAmp*cout.x(indValStates,:);
legs.C0   = systpEst.xShNames(indValStates);
legs.C    = {'x_f','x_s','\delta_+','\delta_-','l_+','l_-'};

%% PLOT
% Initializations--------------------------------------------------------
close all
cols.u      = cols.default.black;
cols.y      = cols.default.yellow;
cols.z      = cols.grays.gr5;
cols.zMark  = brighten(cols.z,-.5);

cols.states = [cols.default.liBlue; cols.default.purple; cols.cbr.green; cols.cbr.red; cols.cbr.green; cols.cbr.red];
liSt.states = {'-','-','-','-',':',':'}';
liWi.states = num2cell([3 3 3 3 3 3])';
markers.a1  = '^';
markers.a2  = 'v';
labels.y    = {'perturbation magnitude', '(a.u.)'};
labels.x    = {'trial (t)'}; %,'trial after pert. (\tau)'};
legs.A      = {'u(t): perturbation','z(t): motor output','e(t): kin. error'};
legs.B      = {'z(t_{A_1})','z(t_{A_2})'};
% legs.D      = {'\delta(t_{A_1})','z(t_{A_2})'};
legs.D1     = {'\delta_+(t_{A_1})','\delta_+(t_{A_2})'};
legs.D2     = {'\delta_-(t_{A_1})','\delta_-(t_{A_2})'};
sz.dot      = 50; 
figure('Color','white','Units','normalized','Position',[0 0 .9 .9]);

%% S11 (A) perturbation + output + error------------------------------------
sp(1,1)     = subplot(3,1,1);
plot(nanTrans(parad.u,parad.transitions),'color', cols.u), hold on, 
plot(nanTrans(parad.z,parad.transitions),'color', cols.z),
plot(nanTrans(parad.y,parad.transitions),'color', cols.y)
% Dots
tA1     = parad.epochs(names.A1).tstart:parad.epochs(names.A1).tend;
% tA2     = parad.epochs(names.A2).tstart:parad.epochs(names.A2).tend;
tstA1   = tA1(2); 
scatter(tstA1, parad.z(tstA1),sz.dot,cols.z,'filled',...
    markers.a1,'MarkerEdgeColor',cols.zMark,'LineWidth',lw.mark);
% tstA2   = tA2(2); 
% scatter(tstA2, parad.z(tstA2),sz.dot,cols.z,'filled',...
%     markers.a2,'MarkerEdgeColor',cols.zMark,'LineWidth',lw.mark);
% % Extra
leg1h = legend(legs.A,'Location','southeast','Box','off',...
    'AutoUpdate','off','FontSize',fs.sma);
xticklabels('')
ylabel(labels.y)
box off;
axis tight



%% S21 (C) details about the states----------------------------------------
sp(2,1) = subplot(3,1,2);
plh = plot(parad.t,parad.x);
set(plh, {'Color'}, mat2cell(cols.states, ones(1,n.states), 3 ));
set(plh, {'LineWidth'}, liWi.states);
set(plh, {'LineStyle'}, liSt.states);
leg1h = legend(legs.C,'Location','southeast','Box','off',...
    'AutoUpdate','off','FontSize',fs.sma,'NumColumns',3);
xlabel(labels.x)
ylabel(labels.y)
box off
axis tight

%%

sp(3,1) = subplot(3,1,3);
plot(nanTrans(-parad.y,parad.transitions))
hold on 
% plot(nanTrans(parad.u,parad.transitions),'color', cols.u)
% set(plh, {'Color'}, mat2cell(cols.states, ones(1,n.states), 3 ));
% set(plh, {'LineWidth'}, liWi.states);
% set(plh, {'LineStyle'}, liSt.states);
% leg1h = legend(legs.C,'Location','southeast','Box','off',...
%     'AutoUpdate','off','FontSize',fs.sma,'NumColumns',3);
xlabel(labels.x)
ylabel('SLA')
box off
axis tight


%% S12 (B) details about A1 and A2-----------------------------------------
% sp(1,2) = subplot(2,6,6); 
% deltA1  = parad.n.trialsPerEpoch(names.A1);
% tA      = 2:deltA1; 
% plot(tA, parad.z(tA1(2:end)),'color', cols.z), hold on,
% plot(tA, parad.z(tA2(2:end)),'color', cols.z)
% % Add dots
% scatter(2, parad.z(tA1(2)),sz.dot,cols.z,'filled',markers.a1,...
%     'MarkerEdgeColor',cols.zMark,'LineWidth',lw.mark)
% scatter(2, parad.z(tA2(2)),sz.dot,cols.z,'filled',markers.a2,...
%     'MarkerEdgeColor',cols.zMark,'LineWidth',lw.mark)
% % lh = legend(legs.B,'Location','southeast','Box','off','AutoUpdate','off');
% box off;
% % 
% 
% % Add custom legend
% x0        = 350+30;
% addCustomLegend(legs.B, fs.sma, x0+[0 70], .075*[1 2.2], 5, -.02, ...
%             cols.z, lw.lines,...
%             sz.dot, {markers.a1, markers.a2}, cols.z, cols.zMark, lw.mark)
% 
% %% S22 (D) details about what states retain the memory-----------------------
% sp(2,2) = subplot(2,6,12);
% indxD   = [3 4];
% nxD     = length(indxD);
% plh1    = plot(tA,parad.x(indxD,tA1(2:end))); hold on
% plh2    = plot(tA,parad.x(indxD,tA2(2:end)));
% for ch = [plh1, plh2] %1:length(plh)
% %     ch = plh{ih};
%     set(ch, {'Color'},     mat2cell(cols.states(indxD,:), ones(1,nxD), 3));
%     set(ch, {'LineWidth'}, liWi.states(indxD));
%     set(ch, {'LineStyle'}, liSt.states(indxD));
% end
% xlabel(labels.x)
% box off;
% dys = [0 0]; %[+.0001, -.0001];
% for st = 1:2
%     colCState = cols.states(indxD(st),:);
%     darkColCS = brighten(colCState,-.5);
% 
% % A1
% sc1 = scatter(2, parad.x(indxD(st),tA1(2)) + dys(st), sz.dot, colCState, 'filled', markers.a1,...
%     'MarkerEdgeColor', darkColCS,'LineWidth', lw.mark);
% % A2
% sc2 = scatter(2, parad.x(indxD(st),tA2(2)), sz.dot, colCState, 'filled', markers.a2,...
%     'MarkerEdgeColor', darkColCS,'LineWidth', lw.mark);
% if st==2
%     uistack(sc1,'bottom');
%     uistack(sc2,'bottom');
% end
% end
% 
% % Add custom legend
% % delta+
% ccol      = cols.states(indxD(1),:);
% darkccol  = brighten(ccol,-.5);
% y0        = .013;
% dx        = .025;
% dy        = .0225;
% dxtext    = 5;
% dytext    = 0;
% % addCustomLegend(legs.D1, fs.sma, 350+[0 70], y0 + [2*dy 3*dy], dxtext, dytext, ...
% %             ccol, lw.lines,...
% %             sz.dot, {markers.a1, markers.a2}, ccol, darkccol, lw.mark)
% addCustomLegend(legs.D1, fs.sma, x0-100-150 + [0 70], y0 + [0, dy], dxtext, dytext, ...
%     ccol, lw.lines,...
%     sz.dot, {markers.a1, markers.a2}, ccol, darkccol, lw.mark)
% % delta-
% ccol      = cols.states(indxD(2),:);
% darkccol  = brighten(ccol,-.5);
% addCustomLegend(legs.D2, fs.sma, x0+[0 70], y0 + [0, dy], dxtext, dytext, ...
%             ccol, lw.lines,...
%             sz.dot, {markers.a1, markers.a2}, ccol, darkccol, lw.mark)
% 
% %% Finish up figure ----------------------------------------------------
% linkaxes(sp(:,1),'y'),
% % linkaxes(sp(:,2),'y'),
% ylAC = 1*max(abs(parad.u));
% ylB  = [0 ylAC];
% dy   = .01;
% ylD  = [min(parad.x(indxD,tA2(2:end)),[],'all') - dy,...
%         max(parad.x(indxD,tA2(2:end)),[],'all') + dy];
% 
% sp(1,1).YLim = [-ylAC ylAC];
% sp(1,2).YLim = ylB; 
% sp(2,2).YLim = ylD; 
% 
% ym2          = mean(ylB);
% linkaxes(sp(:,1),'x');
% linkaxes(sp(:,2),'x');
% 
% % Add patches---------------------------------------------------------
% cols.Apatch = [.7 .7 .7];
% % To (A)
% addPatches(sp(1,1), repmat(cols.Apatch,2,1),...
%     [parad.epochs([2 4]).tstart], [parad.epochs([2 4]).tend],...
%      -ylAC*[1 1], ylAC*[1 1], .5, .1, .3)
%  
% % To (B)
% addPatches(sp(1,2), cols.Apatch,...
%      2, 599,...
%     -ylAC, ylAC, .5, .1, .3)
% 
% % To (C)
% addPatches(sp(2,1), repmat(cols.Apatch,2,1),...
%     [parad.epochs([2 4]).tstart], [parad.epochs([2 4]).tend],...
%      -ylAC*[1 1], ylAC*[1 1], .5, .1, .3)
% 
%  % To (D)
% addPatches(sp(2,2), cols.Apatch,...
%      2, 599,...
%     -ylAC, ylAC, .5, .1, .3)
% 
% % Add epoch names
% addTexts(sp(1,1), parad.midEpochs,.5*ylAC*ones(1,parad.n.epochs),{parad.epochs.printID},fs.sma)
% addTexts(sp(1,2), 300, mean([ym2,ylB(2)]),{'A'},fs.sma)
% 
% %% Adjust width and x coordinate of second column
% tlong   = parad.n.ttot;
% tshort  = 600;
% p1      = sp(1,1).Position;
% sp(1,2).Position = sp(1,2).Position.*[0 1 0 1] + [.7 0 p1(3)*tshort/tlong 0];
% sp(2,2).Position = sp(2,2).Position.*[0 1 0 1] + [.7 0 p1(3)*tshort/tlong 0];
% 
% %% Adjust height of second row
% sp(2,1).Position = sp(2,1).Position.*[1 0 1 1] + [0 .15 0 0];
% sp(2,2).Position = sp(2,2).Position.*[1 0 1 1] + [0 .15 0 0];
% 
% %% Add panel labels
% spvec = reshape(sp',4,1);
% finalPositions = add_panel_label_v2(spvec);
% 
% %% Export figure
% % export_fig Figure1_v1 -pdf
% 
% %% Functions---------------------------------------------------------------
% % Add names-----------------------------------------------------------
% 
% % S12
% % subplot(2,6,5)
% % figure, plot(parad.u), hold on, plot(parad.u-parad.y)
% 
