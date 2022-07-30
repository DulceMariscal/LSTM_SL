function [paradigm, epOrder] = loadParadigm(paradName, varargin)
if nargin>1
    d = varargin{1};
else
    d = 1;
end
switch paradName
    case 'ANA'
        n            = struct('epochs',4,'trialsPerEpoch',[150,600,600,600]);
        inc          = cell(1,n.epochs);
        epochs       = struct('u',inc,'tstart',inc,'id',{'N','A1','NL','A2'},...
            'printID', {'N','A_1','N_L','A_2'});
        transitions  = nan(1,n.epochs-1);
        midEpochs    = nan(1,n.epochs);
        epochs(1).u  = zeros(1,n.trialsPerEpoch(1));
        epochs(2).u  = d*ones(1,n.trialsPerEpoch(2));
        epochs(3).u  = zeros(1,n.trialsPerEpoch(3));
        epochs(4).u  = d*ones(1,n.trialsPerEpoch(4));
    case 'SHAB'
        n            = struct('epochs',4,'trialsPerEpoch',[150,120,600,600]);
        inc          = cell(1,n.epochs);
        epochs       = struct('u',inc,'tstart',inc,'id',{'N','A1','NL','A2'},...
            'printID', {'N','A_1','N_L','A_2'});
        transitions  = nan(1,n.epochs-1);
        midEpochs    = nan(1,n.epochs);
        epochs(1).u  = zeros(1,n.trialsPerEpoch(1));
        epochs(2).u  = d*ones(1,n.trialsPerEpoch(2));
        epochs(3).u  = zeros(1,n.trialsPerEpoch(3));
        epochs(4).u  = d*ones(1,n.trialsPerEpoch(4));
    case 'GR'
        n            = struct('epochs',4,'trialsPerEpoch',[150,600,600,600]);
        inc          = cell(1,n.epochs);
        epochs       = struct('u',inc,'tstart',inc,'id',{'N','A1','NL','A2'},...
            'printID', {'N','A_1','N_L','A_2'});
        transitions  = nan(1,n.epochs-1);
        midEpochs    = nan(1,n.epochs);
        epochs(1).u  = zeros(1,n.trialsPerEpoch(1));
        epochs(2).u  = [zeros(1,tclamp),linspace(0,d,n.trialsPerEpoch(2)-2*tclamp),zeros(1,tclamp)];
        epochs(3).u  = zeros(1,n.trialsPerEpoch(3));
        epochs(4).u  = d*ones(1,n.trialsPerEpoch(4));
    case 'Switch'
        n            = struct('epochs',25,'trialsPerEpoch',[50,450,300,30,50,30,100,30,20,30,15,30,10,30,10,30,10,30,10,30,10,150,100,30,100]);
        inc          = cell(1,n.epochs);
        epochs       = struct('u',inc,'tstart',inc,'id',{'N','A1','NL','A2','NL2','B1','NL3','A3','NL4','A4','NL5','A5','NL6','A6','NL7','A7','NL8','A8','NL9','A9','NL10','A10','NL11','B2','NL12'},...
            'printID', {'N','A1','NL','A2','NL2','B1','NL3','A3','NL4','A4','NL5','A5','NL6','A6','NL7','A7','NL8','A8','NL9','A9','NL10','A10','NL11','B2','NL12'});
        transitions  = nan(1,n.epochs-1);
        midEpochs    = nan(1,n.epochs);
        epochs(1).u  = zeros(1,n.trialsPerEpoch(1));
        epochs(2).u  = d*ones(1,n.trialsPerEpoch(2));
        epochs(3).u  = zeros(1,n.trialsPerEpoch(3));
        epochs(4).u  = d*ones(1,n.trialsPerEpoch(4));
        epochs(5).u  = zeros(1,n.trialsPerEpoch(5));
        epochs(6).u  = -d*ones(1,n.trialsPerEpoch(6));
        epochs(7).u  = zeros(1,n.trialsPerEpoch(7));
        epochs(8).u  = d*ones(1,n.trialsPerEpoch(8));
        epochs(9).u  = zeros(1,n.trialsPerEpoch(9));
        epochs(10).u  = d*ones(1,n.trialsPerEpoch(10));
        epochs(11).u  = zeros(1,n.trialsPerEpoch(11));
        epochs(12).u  = d*ones(1,n.trialsPerEpoch(12));
        epochs(13).u  = zeros(1,n.trialsPerEpoch(13));
        epochs(14).u  = d*ones(1,n.trialsPerEpoch(14));
        epochs(15).u  = zeros(1,n.trialsPerEpoch(15));
        epochs(16).u  = d*ones(1,n.trialsPerEpoch(16));
        epochs(17).u  = zeros(1,n.trialsPerEpoch(17));
        epochs(18).u  = d*ones(1,n.trialsPerEpoch(18));
        epochs(19).u  = zeros(1,n.trialsPerEpoch(19));
        epochs(20).u  = d*ones(1,n.trialsPerEpoch(20));
        epochs(21).u  = zeros(1,n.trialsPerEpoch(21));
        epochs(22).u  = d*ones(1,n.trialsPerEpoch(22));
        epochs(23).u  = zeros(1,n.trialsPerEpoch(23));
        epochs(24).u  = -d*ones(1,n.trialsPerEpoch(24));
        epochs(25).u  = zeros(1,n.trialsPerEpoch(25));
    case 'Switch2'
        n            = struct('epochs',9,'trialsPerEpoch',[150,30,50,30,50,90,150,450,300]);
        inc          = cell(1,n.epochs);
        epochs       = struct('u',inc,'tstart',inc,'id',{'N','SA1','N2','SB1','N3','Ramp','N4','A1','W'},...
            'printID', {'N','SA1','N2','SB1','N3','Ramp','N4','A1','W'});
        transitions  = nan(1,n.epochs-1);
        midEpochs    = nan(1,n.epochs);
        epochs(1).u  = zeros(1,n.trialsPerEpoch(1));
        epochs(2).u  = d*ones(1,n.trialsPerEpoch(2));
        epochs(3).u  = zeros(1,n.trialsPerEpoch(3));
        epochs(4).u  = -d*ones(1,n.trialsPerEpoch(4));
        epochs(5).u  = zeros(1,n.trialsPerEpoch(5));
        epochs(6).u  = [linspace(0,d,30) d*ones(1,60)];
        epochs(7).u  = zeros(1,n.trialsPerEpoch(7));
        epochs(8).u  = d*ones(1,n.trialsPerEpoch(8));
        epochs(9).u  = zeros(1,n.trialsPerEpoch(9));
        
    case 'Switch3'
        n            = struct('epochs',9,'trialsPerEpoch',[150,30,50,30,50,150,450,300,90]);
        inc          = cell(1,n.epochs);
        epochs       = struct('u',inc,'tstart',inc,'id',{'N','SA1','N2','SB1','N3','N4','A1','W','Ramp'},...
            'printID', {'N','SA1','N2','SB1','N3','N4','A1','W','Ramp'});
        transitions  = nan(1,n.epochs-1);
        midEpochs    = nan(1,n.epochs);
        epochs(1).u  = zeros(1,n.trialsPerEpoch(1));
        epochs(2).u  = d*ones(1,n.trialsPerEpoch(2));
        epochs(3).u  = zeros(1,n.trialsPerEpoch(3));
        epochs(4).u  = -d*ones(1,n.trialsPerEpoch(4));
        epochs(5).u  = zeros(1,n.trialsPerEpoch(5));
       
        epochs(6).u  = zeros(1,n.trialsPerEpoch(6));
        epochs(7).u  = d*ones(1,n.trialsPerEpoch(7));
        epochs(8).u  = zeros(1,n.trialsPerEpoch(8));
         epochs(9).u  = [linspace(0,d,30) d*ones(1,60)];
        
    otherwise
        warning('unrecognized paradigm')
end
n.ttot = sum(n.trialsPerEpoch);
t      = 1:n.ttot;
% Store epoch start
e = 1;
epochs(e).tstart = 1;
epochs(e).tend   = epochs(e).tstart + n.trialsPerEpoch(e) - 1;
midEpochs(e)     = epochs(e).tstart + round(.5*n.trialsPerEpoch(e) - 1);
for e = 2:n.epochs
    epochs(e).tstart = epochs(e-1).tstart + n.trialsPerEpoch(e-1);
    epochs(e).tend   = epochs(e).tstart   + n.trialsPerEpoch(e) - 1;
    transitions(e-1) = epochs(e).tstart;
    midEpochs(e)     = epochs(e).tstart + round(.5*n.trialsPerEpoch(e) - 1);
end

% Concatenate u
u = [epochs(:).u];
epOrder  = cell2struct(num2cell(1:n.epochs)',{epochs.id}');
paradigm = struct('epochs',epochs,'transitions',transitions,...
    'midEpochs',midEpochs,'n',n,...
    'u',u,'t',t);
end