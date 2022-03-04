%% PUPIL WIDTH BEHAVIOURAL ANALYSIS
clc
close all
clear all

%Change directory
cd('C:\Users\Alessandro S\Dropbox\BitBrainShared\BitBrain')

%% Select subject and experiment
subjectsToConsider=[1 2 6 7]; %1 2 6 7 better signal
experiment=2;

nsubjs=length(subjectsToConsider);
expstr=num2str(experiment);
current_dir=pwd;

%Initialize variables
nch=1; %TODO nch=3, considering also x and y position of eye tracker

%Parameters
nanVal=-1;
baselineInt=[-500 0]*10^-3;
baselineActivityInt=[-1.38 , 0]; %The median filter added noise at the border!! Remove it by adding padding! CHECK IF I'M NOT NORMALIZING BEFORE COMPUTING THIS VALUE
alpha=0.05;
nbOfTrialsToSmooth=5;
ntrialsSameState=4;
analysis=1;  
%1--> MW vs OT
%2--> Increase of PW after response (1.2 sec)

%% SELECT TRIALS BEFORE TARGET ERROR AND TARGET CORRECT-----------------------------------
%Analysis' parameters
%Load experiment's parameters
path=[pwd '\Storage Data SART experiment\' num2str(experiment) '\Parameters\param1.mat'];
load(path);
load('eventsDefinition.mat');
load('logsDefinition.mat');




%Initialize variables
Mbtctot=cell(1,nch); %Will contain, for each channel, all the trials before target correct
Mbtetot=cell(1,nch);
MbtctotSUB=cell(nsubjs,nch);
MbtetotSUB=cell(nsubjs,nch);
LATENCIESot=cell(1,nsubjs);
LATENCIESmw=cell(1,nsubjs);
LATENCIESotTot=[];
LATENCIESmwTot=[];
ERsubsOT=cell(1,nch);
ERsubsMW=cell(1,nch);
baselineActivityValueOT=[];
baselineActivityValueMW=[];

% Analysis loop
for isub=1:nsubjs
    isub
    csub=subjectsToConsider(isub);
    substr=num2str(csub);
    
    %Load subject's trials-------------------------------------------------
    pathname=[current_dir '\Storage Data SART experiment\' expstr '\S' substr '\'];
    if analysis==1
        fname=strcat('widthData_s',substr,'.mat'); %Contains 'subjectPWallblocks' with kept clean width data
    elseif analysis==2
        fname=strcat('RESPwidthData_s',substr,'.mat'); %Contains 'subjectPWallblocks' with kept clean width data
    end
    dataFileName=[pathname fname];
    load(dataFileName);
    
    %Upload subject's raw data---------------------------------------------
    fname1=strcat('\rawData_s',substr,'.mat');
    rawDataFileName=[pathname fname1];
    load(rawDataFileName);           %E, L, ET, EEG, (BVP, GSR)    (1 cell for each block)
    clear('EEG','BVP','GSR');
    
    %Upload behavioural matrix--------------------------------------------
    fnameBmat=strcat('Bmat_s',substr,'.mat');
    BmatFileName=[pathname fnameBmat];
    load(BmatFileName); %Upload B matrix
    
    %Select trials preceding a target correct (BTC) and those preceding target errors (BTE)
    [totalBTC, totalBTE, latenciesBTC, latenciesBTE,caBTC, caBTE ] = ...
        findERPtargIndicesLoopVersion(B,blocksToAverage{csub},CorrTrial,ShownNumber,mystate,ntrialsSameState);
    
    %Select ET periods
    tobiiData=ET{1,block};
    M=cleanETdata(tobiiData,nanVal); %Average valid rigth-left values
    
    
    %Filter out selected data whose trials have been discarded
    [filtBTC, indbtc, indbtc2]=intersect(subjectPWallblocks{2},totalBTC); %Indices of original trials
    [filtBTE, indbte, indbte2]=intersect(subjectPWallblocks{2},totalBTE);
    
    %Store subject latencies
    LATENCIESot{1,isub}=latenciesBTC(indbtc2);
    LATENCIESmw{1,isub}=latenciesBTE(indbte2);
    if ~isempty(filtBTC) && ~isempty(filtBTE)
        
        LATENCIESotTot=[LATENCIESotTot; latenciesBTC(indbtc2)];
        LATENCIESmwTot=[LATENCIESmwTot; latenciesBTE(indbte2)];
    end
    
    % Load block delay information
    blockFilename=[pathname 'delays.mat'];
    load(blockFilename);
    
    %% Select channels' data
    for channel=1:nch
        
        %Compute ER response for each subject and condition, only if the
        %subject has data for both conditions
        if ~isempty(filtBTC) && ~isempty(filtBTE)
            
            %Extract current subject´s data
            Mbtc=subjectPWallblocks{1}(:,indbtc);
            Mbte=subjectPWallblocks{1}(:,indbte);
            % Mbtc=KEPTtrialsSUB{1,1}(channel,:,indbtc);  d=size(Mbtc);   Mbtc=reshape(Mbtc,d(2),d(3));
            % Mbte=KEPTtrialsSUB{1,1}(channel,:,indbte);  d=size(Mbte);   Mbte=reshape(Mbte,d(2),d(3));
            
            %Compute average value in period previous to number appearance
            baselineActivityValueOT=[baselineActivityValueOT ComputeBaselineActivity(Mbtc,mytimes,baselineActivityInt)]; 
            baselineActivityValueMW=[baselineActivityValueMW ComputeBaselineActivity(Mbte,mytimes,baselineActivityInt)];
            
            %Remove baseline
            MbtcBC=removeBaseline(Mbtc,mytimes,baselineInt);
            MbteBC=removeBaseline(Mbte,mytimes,baselineInt);
            
            %Store data for each subject
            MbtctotSUB{isub,channel}=MbtcBC;
            MbtetotSUB{isub,channel}=MbteBC;
            
            %Merge data of different subjects (not really useful!)
            Mbtctot{1,channel}=[Mbtctot{1,channel} MbtcBC];
            Mbtetot{1,channel}=[Mbtetot{1,channel} MbteBC];
            
            ERsubsOT{1,channel}=[ERsubsOT{1,channel}; mean(MbtcBC,2)'];
            ERsubsMW{1,channel}=[ERsubsMW{1,channel}; mean(MbteBC,2)'];
            
            %ERP plot of all kept trials of the current subject
            figure
            mystr=['Clean trials. Delays: ' num2str(delays(2:end)) '. Subject= ' substr];
            set(gcf,'name',mystr);
            [~, indicesLastTrialBlock]=findBlock(subjectPWallblocks{2},blocksInfo);
            erpimage(subjectPWallblocks{1},[],mytimes,mystr,[],[],'erp','horz',indicesLastTrialBlock);
        end
    end
    
end



%% Statistical analysis
selChannels={'Width'};
indSelChannels=1;

interval=[700 1500]*10^-3; %Interval in which to compare the signals
alpha=0.05;
if indSelChannels ~= 1 %Two ways ANOVA for repeated measures
    [stats, Y, S, F1, F2] = myAnovaTwoWaysRM(MbtctotSUB,MbtetotSUB,selChannels,interval,indSelChannels,mytimes);
else              %paired t-test
    stats= myPairedTtest(MbtctotSUB,MbtetotSUB,interval,mytimes,alpha);
end

%Baseline differences
% [h,p]=ttest(baselineActivityValueOT,baselineActivityValueMW,'Alpha',alpha);

% D=[];
% nselch=length(indSelChannels);
% for sub=1:nsubjs
%    indsub=S==sub;
%    for ch=1:nselch
%        indch=F2==indSelChannels(ch);
%        indc1=F1==1;
%        indc2=F1==2;
%
%        ind1 = indsub & indch & indc1;
%        ind2 = indsub & indch & indc2;
%        d = Y(ind1) - Y(ind2);
%        D=[D d];
%    end
% end

%% ER plot of the two different conditions
X=[-500 2000]*10^-3;
lastSubjectTrial_btc=cumsum(findDims(MbtctotSUB));
lastSubjectTrial_bte=cumsum(findDims(MbtetotSUB));

nt=length(mytimes);
ERPot=zeros(nch,nt);
ERPmw=zeros(nch,nt);
avgsig=zeros(nch,nt);
mydiff=zeros(1,nch);

n1=size(Mbtctot{1,1},2); %Total number of OT trials
n2=size(Mbtetot{1,1},2); %Total number of MW trials
for channel=1:nch
    %     erpOT=mean(Mbtctot{1,channel},2)'; %Erp current electrode
    %     erpMW=mean(Mbtetot{1,channel},2)';
    
    erpOT=mean(ERsubsOT{1,channel},1); %Er signal current channel (average all aubjects)
    erpMW=mean(ERsubsMW{1,channel},1); %Er signal current channel (average all aubjects)
    
    %% Store channel data
    ERPot(channel,:)=erpOT;
    ERPmw(channel,:)=erpMW;
    
    %Compute average signal
    
    %     avgsig(channel,:)=(n1*erpOT + n2*erpMW)/(n1+n2);
    avgsig(channel,:)=(erpOT + erpMW)/2;
    %     if sum(indSelChannels==channel)==1
    %Plot trials of all subjects--------------------------------------
    figure
    %Plot before target correct
    subplot(2,1,1)
    mytitle=['ER Pupil dilation '];
    set(gcf,'name', mytitle,'NumberTitle','off');
    mystr2=['Trials before target correct (ON TASK)'];
    [~,~,~,~,~,erp11]=erpimage(Mbtctot{1,channel},[],mytimes,mystr2,nbOfTrialsToSmooth,[],'erp','horz',lastSubjectTrial_btc); % Image trials in input order
    
    %Plot before target correct
    subplot(2,1,2)
    mystr2=['Trials before target error (MIND WANDERING). '];
    [~,~,~,~,~,erp22]=erpimage(Mbtetot{1,channel},[],mytimes,mystr2,nbOfTrialsToSmooth,[],'erp','horz',lastSubjectTrial_bte); % Image trials in input order
    
    %Plot trials of all subjects ordered by latency--------------------
    figure
    %Plot before target correct
    subplot(2,1,1)
    mytitle=['Latency-ordered ER pupil dilation '];
    set(gcf,'name', mytitle,'NumberTitle','off');
    mystr2=['Trials before target correct (ON TASK)' ];
    [~,~,~,~,~,~,~,~,~,~,~,~,~,sortot]=erpimage(Mbtctot{1,channel},LATENCIESotTot/10^3,mytimes,mystr2,nbOfTrialsToSmooth,[],'erp'); % Image trials in input order
    
    %Plot before target correct
    subplot(2,1,2)
    mystr2=['Trials before target error (MIND WANDERING). '];
    [~,~,~,~,~,~,~,~,~,~,~,~,~,sortmw]=erpimage(Mbtetot{1,channel},LATENCIESmwTot/10^3,mytimes,mystr2,nbOfTrialsToSmooth,[],'erp'); % Image trials in input order
    %% ERPs of the 2 conditions
    
    figure
    mytitle=['ER pupil dilation for MW and OT conditions'];
    set(gcf,'name',mytitle);
    plot(mytimes,avgsig(channel,:),mytimes,erpOT,mytimes,erpMW);
    legend('Grand Average','On task','Mind wandering')
    grid on
    xlim(X);
    
    %     end
    % %     %% Compute mean value in a specific interval
    % %     [mina,ia]=min(abs(mytimes-interval(1)));
    % %     [minb,ib]=min(abs(mytimes-interval(2)));
    % %     amp1=mean(erpOT(ia:ib));
    % %     amp2=mean(erpMW(ia:ib));
    % %     mydiff(channel)=amp1-amp2;
    %
    %     %     title(['Difference in interval I = ( ' num2str(interval(1)) ' , ' num2str(interval(2)) ' )' ]);
end
% kmax=4;
% hmax=8;
% mystr='MW vs OT ERPs - Grand average all subjects';
% myERPplot(avgsig,kmax,hmax,[],allChs,mystr,mytimes,X,ERPot,ERPmw,mydiff);
%
%

%% PLOT BASELINE ACTIVITY VALUE
mystr=['Baseline activity values in the interval ' num2str(baselineActivityInt)];
figure
set(gcf,'name',mystr)
bar([[mean(baselineActivityValueOT) baselineActivityValueOT]' [mean(baselineActivityValueMW) baselineActivityValueMW]']);
grid on
mytick=cell(1,nsubjs);
mytick{1}='Averages';
for isub=1:nsubjs
    csub=subjectsToConsider(isub);
    substr=num2str(csub);
    mytick{1+isub}=['S' substr];
end
set(gca,'XTickLabel',mytick)
title(mystr)

legend('On task','Mind wandering');