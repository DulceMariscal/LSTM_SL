%% SEGMENT EYE TRACKER DATA (EVOKED RESPONSE EXTRACTION)
clc
close all
clear all


%Change directory
cd('C:\Users\Alessandro S\Dropbox\BitBrainShared\BitBrain')

%Select experiment and subject
experiment=2;
subjects=[1:8];

%Define strings useful to upload/download
expstr=num2str(experiment);
current_dir = pwd;

%Load events, logs and experiment parameters
load('eventsDefinition.mat');                      %Contains events' definitions (StartSquare...)
load('logsDefinition.mat');                        %Contains logs' definitions   (data relative to each trial - i.e. number shown...)
expParPath=[pwd '\Storage Data SART experiment\' expstr '\Parameters\param1.mat'];
load(expParPath);

%Parameters
nblocks=5;
etsr=32;
tstep=1/etsr;
nanVal=-1;
nchET=4;
nsubs=length(subjects);
pw=3; %Pupil width
vd=4; %Valid data mask
numBlockTrials=360;
smallwoodPreproc=1;
grandchampPreproc=1;
minBlinkDur=50*10^-3;
maxBlinkDur=1500*10^-3; %However blink dilation is in the range 50-500ms; episodes of eyelid closure of longer duration are called microsleep episodes
%Alignment event-----------------------------------------------------------
selectedEvent=startNumber;
% selectedEvent=keyPress;
%--------------------------------------------------------------------------
SAVEFLAG=0;

%Epoch interval
t1=-1.38; %-1.19;
t2=2;     %1.19

%Parameter artifact removal
blinkGap= 150 *10^-3; %Blinks separated by a gap strictly lesser than blinkGap are merged
perc=0.5;             %Maximum percentage of bad samples in a trial to keep trial's data
differentValues=10;   %Minimum number of different samples in a trial for it to be kept
medianFilterOrder=5;
highCF=10; %Hertz
[bfilt,afilt] = butter(4,[highCF]*2/etsr,'low');
a_velcoeff=1.5; %Coefficient to individuate abnormal pw velocity changes

%Initialize data structures
ETer=cell(1,nblocks-1); %Event related eye tracker data
ETerAllblocks=[];

%Will keep the data for each subject
subjectPW=cell(2,nblocks-1); %Block by block
% subjectPWallblocks=[];       %Concatenated

%Will keep the data of every subject
subjectsPWallblocks=[];      %Concatenated
rawSubjectData=[];
blinksData = repmat(struct('rates',[],'durations',[]),1, nsubs);


f1=figure;
f2=figure;
mycolors={'r','g','c'};

%Start analysis
for sub=1:nsubs
    %Initialize variables
    subjectPWallblocks=cell(1,2);       %Concatenated
    offs=0;
    ETerAllblocks=[];
    blinkRates=[];         blinkRate=[];
    meanBlinkDurations=[]; meanBlinkDuration=[];
    
    %Extract subject identifier
    csub=subjects(sub);
    substr=num2str(csub);
    
    %Selection of the directory from which to upload data
    pathname=[current_dir '\Storage Data SART experiment\' expstr  '\S' substr];
    fname1=strcat('\rawData_s',substr,'.mat');
    fname2=strcat('\Bmat_s',substr,'.mat');
    rawDataFileName=[pathname fname1];
    BmatFileName=[pathname fname2];
    
    %Upload data
    load(rawDataFileName);           %E, L, ET, EEG, (BVP, GSR)    (1 cell for each block)
    clear('EEG','BVP','GSR');
    load(BmatFileName);              %B matrix containts log       (1 cell for each block)
    
    %Segment data
    subjectBlocks=blocksToAverage{csub};
    nbl=length(subjectBlocks);
    for iblock=1:nbl
        %Select current block
        block=subjectBlocks(iblock);
%       block
        %Find trials in which there has been a response
        indNoResponse=B{response_time,block}==-1;
        
        %Clean and extract data
        tobiiData=ET{1,block};
        M=cleanETdata(tobiiData,nanVal); %Average valid rigth-left values
        [blockData,mytimes]=extractERdata(M,E{block},L{block},selectedEvent,etsr,nchET,t1,t2,indNoResponse); %Extract event related data 4 current block
        
        %Store data current block
        ETer{block-1}=blockData;
        ETerAllblocks=cat(3,ETerAllblocks,blockData);
        
        %Extract only needed (width and valid) data
        dim=size(blockData);
        widthData=reshape(blockData(pw,:,:),dim(2),dim(3));
        validData=reshape(blockData(vd,:,:),dim(2),dim(3));
        
        %Remove bad trials (those with almost all nonvalid values or always same value)
        keptWidthTrials=removeBadTrials(widthData,validData,perc,differentValues);
        redValid=validData(:,keptWidthTrials{2});
        redNotValid=~redValid;
        
        %ARTIFACT INTERPOLATION
        %Find invalid values, merge close intervals of invalid values, linear interpolation
        [cleanWidth, newNotValidData]=ETartifactRemover(keptWidthTrials{1},substr,block,blinkGap,etsr,mytimes,redNotValid);
        
        if smallwoodPreproc
            %MEDIAN FILTER (to remove spikes)
            cleanWidth=myMedianFilter(cleanWidth,medianFilterOrder); %Filters along the columns
            %LOW PASS FILTER
            cleanWidth = myFiltfilt(bfilt,afilt,cleanWidth);
            %NORMALIZATION
            cleanWidth = zscore(cleanWidth,[],1);
        end
        
        if grandchampPreproc
            %Compute velocity
            velocities=diff(cleanWidth)./tstep;
            avgVel=mean(velocities(:));
            avgStd=std(velocities(:));
            
            %Find artifacts
            AANonValidData=findVelArtifacts(velocities,keptWidthTrials,avgVel,avgStd,newNotValidData,a_velcoeff); %Valid data after velocity artifacts detection
            
            %Merge and interpolate (Original dataset and new mask of not valid values)
            [cleanWidth, finalNonValidMask,blinksCellarr]=ETartifactRemover(keptWidthTrials{1},substr,block,blinkGap,etsr,mytimes,AANonValidData);
            
            %Find blinks
            nsamples=size(cleanWidth,1);
            [blinkRate,meanBlinkDuration]=findBlinks(blinksCellarr,minBlinkDur,maxBlinkDur,etsr,nsamples); %Should be done on continuous epochs and not on trials!!
            
        end
        
        %Store data current block
        subjectPW{1,block-1}=cleanWidth;
        subjectPW{2,block-1}=keptWidthTrials{2};
        
        subjectPWallblocks{1}=[subjectPWallblocks{1,1}, cleanWidth];
        subjectPWallblocks{2}=[subjectPWallblocks{2}, keptWidthTrials{2} + offs];
        
        blinkRates=[blinkRates blinkRate];
        meanBlinkDurations=[meanBlinkDurations meanBlinkDuration ];
        
        offs = offs + numBlockTrials;
    end
    subjectsPWallblocks=[subjectsPWallblocks, subjectPWallblocks{1}];
    
    %Clean data (interpolate -1 values and delete trials with too many interpolation )
    %Save data
    subjectMeanResponse=mean(subjectPWallblocks{1},2);
    if SAVEFLAG
        if selectedEvent==keyPress
            fname3=strcat('\RESPwidthData_s',substr,'.mat');
        else
            fname3=strcat('\widthData_s',substr,'.mat');
        end
        widthDataFileName=[pathname fname3];
        save(widthDataFileName,'subjectMeanResponse','subjectPWallblocks','mytimes');
    end
    if sub<=3
        %Plot subject data
        figure(f1)
        hold on
        plot(mytimes,subjectMeanResponse,mycolors{sub});
        grid on
        
        %Plot average of old and clean data
        figure(f2)
        subplot(3,1,sub)
        plot(mytimes,subjectMeanResponse,'b')
        hold on
        originalWidth=ETerAllblocks(pw,:,:);
        
        nanOW=originalWidth;
        nanOW(nanOW==-1)=nan;
        valMean=nanmean(nanOW,3);
        plot(mytimes,valMean,'r') %To see how it was before artifact removal
        hold on
        plot(mytimes,mean(originalWidth,3),'g') %To see how it was before artifact removal
    end
    %Clear non needed data
    disp(['Subject num ' substr ' has a total of ' num2str(size(subjectPWallblocks{1},2))  ' trials with valid ET data'])
    disp(['Subject num ' substr ' has a total of ' num2str(sum(sum(subjectPWallblocks{1}==-1))) ' remaining invalid samples '])
    %     subjectPWallblocks=[];
    
    blinksData(sub).rates=blinkRates;
    blinksData(sub).durations=meanBlinkDurations;
end
