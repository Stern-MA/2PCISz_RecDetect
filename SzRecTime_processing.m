%MASTER PROCESSING SCRIPT
%
%
%Input .mat file from Suite2P
%Output is .mat file structure with fields to work with data for anaysis


%Navigate to parent suite2P directory for seizure processing
%Directory Tree ch1 is in suite2P from directory 1 and ch 2 has its own

%Outline
%  0) Define Variables
%  1) Open EEG from EDF and annotation files
%  2) Open F and Fneu values for cells
%  3) Sort the F values into true cells and generate sorting indicies
%  4) Subtract background signal and neuropil contamination
%  5) Determine dF/F
%  6) Generate population mean traces
%  7) Correct for EEG Bias (time shift) between 2p and eeg
%  8) Power Spectral Density Determinations
%  9) Death Determination
% 10) Mean Seziure and CSD Seed Time Determination
% 11) Mean Pre-ictal Spike Seed Time Determiantion ---EEG_PIS_times
% 12) EEG Kernal Generation and Refinement of EEG PIS detection ---EEG_PIS_times2
% 13) First PIS determination and Final PIS SEED adjustment ---EEG_PIS_times2
% 14) Filtering of Individual Traces to 1 Hz
% 15) Individual trace seizure recruitment determination (f') ---RecTimeG_sz2
% 16) Define if a cell is recruited during seizure ---isRecruitedG
% 17) Individual trace PIS recruitment determination (f') ---RecTimeG_PIS.v3
% 18) Define if spike included in analysis (>10% cell rec) ---EEG_PIS_times3, recTimeG_PIS.true
%
%
%
%
%

%%
szID='S75ptz2';
plotON=0;

%% Variables
if strcmp(szID,'S75ptz1')
EEGbias=0.2;
else
EEGbias=0; %S75ptz1:.2 ;if EEG seems shifted relative to the Ca signal(+ is EEG forward/ahead relative to calcium); default is 0
end

fs_eeg=2000;
fs_2p=30;

%Manually constrain Pre-Ictal times
limit1=0; %varaible to set to 1 if need to manually limit pre-ictal period

%% Loader EEG from EDF file
[EEG_all]=ParseEDF('export1.edf','annotations1.txt',fs_eeg,plotON);
EEG=EEG_all.part{2};
EEG_ts=[1/fs_eeg:1/fs_eeg:length(EEG)/fs_eeg];

%% Loader (Green (ch 1))
load('suite2p/plane0/Fall.mat')
F_G=F; Fneu_G=Fneu; iscell_G=iscell; ops_G=ops; spks_G=spks; stat_G=stat;
clear F Fneu iscell ops redcell spks stat


% %% Loader (Red (ch 2))
% load('chan2/suite2p/plane0/Fall.mat')
% F_R=F; Fneu_R=Fneu; iscell_R=iscell; ops_R=ops; spks_R=spks; stat_R=stat;
% clear F Fneu iscell ops redcell spks stat


%% Select Only Selected ROIs
% GREEN CHANNEL
F1G=nan(sum(iscell_G(:,1)),size(F_G,2)-1); %generate matrix to fill with only cell ROIs F
Fneu1G=nan(sum(iscell_G(:,1)),size(Fneu_G,2)-1); %generate matrix to fill with only cell ROIs neu
xCoorG=nan(sum(iscell_G(:,1)),1); %generate matrix to fill with x coordinate values for only cells green
yCoorG=nan(sum(iscell_G(:,1)),1); %generate matrix to fill with y coordinate values for only cells green
iscellIG=nan(size(iscell_G)); %generate matrix for a key of the iscell ROIs' original ROI ID numbers

jj=1;%counter total will equal number of true cells
for ii=1:size(F_G,1)%number of ROIs (both true and not true cells)
    if iscell_G(ii)==1
        F1G(jj,:)=F_G(ii,1:end-1); %populate new cell only matrix with traces
        Fneu1G(jj,:)=Fneu_G(ii,1:end-1);
        xCoorG(jj)=stat_G{ii}.med(2); %MODIFIED TO FIX X AND Y FLIP
        yCoorG(jj)=513-stat_G{ii}.med(1);%Modified to flip the y axis to match the images without needign to change axis orientation 
        iscellIG(ii,:)=[ii,jj];
        jj=jj+1;
    end
end
clear ii
clear jj

[~,xCoorIG]=sort(xCoorG); %get index values for a sorted dataset with respect x coordinate
%can sort cells by this index using F(xCorrI,:); %reorder F based upon x
%coordinant
[~,yCoorIG]=sort(yCoorG); %same for y

% % RED CHANNEL
% F1R=nan(sum(iscell_R(:,1)),size(F_R,2)-1);
% Fneu1R=nan(sum(iscell_R(:,1)),size(Fneu_R,2)-1);
% xCoorR=nan(sum(iscell_R(:,1)),1);
% yCoorR=nan(sum(iscell_R(:,1)),1);
% iscellIR=nan(size(iscell_R));
% 
% jj=1;
% for ii=1:size(F_R,1)
%     if iscell_R(ii)==1
%         F1R(jj,:)=F_R(ii,1:end-1);
%         Fneu1R(jj,:)=Fneu_R(ii,1:end-1);
%         xCoorR(jj)=stat_R{ii}.med(2);
%         yCoorR(jj)=513-stat_R{ii}.med(1);
%         iscellIR(ii,:)=[ii,jj];
%         jj=jj+1;
%     end
% end
% clear ii
% clear jj
% 
% [~,xCoorIR]=sort(xCoorR);
% [~,yCoorIR]=sort(yCoorR);

%Write a timestamp for plotting
F1_ts=[1/fs_2p:1/fs_2p:length(F1G)/fs_2p];


%% Baseline shift and neuropil subtraction
%must subtract the baseline or at least get transients to roughly same
%positive range
%Fb: indicates background subtracted; Fc neuropil corrected
Fb1G=F1G-min(min(Fneu1G));%shifts F trace by the global Fneu minimum (assumed as background)
% Fb1R=F1R-min(min(Fneu1R));
Fneub1G=Fneu1G-min(min(Fneu1G));% must shift the neuropil too so it can be properly scaled
% Fneub1R=Fneu1R-min(min(Fneu1R));

Fc1G=Fb1G-(0.7.*Fneub1G);%gives neuropil subtracted data (serves as pure F)
% Fc1R=Fb1R-(0.7.*Fneub1R);


%% Generate Mean Population Traces
F1G_Mean=mean(F1G,1);%population mean of raw F
% F1R_Mean=mean(F1R,1);
Fneu1G_Mean=mean(Fneu1G,1);%population mean of raw F
% Fneu1R_Mean=mean(Fneu1R,1);

Fb1G_Mean=mean(Fb1G,1);%cell population mean of background subtracted data 
% Fb1R_Mean=mean(Fb1R,1);
Fneub1G_Mean=mean(Fneub1G,1);%neuropil population mean of background subtracted data 
% Fneub1R_Mean=mean(Fneub1R,1);

Fc1G_Mean=mean(Fc1G,1);%population mean of background and neuropil subtracted data 
% Fc1R_Mean=mean(Fc1R,1);


%% Generate dF/F
%(F-F0)/F0 with F0 being the firat 30 seconds of data
F0sec=30; %time period for F0 in seconds
Fc1Gdff=(Fc1G-mean(Fc1G(:,1:round(F0sec*fs_2p)),2))./mean(Fc1G(:,1:round(F0sec*fs_2p)),2);%dF/F
% Fc1Rdff=(Fc1R-mean(Fc1R(:,1:round(F0sec*fs_2p)),2))./mean(Fc1R(:,1:round(F0sec*fs_2p)),2);
Fc1Gdff_Mean=mean(Fc1Gdff,1);%population mean of dF/F
% Fc1Rdff_Mean=mean(Fc1Rdff,1);
Fc1Gdff_std=std(Fc1Gdff);%standard deviation of the population mean of dF/F
% Fc1Rdff_std=std(Fc1Rdff);

%background subtracted soma F
Fb1Gdff=(Fb1G-mean(Fb1G(:,1:round(F0sec*fs_2p)),2))./mean(Fb1G(:,1:round(F0sec*fs_2p)),2);
% Fb1Rdff=(Fb1R-mean(Fb1R(:,1:round(F0sec*fs_2p)),2))./mean(Fb1R(:,1:round(F0sec*fs_2p)),2);
Fb1Gdff_Mean=mean(Fb1Gdff,1);
% Fb1Rdff_Mean=mean(Fb1Rdff,1);
Fb1Gdff_std=std(Fb1Gdff);
% Fb1Rdff_std=std(Fb1Rdff);

%background subtracted neuropil Fneu
Fneub1Gdff=(Fneub1G-mean(Fneub1G(:,1:round(F0sec*fs_2p)),2))./mean(Fneub1G(:,1:round(F0sec*fs_2p)),2);
% Fneub1Rdff=(Fneub1R-mean(Fneub1R(:,1:round(F0sec*fs_2p)),2))./mean(Fneub1R(:,1:round(F0sec*fs_2p)),2);
Fneub1Gdff_Mean=mean(Fneub1Gdff,1);
% Fneub1Rdff_Mean=mean(Fneub1Rdff,1);
Fneub1Gdff_std=std(Fneub1Gdff);
% Fneub1Rdff_std=std(Fneub1Rdff);


%% Filter (low pass) version 1: on dffs calculated from raw unfiltered data
lpfilt=5;
%hpfilt=.1;
%green
for ii=1:size(F1G,1)
    Fb1Gdff_flt1(ii,:)=lofi(Fb1Gdff(ii,:),10^6/fs_2p,lpfilt);
    Fneub1Gdff_flt1(ii,:)=lofi(Fneub1Gdff(ii,:),10^6/fs_2p,lpfilt);
    Fc1Gdff_flt1(ii,:)=lofi(Fc1Gdff(ii,:),10^6/fs_2p,lpfilt);
    %Fc1Gdff_flt1hp(ii,:)=hifi(Fc1Gdff_flt1(ii,:),10^6/F1_fs,hpfilt);
end

% %red
% for jj=1:size(F1R,1)
%     Fb1Rdff_flt1(jj,:)=lofi(Fb1Rdff(jj,:),10^6/fs_2p,lpfilt);
%     Fneub1Rdff_flt1(jj,:)=lofi(Fneub1Rdff(jj,:),10^6/fs_2p,lpfilt);
%     Fc1Rdff_flt1(jj,:)=lofi(Fc1Rdff(jj,:),10^6/fs_2p,lpfilt);
%     %Fc1Rdff_flt1hp(jj,:)=hifi(Fc1Rdff_flt1(jj,:),10^6/F1_fs,hpfilt);
% end

clear lpfilt hpfilt
clear ii jj

% Generate Means for filtered data version 1
Fb1Gdff_flt1_Mean=mean(Fb1Gdff_flt1,1);
% Fb1Rdff_flt1_Mean=mean(Fb1Rdff_flt1,1);
Fneub1Gdff_flt1_Mean=mean(Fneub1Gdff_flt1,1);
% Fneub1Rdff_flt1_Mean=mean(Fneub1Rdff_flt1,1);
Fc1Gdff_flt1_Mean=mean(Fc1Gdff_flt1,1);
% Fc1Rdff_flt1_Mean=mean(Fc1Rdff_flt1,1);

%Fc1Gdff_flt1hp_Mean=mean(Fc1Gdff_flt1hp,1);
%Fc1Rdff_flt1hp_Mean=mean(Fc1Rdff_flt1hp,1);


%% Correction of EEG bias
if EEGbias>0
    EEGa=[EEG(fs_eeg*EEGbias+1:end),zeros(1,fs_eeg*EEGbias)];
    EEG=EEGa;
elseif EEGbias<0
    EEGa=[zeros(1,fs_eeg*EEGbias),EEG(1:end-fs_eeg*EEGbias)];
    EEG=EEGa;
end
EEGb=circshift(EEG,-fs_eeg*EEGbias);

clear EEGa EEGb


%% PSD initial processing for feature detection
%used for finding pre-ictal spikes and also death following seizure death

PSDbin=1024;
[~,PSD_F]=pwelch(EEG,[],[],PSDbin,fs_eeg);
PSDwin=1;%in seconds

EEGstart=[1:PSDwin*fs_eeg/2:length(EEG)];%sliding window width 2x PSDwin shifting by 1/2 PSDwin
PSD_P=nan(length(PSD_F),length(EEGstart));
for ii = 1:numel(EEGstart)-2
    EEG1=EEG(EEGstart(ii):EEGstart(ii)+PSDwin*fs_eeg-1);
    [PSD_P(:,ii),~]=pwelch(EEG1,[],[],PSDbin,fs_eeg);
    clear EEG1
end

PSD_totP=sum(PSD_P,1);
PSD_theta=sum(PSD_P(and(PSD_F>3,PSD_F<15),:),1);
PSD_gammaL=sum(PSD_P(and(PSD_F>20,PSD_F<55),:),1);
PSD_100plus=sum(PSD_P(PSD_F>=100 & PSD_F<250 ,:),1);
PSD_100less=sum(PSD_P(PSD_F<100,:),1);
Ratio100=PSD_100plus./PSD_100less;
RA_Ratio100=SternRollAvg(Ratio100,60/PSDwin); % +/- 30 second window roll avg (time step here is 0.5s)


%% Deriving recording Features (Death)
if max(RA_Ratio100) > 1
    death=1;
    deathTime=find(Ratio100>1,1)*PSDwin/2;
    deathTime_eeg=deathTime*fs_eeg;
    deathTime_2p=deathTime*fs_2p;
else
    death=0;
    [deathTime,deathTime_2p,deathTime_eeg]=deal(NaN); 
end


%% Mean Trace Seizure and CSD Time Seeds (max first derivative block method)
%using neuropil signal

%Mean Trace Filter
lpfilt=1;

%green
FGdff_MeanSm=lofi(Fneub1Gdff_Mean,10^6/fs_2p,lpfilt); %smoothed

%red
% FRdff_MeanSm=lofi(Fneub1Rdff_Mean,10^6/fs_2p,lpfilt); %smoothed

clear lpfilt

% determine all the points that are greater than the average of the min and
% max values of the trace

if death==1 %to account for an animal dying the data 1 min beyond death is excluded (1 min is to account for hypoxic wave to propogate uncertatinty in exact time of death)
    postDeathT=deathTime_2p+fs_2p*60;
    temp1G=FGdff_MeanSm(1:postDeathT)>=mean([min(FGdff_MeanSm(1:postDeathT)),max(FGdff_MeanSm(1:postDeathT))]);
%     temp1R=FRdff_MeanSm(1:postDeathT)>=mean([min(FRdff_MeanSm(1:postDeathT)),max(FRdff_MeanSm(1:postDeathT))]);
else
    temp1G=FGdff_MeanSm>=mean([min(FGdff_MeanSm),max(FGdff_MeanSm)]);
%     temp1R=FRdff_MeanSm>=mean([min(FRdff_MeanSm),max(FRdff_MeanSm)]);   
end

% find block bounds
temp2G=diff(temp1G,1,2);
% temp2R=diff(temp1R,1,2);

temp3G=cell(1,2);
% temp3R=cell(1,2);
temp3G{1}=find(temp2G==1);
temp3G{2}=find(temp2G==-1);
% temp3R{1}=find(temp2R==1);
% temp3R{2}=find(temp2R==-1);

%integrate over blocks
temp4G=cell(1,1);
% temp4R=cell(1,1);
if temp3G{1}(1)>temp3G{2}(1)%correct for first slope value index being negative
    temp3G{2}=temp3G{2}(2:end);
end

if temp3G{1}(end)>temp3G{2}(end)%correct for extra final slope value being positive
    temp3G{1}=temp3G{1}(1:end-1);
end

for jj=1:length(temp3G{1}) %find the area under the curve of each region
    temp4G{1}(jj)=trapz(FGdff_MeanSm(1,[temp3G{1}(jj):temp3G{2}(jj)]));
end
clear jj

% if temp3R{1}(1)>temp3R{2}(1)%correct for first slope value index being negative
%     temp3R{2}=temp3R{2}(2:end);
% end
% 
% if temp3R{1}(end)>temp3R{2}(end)%correct for extra final slope value being positive
%     temp3R{1}=temp3R{1}(1:end-1);
% end
% 
% for jj=1:length(temp3R{1}) %find the area under the filtered curve of each positive region
%     temp4R{1}(jj)=trapz(FRdff_MeanSm(1,[temp3R{1}(jj):temp3R{2}(jj)]));
% end

clear ii jj


% find maximum block and pull boundary indices of that block
% find each largest value in each cell array and then use this to index into temp3 to find the index (time point) of the seizure in the origional trace
temp5G=cell(1,2);%gives the block in each cell that is seizure
% temp5R=cell(1,2);%gives the block in each cell that is seizure

[temp5G{1},temp5G{2}] = sort(temp4G{1},2,'descend');%sort blocks by size (block index is cell 2)

%determine the block bounds of the seizure and CSD blocks
MeanSzRecTimeG = min(temp3G{1}(temp5G{2}(1:2))+1)/fs_2p;%defines the Sz block as the first block of the top 2 blocks; get the first bound of that block which is the max recruitment time (add 1 due to use of diff fn)
MeanCSDRecTimeG = max(temp3G{1}(temp5G{2}(1:2))+1)/fs_2p; %defines the CSD block as the last block of the top 2 blocks
%MeanCSDRecTimeG = max((temp3G{1}(temp5G{2}(1:ceil(0.5*length(temp5G{2}))))+1)/fs_2p); %defines the CSD block as the last block of the top 20% of blocks

% [temp5R{1},temp5R{2}] = sort(temp4R{1},2,'descend');
% 
% %determine the block bounds of the seizure and CSD blocks
% MeanSzRecTimeR = min(temp3R{1}(temp5R{2}(1:2))+1)/fs_2p;
% MeanCSDRecTimeR = max(temp3R{1}(temp5R{2}(1:2))+1)/fs_2p;
% %MeanCSDRecTimeR = max((temp3R{1}(temp5R{2}(1:ceil(0.5*length(temp5R{2}))))+1)/fs_2p);

newMaxG = max(FGdff_MeanSm((MeanSzRecTimeG-2)*fs_2p:(MeanSzRecTimeG+5)*fs_2p));
% newMaxR = max(FRdff_MeanSm((MeanSzRecTimeR-2)*fs_2p:(MeanSzRecTimeR+5)*fs_2p));

if death==1 %to account for an animal dying the data 1 min beyond death is excluded (1 min is to account for hypoxic wave to propogate uncertatinty in exact time of death)
    postDeathT=deathTime_2p+fs_2p*60;
    temp1G=FGdff_MeanSm(1:postDeathT)>=mean([min(FGdff_MeanSm(1:postDeathT)),newMaxG]);
%     temp1R=FRdff_MeanSm(1:postDeathT)>=mean([min(FRdff_MeanSm(1:postDeathT)),newMaxR]);
else
    temp1G=FGdff_MeanSm>=mean([min(FGdff_MeanSm),newMaxG]);
%     temp1R=FRdff_MeanSm>=mean([min(FRdff_MeanSm),newMaxR]);   
end

% find block bounds
temp2G=diff(temp1G,1,2);
% temp2R=diff(temp1R,1,2);

temp3G=cell(1,2);
% temp3R=cell(1,2);
temp3G{1}=find(temp2G==1);
temp3G{2}=find(temp2G==-1);
% temp3R{1}=find(temp2R==1);
% temp3R{2}=find(temp2R==-1);

%integrate over blocks
temp4G=cell(1,1);
% temp4R=cell(1,1);
if temp3G{1}(1)>temp3G{2}(1)%correct for first slope value index being negative
    temp3G{2}=temp3G{2}(2:end);
end

if temp3G{1}(end)>temp3G{2}(end)%correct for extra final slope value being positive
    temp3G{1}=temp3G{1}(1:end-1);
end

for jj=1:length(temp3G{1}) %find the area under the curve of each region
    temp4G{1}(jj)=trapz(FGdff_MeanSm(1,[temp3G{1}(jj):temp3G{2}(jj)]));
end
clear jj

% if temp3R{1}(1)>temp3R{2}(1)%correct for first slope value index being negative
%     temp3R{2}=temp3R{2}(2:end);
% end
% 
% if temp3R{1}(end)>temp3R{2}(end)%correct for extra final slope value being positive
%     temp3R{1}=temp3R{1}(1:end-1);
% end
% 
% for jj=1:length(temp3R{1}) %find the area under the filtered curve of each positive region
%     temp4R{1}(jj)=trapz(FRdff_MeanSm(1,[temp3R{1}(jj):temp3R{2}(jj)]));
% end

clear ii jj


% find maximum block and pull boundary indices of that block
% find each largest value in each cell array and then use this to index into temp3 to find the index (time point) of the seizure in the origional trace
temp5G=cell(1,2);%gives the block in each cell that is seizure
% temp5R=cell(1,2);%gives the block in each cell that is seizure

[temp5G{1},temp5G{2}] = sort(temp4G{1},2,'descend');%sort blocks by size (block index is cell 2)

%determine the block bounds of the seizure
MeanSzRecTimeG = min(temp3G{1}(temp5G{2}(1:2))+1)/fs_2p;%defines the Sz block as the first block of the top 2 blocks; get the first bound of that block which is the max recruitment time (add 1 due to use of diff fn)

% [temp5R{1},temp5R{2}] = sort(temp4R{1},2,'descend');
% 
% %determine the block bounds of the seizure
% MeanSzRecTimeR = min(temp3R{1}(temp5R{2}(1:2))+1)/fs_2p;

%% Preictal Mean RecTime Generation

PreSpkRecTimeG=MeanRecTimes3pre(FGdff_MeanSm,fs_2p,15);
% PreSpkRecTimeR=MeanRecTimes3pre(FRdff_MeanSm,fs_2p,15);

%calculate an integral of the time 1 second following recruitment to filter
%out false positives
for ii=1:length(PreSpkRecTimeG{1,1})
    tempTrace=FGdff_MeanSm(PreSpkRecTimeG{1,1}(ii)*fs_2p:(PreSpkRecTimeG{1,1}(ii)+1)*fs_2p);
    valPIS1sG(ii)=trapz(tempTrace-min(tempTrace));
    clear tempTrace
end
clear ii

% for ii=1:length(PreSpkRecTimeR{1,1})
%     tempTrace=FRdff_MeanSm(PreSpkRecTimeR{1,1}(ii)*fs_2p:(PreSpkRecTimeR{1,1}(ii)+1)*fs_2p);
%     valPIS1sR(ii)=trapz(tempTrace-min(tempTrace));
%     clear tempTrace
% end
% clear ii

CaisLrgG=valPIS1sG>std(FGdff_MeanSm(1:MeanSzRecTimeG*fs_2p))*fs_2p/3/0.3989; %gaussian, height=std of preictal period, sigma=1/3 second in frames (1s width at 1.5 sigma)
% CaisLrgR=valPIS1sR>std(FRdff_MeanSm(1:MeanSzRecTimeR*fs_2p))*fs_2p/3/0.3989;


for ii=1:4
     PreSpkRecTimeG{1,ii}=PreSpkRecTimeG{1,ii}(CaisLrgG);
%      PreSpkRecTimeR{1,ii}=PreSpkRecTimeR{1,ii}(CaisLrgR);
end

S_th = 4; %power threshold based on median for PIS detection

winW = .7; %seconds for spike detection window width

CaDetectWin=.2;

EEG_spk_prelim = BO_IIS_Detect_Power_Rolston_v5(EEG',S_th,fs_eeg,winW,0); %use v4 to not use theta gamma power ratio component
EEG_spkTimes_prelim_old=EEG_spk_prelim.time_Sp_final{1};%grab eeg spikes indexed by largest change in feature point
EEG_spkTimes_prelim=EEG_spk_prelim.time_Sp_minPeak{1};%grab all possible EEG spike times indexed by the min peak (spike) time point

%find Ca spikes in mean traces that coorepsond to EEG spikes
DiffPreEEG_G=abs(PreSpkRecTimeG{1,1}-(EEG_spkTimes_prelim'));%generates a Ca spk by EEG spk matrix of time diff
% DiffPreEEG_R=abs(PreSpkRecTimeR{1,1}-(EEG_spkTimes_prelim'));
MeanPreIsRecG=min(DiffPreEEG_G,[],1)<CaDetectWin & PreSpkRecTimeG{1,1}<MeanSzRecTimeG;%logical indexing across Ca spk min to indicate if recruited
% MeanPreIsRecR=min(DiffPreEEG_R,[],1)<CaDetectWin & PreSpkRecTimeR{1,1}<MeanSzRecTimeR;
MeanCaG_PIS_times=PreSpkRecTimeG{1,1}(MeanPreIsRecG);
% MeanCaR_PIS_times=PreSpkRecTimeR{1,1}(MeanPreIsRecR);

%refine EEG spike times based on Ca recruitment and the seizure recruitment
%time (smaller of the red or green channel)...may need to make it the
%larger of the two instead, REVISIT
EEGisTrue=and(min(DiffPreEEG_G,[],2)<CaDetectWin | min(DiffPreEEG_R,[],2)<CaDetectWin , EEG_spkTimes_prelim'<MeanSzRecTimeG);%final EEG spikes times i.e. correspond to at least one channel of Ca recruitment
% EEGisTrue=and(min(DiffPreEEG_G,[],2)<CaDetectWin | min(DiffPreEEG_R,[],2)<CaDetectWin , EEG_spkTimes_prelim'<max(MeanSzRecTimeR,MeanSzRecTimeG));%final EEG spikes times i.e. correspond to at least one channel of Ca recruitment

EEG_PIS_times=EEG_spkTimes_prelim(EEGisTrue);%pre-ictal spike times 


%% EEG spike wave kernel and PIS refinement
EEGswd_KernSet=nan(numel(EEG_PIS_times),1*fs_eeg);
for ii=1:numel(EEG_PIS_times)
    EEGswd_KernSet(ii,:)=EEG((EEG_PIS_times(ii)-.5)*fs_eeg+1:(EEG_PIS_times(ii)+.5)*fs_eeg);
end

EEGswd_kern=mean(EEGswd_KernSet,1);

% Convolve template with eeg
EEGswdConv=conv(flip(EEGswd_kern), EEG);
[~,EEGswd_spikeTimes]=findpeaks(EEGswdConv(1:MeanSzRecTimeG*fs_eeg),'MinPeakHeight',5*std(EEGswdConv(1:MeanSzRecTimeG*fs_eeg)),'MinPeakDistance',fs_eeg);
EEGswd_spikeTimes=EEGswd_spikeTimes/fs_eeg-0.5;

%combine the new detected spikes into the calcium contstrained spikes
EEGswdPISdiff=EEG_PIS_times-EEGswd_spikeTimes';
EEGswdPISmin=min(abs(EEGswdPISdiff),[],2);
EEG_PIS_new=EEGswd_spikeTimes(EEGswdPISmin>1);
EEG_PIS_times2=sort([EEG_PIS_times,EEG_PIS_new]);
EEG_PIS_times2=EEG_PIS_times2(EEG_PIS_times2<MeanSzRecTimeG; %make sure to grab sentinal spike but not seizure
% EEG_PIS_times2=EEG_PIS_times2(EEG_PIS_times2<max(MeanSzRecTimeG,MeanSzRecTimeR)); %make sure to grab sentinal spike but not seizure

if plotON==1
figure
for ii=1:size(EEGswd_KernSet,1)
    subplot(ceil(sqrt(size(EEGswd_KernSet,1))),ceil(sqrt(size(EEGswd_KernSet,1))),ii)
    plot(EEGswd_KernSet(ii,:))
    title(ii)
end
suptitle('Ca EEG Concordant Spks')
clear ii

figure
plot(EEGswd_kern)
suptitle('EEG PIS Kernel')
end

%% interspike interval and first spike determination
%looks for a rolling average of 10 to designate consistent spiking but then
%looks for the first spike in this train that has an ISI within 15 seconds
%as an upper limit on the allowable ISI
ISI=diff(EEG_PIS_times2);
ISI_RAwin=5;
RAtolerance=10;
IndivTolerance=15;
if numel(ISI)<ISI_RAwin
    firstSpkI=find(ISI<IndivTolerance,1);%if less than 5 spikes allow a slightly larger tolerance for ISI to grab spieks to maximize spikes collected
else
    ISI_RA=SternRollAvg(ISI,ISI_RAwin);
    
    ISI_RAthresh=find(ISI_RA<RAtolerance,1);%finds the spike at the beginning of a train where the roll avg is <10
    if ISI_RAthresh>ISI_RAwin
        firstSpkI=find(ISI(ISI_RAthresh-ISI_RAwin:ISI_RAthresh)<IndivTolerance,1);%finds the first spike with an ISI <10 in the first window where the roll avg is <10.
        firstSpkI=firstSpkI+ISI_RAthresh-ISI_RAwin-1;%neg 1 is needed for fence post (index 1 needs to coorespond to first spk in range considered)
    else
        firstSpkI=find(ISI(1:ISI_RAthresh+ISI_RAwin)<IndivTolerance,1);
    end
end

firstSpkT=EEG_PIS_times2(firstSpkI);

EEG_PIS_times2=EEG_PIS_times2(firstSpkI:end);

    
%% Filtering for Individual Trace Recruitment Time Determination
% filtering to smooth for recruitment detection
lpfilt=1;
%hpfilt=0.1;
Fc1Gdff_flt2=zeros(size(Fc1Gdff));
%Fc1Gdff_flt2_hp=zeros(size(Fc1Gdff_flt2));
% Fc1Rdff_flt2=zeros(size(Fc1Rdff));
%Fc1Rdff_flt2_hp=zeros(size(Fc1Rdff_flt2));
%green
for ii=1:size(Fc1Gdff_flt1,1)
    Fc1Gdff_flt2(ii,:)=lofi(Fc1Gdff(ii,:),10^6/fs_2p,lpfilt);
    %Fc1Gdff_flt2_hp(ii,:)=hifi(Fc1Gdff_flt2(ii,:),10^6/fs_2p,hpfilt);
end

% %red
% for jj=1:size(Fc1Rdff_flt1,1)
%     Fc1Rdff_flt2(jj,:)=lofi(Fc1Rdff(jj,:),10^6/fs_2p,lpfilt);
%     %Fc1Rdff_flt2_hp(jj,:)=hifi(Fc1Rdff_flt2(jj,:),10^6/fs_2p,hpfilt);
% end

clear lpfilt
clear ii jj


%% Individual seizure recruitment detection (Version 3)
%Did sigma 3s origionally but that grabs the sentinal spike more often in
%somw of the data so using a narrower signam is better for elimianting that
%(e.g. s69), but in some of the slower seizures that grabs the incorrect
%times more often hence needing a larger sigma (e.g. s75)
%s61...1s; s71...either; s75-1...2s; s75-2...
RecTimeG_sz3_struct=IndivRecTimes3(Fc1Gdff_flt2,MeanSzRecTimeG,fs_2p,1); %selecting a 1s sigma for gaussian to prevent sentinal spike from being grabbed and also if the average speed of the seizure is 200um/s then it should not take more than 1.5 seconds to cross the field so a 1s sigma gives us a 3 second window at at least 33% (1.5sigma) amplification of feature occurs
% RecTimeR_sz3_struct=IndivRecTimes3(Fc1Rdff_flt2,MeanSzRecTimeR,fs_2p,1);
RecTimeG_sz2=RecTimeG_sz3_struct.time;
% RecTimeR_sz2=RecTimeR_sz3_struct.time;


%% Recruited cell or not
%trace needs to have atleast a 20% increase in signal during seizure and the
%seizure detected time needs to be determined to be within 3 seconds of
%the seed time (if average speed is 100um/s and field is <300um all values
%should be within 3 seconds of mean giving 6 second window)

tprd=10;%time period (s)
gthresh=0.2; %green threshold typical G threshold .3
% rthresh=0.2; %red threshold typical R threshold .3
preSzMeanDffG=mean(Fc1Gdff_flt2(:,(MeanSzRecTimeG-tprd)*fs_2p:MeanSzRecTimeG*fs_2p),2);
postSzMeanDffG=mean(Fc1Gdff_flt2(:,MeanSzRecTimeG*fs_2p:(MeanSzRecTimeG+tprd)*fs_2p),2);
SzMeanDffDiffG=postSzMeanDffG-preSzMeanDffG;
isRecruitedG=and(SzMeanDffDiffG>gthresh,abs(RecTimeG_sz2-MeanSzRecTimeG)<3); 
%disp(sum(isRecruitedG))
% 
% preSzMeanDffR=mean(Fc1Rdff_flt2(:,(MeanSzRecTimeR-tprd)*fs_2p:MeanSzRecTimeR*fs_2p),2);
% postSzMeanDffR=mean(Fc1Rdff_flt2(:,MeanSzRecTimeR*fs_2p:(MeanSzRecTimeR+tprd)*fs_2p),2);
% SzMeanDffDiffR=postSzMeanDffR-preSzMeanDffR;
% isRecruitedR=and(SzMeanDffDiffR>rthresh,abs(RecTimeR_sz2-MeanSzRecTimeR)<3);
% %disp(sum(isRecruitedR))
clear tprd gthresh rthresh


%% Individual CSD / Death recruitment detection
%generate signal to noise ratio for weighing the traces
sz_snrG=nan(size(RecTimeG_sz2));
% sz_snrR=nan(size(RecTimeR_sz2));
for ii=1:length(sz_snrG)
     SzRange=[(RecTimeG_sz2(ii)-3)*fs_2p:(RecTimeG_sz2(ii)+3)*fs_2p];
     sz_snrG(ii)=max(Fc1Gdff_flt1(ii,SzRange))/std(Fc1Gdff_flt1(ii,:));
end
clear ii SzRange
% 
% for ii=1:length(sz_snrR)
%      SzRange=[(RecTimeR_sz2(ii)-3)*fs_2p:(RecTimeR_sz2(ii)+3)*fs_2p];
%      sz_snrR(ii)=max(Fc1Rdff_flt1(ii,SzRange))/std(Fc1Rdff_flt1(ii,:));
% end
% clear ii SzRange

% RecTimeG_csd_struct=IndivRecTimes2IP(Fc1Gdff_flt2,MeanCSDRecTimeG,fs_2p,mean(std(Fc1Gdff_flt2))*20);
% RecTimeR_csd_struct=IndivRecTimes2IP(Fc1Rdff_flt2,MeanCSDRecTimeR,fs_2p,mean(std(Fc1Rdff_flt2))*20);
% RecTimeG_csd=RecTimeG_csd_struct.time;
% RecTimeR_csd=RecTimeR_csd_struct.time;
% isRecruitedG_csd=and(RecTimeG_csd_struct.value./sz_snrG>0.2,abs(RecTimeG_csd-MeanCSDRecTimeG)<5);
% isRecruitedR_csd=and(RecTimeR_csd_struct.value./sz_snrR>0.2,abs(RecTimeR_csd-MeanCSDRecTimeR)<5);

RecTimeG_csd_struct=IndivRecTimes3(Fc1Gdff_flt2,MeanCSDRecTimeG,fs_2p,5);%average speed of 40um/s gives means it takes 7s to cross 280um feidl of view, therefore a sigma of 5s gives us a 10s window over which 
% RecTimeR_csd_struct=IndivRecTimes3(Fc1Rdff_flt2,MeanCSDRecTimeR,fs_2p,5);
RecTimeG_csd=RecTimeG_csd_struct.time;
% RecTimeR_csd=RecTimeR_csd_struct.time;
%isRecruitedG_csd=and(RecTimeG_csd_struct.value./sz_snrG>0.2,abs(RecTimeG_csd-MeanCSDRecTimeG)<10);
%isRecruitedR_csd=and(RecTimeR_csd_struct.value./sz_snrR>0.2,abs(RecTimeR_csd-MeanCSDRecTimeR)<10);

%% Recruited cell or not (CSD)
%trace needs to have atleast a 20% increase in signal during CSD relative to before and the
%CSD detected time needs to be determined to be within 10 seconds of
%the seed

tprd=5;%time period (s)
gthresh=0.2; 
% rthresh=0.2; 
CSDMeanDffDiffG=nan(size(RecTimeG_csd));
for ii=1:length(RecTimeG_csd)
    preCSDrange=floor([(RecTimeG_csd(ii)-tprd)*fs_2p:RecTimeG_csd(ii)*fs_2p]);
    postCSDrange=floor([RecTimeG_csd(ii)*fs_2p:(RecTimeG_csd(ii)+tprd)*fs_2p]);
    preCSDMeanDff=mean(Fc1Gdff_flt2(ii,preCSDrange),2);
    postCSDMeanDff=mean(Fc1Gdff_flt2(ii,postCSDrange),2);
    CSDMeanDffDiffG(ii)=postCSDMeanDff-preCSDMeanDff;
end
isRecruitedG_csd=and(CSDMeanDffDiffG>gthresh,abs(RecTimeG_csd-MeanCSDRecTimeG)<10); 
clear preCSDrange postCSDrange preCSDMeanDff postCSDMeanDff
% 
% CSDMeanDffDiffR=nan(size(RecTimeR_csd));
% for ii=1:length(RecTimeR_csd)
%     preCSDrange=floor([(RecTimeR_csd(ii)-tprd)*fs_2p:RecTimeR_csd(ii)*fs_2p]);
%     postCSDrange=floor([RecTimeR_csd(ii)*fs_2p:(RecTimeR_csd(ii)+tprd)*fs_2p]);
%     preCSDMeanDff=mean(Fc1Rdff_flt2(ii,preCSDrange),2);
%     postCSDMeanDff=mean(Fc1Rdff_flt2(ii,postCSDrange),2);
%     CSDMeanDffDiffR(ii)=postCSDMeanDff-preCSDMeanDff;
% end
% isRecruitedR_csd=and(CSDMeanDffDiffR>rthresh,abs(RecTimeR_csd-MeanCSDRecTimeR)<10);
% clear preCSDrange postCSDrange preCSDMeanDff postCSDMeanDff
clear tprd gthresh rthresh


%% Individual Cell PIS Recruitment Time Detection
%taking all PIS including sentintal spike, and using f'-steepest slope for
%time indexing
tempLimSz=MeanSzRecTimeG;
% tempLimSz=max(MeanSzRecTimeG,MeanSzRecTimeR);
RecTimeG_PIS=IndivRecTimes3pre(Fc1Gdff_flt2(:,1:tempLimSz*fs_2p),EEG_PIS_times2,fs_2p,20,1,4);%fed the EEG PIS times
% RecTimeR_PIS=IndivRecTimes3pre(Fc1Rdff_flt2(:,1:tempLimSz*fs_2p),EEG_PIS_times2,fs_2p,20,1,4);
clear tempLimSz


%% Refine Spikes with min 10% recruitment of cells in each channel
   %And remove spikes where the mean is clearly before or after the EEG
   %event window of 0.75 seconds (window used for spike seperation)
PISisRecG=nansum(~isnan(RecTimeG_PIS.v3),1)/size(RecTimeG_PIS.v3,1)>.1;
% PISisRecR=nansum(~isnan(RecTimeR_PIS.v3),1)/size(RecTimeR_PIS.v3,1)>.1;
%PISisRec=(PISisRecG+PISisRecR)==2;
PISisRecG2=abs(nanmean(RecTimeG_PIS.v3,1)-EEG_PIS_times2)<.75;
% PISisRecR2=abs(nanmean(RecTimeG_PIS.v3,1)-EEG_PIS_times2)<.75;
PISisRec=(PISisRecG+PISisRecG2)==2;
% PISisRec=(PISisRecG+PISisRecR+PISisRecG2+PISisRecR2)==4;
RecTimeG_PIS.true=RecTimeG_PIS.v3(:,PISisRec);
%RecTimeR_PIS.true=RecTimeR_PIS.v3(:,PISisRec);
EEG_PIS_times3=EEG_PIS_times2(PISisRec);

%% RecTime PLOTS
if plotON==1
tempG=Fc1Gdff_flt2; %tempR=Fc1Rdff_flt2;

%Mean PIS, Sz and CSD traces and times 
figure
suptitle('Means')
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
plot(F1_ts,Fc1Gdff_flt1_Mean,'color', [0 1 0]);
% plot(F1_ts,Fc1Rdff_flt1_Mean-1,'color',[1 0 1]);
line([MeanSzRecTimeG;MeanSzRecTimeG],repmat([-4;5],1,2),'color',[0 0.5 0],'LineStyle','--')
% line([MeanSzRecTimeR;MeanSzRecTimeR],repmat([-4;5],1,2),'color',[0.5 0 0.5],'LineStyle','--')
line([MeanCSDRecTimeG;MeanCSDRecTimeG],repmat([-4;5],1,2),'color',[0 0.5 0],'LineStyle','-.')
% line([MeanCSDRecTimeR;MeanCSDRecTimeR],repmat([-4;5],1,2),'color',[0.5 0 0.5],'LineStyle','-.')
line([EEG_PIS_times2;EEG_PIS_times2],repmat([-4;5],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)
if death==1
    line([deathTime;deathTime],[-4;5],'color',[.5 0 0],'linewidth',2)
end


% Individual cell sz recruitment times

figure
suptitle('Sz')
subplot(1,2,1)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt2,1)
    if isRecruitedG(ii)==1 %incldG(ii)==1
        plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 0]);
    else
        plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 1]);
    end
    line([RecTimeG_sz2(ii);RecTimeG_sz2(ii)],repmat([ii*10;ii*10+1],1,2),'color',[0 0 1],'LineWidth',2)
    line([MeanSzRecTimeG;MeanSzRecTimeG],repmat([0;(size(Fc1Gdff_flt2,1)+1)*10],1,2),'color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',.5)
end
% 
% subplot(1,2,2)
% plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
% hold on
% for jj=1:size(Fc1Rdff_flt2,1)
%     if isRecruitedR(jj)==1
%         plot(F1_ts,tempR(jj,:)+jj*6,'color',[1 0.25 1]);
%     else
%         plot(F1_ts,tempR(jj,:)+jj*6,'color',[1 0.75 1]);
%     end
%     line([RecTimeR_sz2(jj);RecTimeR_sz2(jj)],repmat([jj*6;jj*6+1],1,2),'color',[0 0 1],'LineWidth',2)
%     line([MeanSzRecTimeR;MeanSzRecTimeR],repmat([0;(size(Fc1Rdff_flt2,1)+1)*10],1,2),'color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',.5)
% end
%ylim([-6.5 5])
%xlim([6 7])
hold off
clear ii jj

%% Plot CSD indiv rec times

% Individual cell csd/death recruitment times
figure
suptitle('CSD')
subplot(1,2,1)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt2,1)
    if isRecruitedG_csd(ii)==1 %incldG(ii)==1
        plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 0]);
    else
        plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 1]);
    end
    line([RecTimeG_csd(ii);RecTimeG_csd(ii)],repmat([ii*10;ii*10+1],1,2),'color',[0 0 1],'LineWidth',2)
    line([MeanCSDRecTimeG;MeanCSDRecTimeG],repmat([0;(size(Fc1Gdff_flt2,1)+1)*10],1,2),'color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',.5)
end

% subplot(1,2,2)
% plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
% hold on
% for jj=1:size(Fc1Rdff_flt2,1)
%     if isRecruitedR_csd(jj)==1
%         plot(F1_ts,tempR(jj,:)+jj*10,'color',[1 0.25 1]);
%     else
%         plot(F1_ts,tempR(jj,:)+jj*10,'color',[1 0.75 1]);
%     end
%     line([RecTimeR_csd(jj);RecTimeR_csd(jj)],repmat([jj*10;jj*10+1],1,2),'color',[0 0 1],'LineWidth',2)
%     line([MeanCSDRecTimeR;MeanCSDRecTimeR],repmat([0;(size(Fc1Rdff_flt2,1)+1)*10],1,2),'color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',.5)
% end
%ylim([-6.5 5])
%xlim([6 7])
hold off
clear ii jj


%% pre ictal spikes

figure
suptitle('PIS')
subplot(1,2,1)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt2,1)
    plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 0]);
    line([RecTimeG_PIS.v1(ii,:);RecTimeG_PIS.v1(ii,:)],repmat([ii*10+1;ii*10+1.75],1,size(RecTimeG_PIS.v1,2)),'color',[0 0 0],'LineWidth',1)
    line([RecTimeG_PIS.v3(ii,:);RecTimeG_PIS.v3(ii,:)],repmat([ii*10+3;ii*10+3.75],1,size(RecTimeG_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)

end
line([EEG_PIS_times2;EEG_PIS_times2],repmat([0;(size(Fc1Gdff_flt2,1)+1)*10],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)
% 
% subplot(1,2,2)
% plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
% hold on
% for ii=1:size(Fc1Rdff_flt2,1)
%     plot(F1_ts,tempR(ii,:)+ii*10,'color', [1 0 1]);
%     line([RecTimeR_PIS.v1(ii,:);RecTimeR_PIS.v1(ii,:)],repmat([ii*10+1;ii*10+1.75],1,size(RecTimeR_PIS.v1,2)),'color',[0 0 0],'LineWidth',1)
%     line([RecTimeR_PIS.v3(ii,:);RecTimeR_PIS.v3(ii,:)],repmat([ii*10+3;ii*10+3.75],1,size(RecTimeR_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)
% 
% end
% line([EEG_PIS_times2;EEG_PIS_times2],repmat([0;(size(Fc1Rdff_flt2,1)+1)*10],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)

clear ii tempG tempR

end

%% Generate Data Structure
%General Variables
StructFieldNames={'szID'
    %Time/Sampling
    'death'
    'deathTime'
    'deathTime_2p'
    'deathTime_eeg'
    'firstSpkT'
    'firstSpkI'
    'fs_2p'
    'fs_eeg'
    'EEGbias'
    'EEG_ts'
    'F1_ts'
    %EEG
    'EEG'
    'PSD_F'
    'PSD_P'
    %F1 processing
    'F1G'
    'Fneu1G'
    'stat_G'
    'xCoorG'
    'yCoorG'
    'xCoorIG'
    'yCoorIG'
    'iscell_G'
    'iscellIG'
    'Fb1G'
    'Fneub1G'
    'Fc1G'
    'Fb1Gdff'
    'Fneub1Gdff'
    'Fc1Gdff'
    'Fneub1Gdff_flt1'
    'Fc1Gdff_flt1'
    'Fneub1Gdff_flt1'
    'Fc1Gdff_flt2'
%     'F1R'
%     'Fneu1R'
%     'stat_R'
%     'xCoorR'
%     'yCoorR'
%     'xCoorIR'
%     'yCoorIR'
%     'iscell_R'
%     'iscellIR'
%     'Fb1R'
%     'Fneub1R'
%     'Fc1R'
%     'Fb1Rdff'
%     'Fneub1Rdff'
%     'Fc1Rdff'
%     'Fneub1Rdff_flt1'
%     'Fc1Rdff_flt1'
%     'Fneub1Rdff_flt1'
%     'Fc1Rdff_flt2'
    %Features
    'EEG_spkTimes_prelim'
    'EEG_PIS_times2'
    'EEG_PIS_times3'
    'MeanCaG_PIS_times'
%     'MeanCaR_PIS_times'
    'PISisRec'
    'RecTimeG_PIS'
%     'RecTimeR_PIS'
    'MeanSzRecTimeG'
%     'MeanSzRecTimeR'
    'MeanCSDRecTimeG'
%     'MeanCSDRecTimeR'
    'RecTimeG_sz2'
%     'RecTimeR_sz2'
    'RecTimeG_csd'
%     'RecTimeR_csd'
    'isRecruitedG'
%     'isRecruitedR'
    'isRecruitedG_csd'
%     'isRecruitedR_csd'
    };    
for kk=1:length(StructFieldNames)
    eval(strcat('szStruct.',StructFieldNames{kk},'=',StructFieldNames{kk},';'));
end

%% SAVE Structure Varaible as the SeizureID (szID)
assignin('base', szID, szStruct);
save(strcat(szID,'.mat'),szID,'-v7.3')
%modified to save in the compressed form of .mat to save disc space





