%MASTER PROCESSING SCRIPT for Mouse Seizure Event Detection in 2PCI data
%
%
%Input .mat file from Suite2P
%Output is .mat file structure with fields to work with data for anaysis

%Outline
%  0) Define Variables/Parameters
%  1) Open EEG from EDF and annotation files
%  2) Open F and Fneu values for cells
%  3) Sort the F values into true cells and generate sorting indicies
%  4) Subtract background signal and neuropil contamination
%  5) Determine dF/F
%  6) Generate population mean traces
%  7) Correct for EEG Bias (time shift) between 2p and eeg
%  8) Power Spectral Density Determinations
%  9) Death Determination
% 10) Mean Seziure and TSW Seed Time Determination
% 11) Mean Pre-ictal Spike Seed Time Determiantion ---EEG_PIS_times
% 12) EEG Kernal Generation and Refinement of EEG PIS detection ---EEG_PIS_times2
% 13) First PIS determination and Final PIS SEED adjustment ---EEG_PIS_times2
% 14) Filtering of Individual Traces to 1 Hz
% 15) Individual trace seizure recruitment determination (f') ---RecTimeG_sz2, isRecruitedG
% 16) Individual trace TSW recruitment determination (f') ---RecTimeG_csd, isRecruitedG_csd
% 17) Individual trace PIS recruitment determination (f') ---RecTimeG_PIS.v3
% 18) Define if spike included in analysis (>10% cell rec) ---EEG_PIS_times3, recTimeG_PIS.true
% 19) Save data as a matlab structure

% Non Matlab Native Functions required:
% BO_IIS_Detect_Power_Rolston_v5
% edfread
% hifi
% IndivRecTimes3
% IndivRecTimes3pre
% lofi
% MeanRecTimes3pre
% ParseEDF
% SternRollAvg
%
% Version 231116 Matthew A. Stern and Eric R. Cole, Emory University
% Contact: matthew.a.stern@emory.edu or matt@matthewastern.com
%
% If using this code please cite our work.

%% Parameters 1
%Enter your parameters for the recording below
szID='S77ptz2'; %Enter seizure recording identification name
PIS_On=1; %Turn on (=1) if your recording has pre-ictal spikes
Sz_On=1; %Turn on (=1) if your recording has a seizure
TSW_On=1; %Turn on (=1) if your recording has a terminal spreading wave (e.g. spreading depolarization)
plotON=0; %Turn on (=1) if you want summary plots to be printed to screen
fs_eeg=2000; %Enter EEG sampling frequency
fs_2p=30; %Enter 2-photon calcium imaging sampling frequency
EEGbias=0;%if EEG seems shifted relative to the Ca signal(+ is EEG forward/ahead relative to calcium); default is 0
EEGLoader_On=1; %turn on (=1) if you want to use the EDF loader included
EEG=nan; %Enter your EEG data as a row vector if you choose to not use the EEG loader included

%% Loader for EEG
%This script is written to take an EDF file and load the first channel of 
%the recording and then parse it into segments based upon TTL pulse 
%annotation from the included txt file output from a pinnacle technologies
%recording.

%Make sure your EEG and annotation files are in the directory for analysis 
%and named export1.edf and annotations1.txt

if EEGLoader_On==1 %takes the second EEG segment in recording
    [EEG_all]=ParseEDF('export1.edf','annotations1.txt',fs_eeg,plotON);
    EEG=EEG_all.part{2};
end

%Write a timestamp for plotting
EEG_ts=[1/fs_eeg:1/fs_eeg:length(EEG)/fs_eeg];
   

%% Loader of Calcium Data (Single Channel)
%This section loads the suite2p output

load('suite2p/plane0/Fall.mat') %Enter data path for the suite2p .mat file
F_G=F; Fneu_G=Fneu; iscell_G=iscell; ops_G=ops; spks_G=spks; stat_G=stat;
clear F Fneu iscell ops redcell spks stat


%% Select Only Selected ROIs
%This section sorts the suite2P output to only included cells that were
%defined as iscell in suite2p. Note that it also translates the suite2p
%cell location to cartesian coordinates for easier matlab plotting
% 
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
        yCoorG(jj)=513-stat_G{ii}.med(1);%Modified to flip the y axis to match the images without needing to change axis orientation 
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

%Write a timestamp for plotting
F1_ts=[1/fs_2p:1/fs_2p:length(F1G)/fs_2p];


%% Baseline shift and neuropil subtraction
%This section performs a baseline shift in the data to remove background
%signal, defined as the global neuropil minima fluorescence. This ensures 
%transients are roughly in the same positive range. This section also subtracts a
%percentage of the local neuropil signal to ensure the somatic signal is 
%not contaminated by out of plane fluorescence, given the typical 5mm point 
%spread function in 2P imaging.    
%Fb: indicates background subtracted; Fc neuropil corrected

Fb1G=F1G-min(min(Fneu1G));%shifts F trace by the global Fneu minimum (defined as background)
Fneub1G=Fneu1G-min(min(Fneu1G));% must shift the neuropil too so it can be properly scaled
Fc1G=Fb1G-(0.7.*Fneub1G);%gives neuropil subtracted data (serves as pure somatic F)


%% Generate Mean Population Traces
%This section generates population mean traces to be used in later sections

F1G_Mean=mean(F1G,1);%population mean of raw F
Fneu1G_Mean=mean(Fneu1G,1);%population mean of raw F
Fb1G_Mean=mean(Fb1G,1);%cell population mean of background subtracted data 
Fneub1G_Mean=mean(Fneub1G,1);%neuropil population mean of background subtracted data 
Fc1G_Mean=mean(Fc1G,1);%population mean of background and neuropil subtracted data 


%% Generate dF/F
%This section generates the dF/F for each cell trace and the means of these
%(F-F0)/F0 with F0 being the mean of the beginning of the data 

F0sec=30; %Enter time period for F0 in seconds

%neuropil corrected somatic signal (Fc)
Fc1Gdff=(Fc1G-mean(Fc1G(:,1:round(F0sec*fs_2p)),2))./mean(Fc1G(:,1:round(F0sec*fs_2p)),2);%dF/F
Fc1Gdff_Mean=mean(Fc1Gdff,1);%population mean of dF/F
Fc1Gdff_std=std(Fc1Gdff);%standard deviation of the population mean of dF/F

%background subtracted soma (Fb)
Fb1Gdff=(Fb1G-mean(Fb1G(:,1:round(F0sec*fs_2p)),2))./mean(Fb1G(:,1:round(F0sec*fs_2p)),2);
Fb1Gdff_Mean=mean(Fb1Gdff,1);
Fb1Gdff_std=std(Fb1Gdff);

%background subtracted neuropil (Fneu)
Fneub1Gdff=(Fneub1G-mean(Fneub1G(:,1:round(F0sec*fs_2p)),2))./mean(Fneub1G(:,1:round(F0sec*fs_2p)),2);
Fneub1Gdff_Mean=mean(Fneub1Gdff,1);
Fneub1Gdff_std=std(Fneub1Gdff);


%% Filter (low pass) version 1: on dffs calculated from raw unfiltered data
%This section performs a lowpass filter to smooth the data; here primarily
%for plotting in the accompanying script for interim analysis plotting

lpfilt=5;
%hpfilt=.1;

Fb1Gdff_flt1=nan(size(Fb1Gdff));
Fneub1Gdff_flt1=nan(size(Fneub1Gdff));
Fc1Gdff_flt1=nan(size(Fc1Gdff));
%Fc1Gdff_flt1hp=nan(size(Fc1Gdff));

for ii=1:size(F1G,1)
    Fb1Gdff_flt1(ii,:)=lofi(Fb1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fneub1Gdff_flt1(ii,:)=lofi(Fneub1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    Fc1Gdff_flt1(ii,:)=lofi(Fc1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    %Fc1Gdff_flt1hp(ii,:)=hifi(Fc1Gdff_flt1(ii,:),10^6/F1_fs,hpfilt,'verbose',0);
end

clear lpfilt hpfilt
clear ii jj

% Generate Means for filtered data version 1
Fb1Gdff_flt1_Mean=mean(Fb1Gdff_flt1,1);
Fneub1Gdff_flt1_Mean=mean(Fneub1Gdff_flt1,1);
Fc1Gdff_flt1_Mean=mean(Fc1Gdff_flt1,1);

%Fc1Gdff_flt1hp_Mean=mean(Fc1Gdff_flt1hp,1);


%% Correction of EEG bias
%This sections corrects for a shift in the time syncing in EEG relative 
%to the Ca data.

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
%This section performs a PSD anaysis to find spectral features in the 
%frequency domain in the EEG data 

PSDbin=1024; %freq resolution
[~,PSD_F]=pwelch(EEG,[],[],PSDbin,fs_eeg); %PSD
PSDwin=1;%in seconds; sliding window

EEGstart=[1:PSDwin*fs_eeg/2:length(EEG)];%sliding window width 2x PSDwin shifting by 1/2 PSDwin
PSD_P=nan(length(PSD_F),length(EEGstart));
for ii = 1:numel(EEGstart)-2
    EEG1=EEG(EEGstart(ii):EEGstart(ii)+PSDwin*fs_eeg-1);
    [PSD_P(:,ii),~]=pwelch(EEG1,[],[],PSDbin,fs_eeg);
    clear EEG1
end

PSD_totP=sum(PSD_P,1);
PSD_theta=sum(PSD_P(and(PSD_F>3,PSD_F<15),:),1); %extended theta range (high delta to low alpha)
PSD_gammaL=sum(PSD_P(and(PSD_F>20,PSD_F<55),:),1); %extended gamma range (beta into gamma)
PSD_100plus=sum(PSD_P(PSD_F>=100 & PSD_F<250 ,:),1); %high frequency power
PSD_100less=sum(PSD_P(PSD_F<100,:),1); %low frequency power
Ratio100=PSD_100plus./PSD_100less; %high to low power ratio
RA_Ratio100=SternRollAvg(Ratio100,60/PSDwin); % +/- 30 second window roll avg (time step here is 0.5s)


%% Deriving recording Features (Death)
%Section is used to determine if the animal died post seizure
%We observed that in the terminal seizure instance there was a period of
%high frequency activity following the seizure, preceding the terminal
%spreading depolarization during death.
%Note that occasionaly this may also be observed following seiures without
%terminal spreading waves and should be manually corrected for

if max(RA_Ratio100) > 1
    death=1;
    deathTime=find(Ratio100>1,1)*PSDwin/2;
    deathTime_eeg=deathTime*fs_eeg;
    deathTime_2p=deathTime*fs_2p;
else
    death=0;
    [deathTime,deathTime_2p,deathTime_eeg]=deal(NaN); 
end


%% Mean Trace Seizure and CSD Time Seeds (half max integral method)
%This section is used to determine the mean seed times for the seizure and
%TSW using the neuropil signal.

if Sz_On==1

%Mean Trace low pass filter
lpfilt=1;%low pass fitlr in Hz
FGdff_MeanSm=lofi(Fneub1Gdff_Mean,10^6/fs_2p,lpfilt,'verbose',0); %smoothed trace
clear lpfilt

% determines all the points that are greater than the average of the min and
% max values of the trace (above half max)

if death==1 %to account for an animal dying the data 1 min beyond death is excluded (1 min is to account for hypoxic wave to propogate uncertatinty in exact time of death)
    postDeathT=deathTime_2p+fs_2p*60;
    temp1G=FGdff_MeanSm(1:postDeathT)>=mean([min(FGdff_MeanSm(1:postDeathT)),max(FGdff_MeanSm(1:postDeathT))]);
else
    temp1G=FGdff_MeanSm>=mean([min(FGdff_MeanSm),max(FGdff_MeanSm)]);
end

% find block bounds of above half max region
temp2G=diff(temp1G,1,2);
temp3G=cell(1,2);
temp3G{1}=find(temp2G==1);
temp3G{2}=find(temp2G==-1);

%integrate over blocks
temp4G=cell(1,1);

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

% find maximum block and pull boundary indices of that block
% find each largest value in each cell array and then use this to index into temp3 to find the index (time point) of the seizure in the origional trace
temp5G=cell(1,2);%gives the block in each cell that is seizure

[temp5G{1},temp5G{2}] = sort(temp4G{1},2,'descend');%sort blocks by size (block index is cell 2)
%determine the block bounds of the seizure and CSD blocks
if TSW_On==1
    MeanSzRecTimeG = min(temp3G{1}(temp5G{2}(1:2))+1)/fs_2p;%defines the Sz block as the first block of the top 2 blocks; get the first bound of that block which is the max recruitment time (add 1 due to use of diff fn)
    MeanCSDRecTimeG = max(temp3G{1}(temp5G{2}(1:2))+1)/fs_2p; %defines the CSD block as the last block of the top 2 blocks
else
    %determine the block bounds of the seizure blocks
    MeanSzRecTimeG = (temp3G{1}(temp5G{2}(1))+1)/fs_2p; %defines the CSD block as the last block of the top 2 blocks
    MeanCSDRecTimeG = NaN;
end

%Note that sometimes the differential between the seizure and TSW max F can
%cause issues in determining the half max appropriate to each event so
%these are recalibrated to the seizure and this process is run again below 
%to ensure correct seizure recruitment time determination
newMaxG = max(FGdff_MeanSm((MeanSzRecTimeG-2)*fs_2p:(MeanSzRecTimeG+5)*fs_2p));

clear temp*

% Iteration 2 (with new max value threshold)

if death==1 %to account for an animal dying the data 1 min beyond death is excluded (1 min is to account for hypoxic wave to propogate uncertatinty in exact time of death)
    postDeathT=round(deathTime_2p+fs_2p*60);
    temp1G=FGdff_MeanSm(1:postDeathT)>=mean([min(FGdff_MeanSm(1:postDeathT)),newMaxG]);
else
    temp1G=FGdff_MeanSm>=mean([min(FGdff_MeanSm),newMaxG]);   
end

% find block bounds
temp2G=diff(temp1G,1,2);
temp3G=cell(1,2);
temp3G{1}=find(temp2G==1);
temp3G{2}=find(temp2G==-1);

%integrate over blocks
temp4G=cell(1,1);
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

% find maximum block and pull boundary indices of that block
% find each largest value in each cell array and then use this to index into temp3 to find the index (time point) of the seizure in the origional trace
temp5G=cell(1,2);%gives the block in each cell that is seizure
[temp5G{1},temp5G{2}] = sort(temp4G{1},2,'descend');%sort blocks by size (block index is cell 2)
%determine the block bounds of the seizure blocks
if TSW_On==1
    MeanSzRecTimeG = min(temp3G{1}(temp5G{2}(1:2))+1)/fs_2p;%defines the Sz block as the first block of the top 2 blocks; get the first bound of that block which is the max recruitment time (add 1 due to use of diff fn)
else
    MeanSzRecTimeG = (temp3G{1}(temp5G{2}(1))+1)/fs_2p; %defines the CSD block as the last block of the top 2 blocks
end

else
MeanSzRecTimeG=NaN;
MeanCSDRecTimeG=NaN;
end


%% Preictal Mean RecTime Generation
%This section determines the preictal spike recruitment seed times using
%a function that takes the array of cells traces and finds all
%the points of the top percentage (TopP: 15%) of their slope integral
%features. This is then cross refrenced against spikes detected in EEG. 

if PIS_On==1 %loop for next few code blocks for PIS

%runs the pre-ictal mean recrutiment function on the 1 hz smoothed mean neuropil data
PreSpkRecTimeG=MeanRecTimes3pre(FGdff_MeanSm,fs_2p,15);

%calculates an integral of the time 1 second following recruitment to filter
%out false positives
valPIS1sG=nan(length(PreSpkRecTimeG{1,1}));
for ii=1:length(PreSpkRecTimeG{1,1})
    tempTrace=FGdff_MeanSm(round(PreSpkRecTimeG{1,1}(ii)*fs_2p):(round(PreSpkRecTimeG{1,1}(ii)+1)*fs_2p));
    valPIS1sG(ii)=trapz(tempTrace-min(tempTrace));
    clear tempTrace
end
clear ii

%Filter the Ca PIS candidate by a theoretical baseline spike threshold 
%for the signal integral
CaisLrgG=valPIS1sG>std(FGdff_MeanSm(1:round(MeanSzRecTimeG*fs_2p)))*fs_2p/3/0.3989; %gaussian, height=std of preictal period, sigma=1/3 second in frames (1s width at 1.5 sigma)
for ii=1:size(PreSpkRecTimeG,2)
     PreSpkRecTimeG{1,ii}=PreSpkRecTimeG{1,ii}(CaisLrgG);
end

%Determines the EEG PIS candidates
S_th = 4; %power threshold scaler for median for PIS detection
winW = .7; %seconds for spike detection window exclusion width

EEG_spk_prelim = BO_IIS_Detect_Power_Rolston_v5(EEG',S_th,fs_eeg,winW,0); %use v4 to not use theta gamma power ratio component
EEG_spkTimes_prelim_old=EEG_spk_prelim.time_Sp_final{1};%grab eeg spikes indexed by largest change in feature point
EEG_spkTimes_prelim=EEG_spk_prelim.time_Sp_minPeak{1};%grab all possible EEG spike times indexed by the min peak (spike) time point

%Finds concordant EEG and Ca PIS within a window
CaDetectWin=.2; %EEG and Ca PIS concordance window in seconds

%find Ca spikes in mean traces that coorepsond to EEG spikes
DiffPreEEG_G=abs(PreSpkRecTimeG{1,1}-(EEG_spkTimes_prelim'));%generates a Ca spk by EEG spk matrix of time diff
MeanPreIsRecG=min(DiffPreEEG_G,[],1)<CaDetectWin & PreSpkRecTimeG{1,1}<MeanSzRecTimeG;%logical indexing across Ca spk min to indicate if recruited
MeanCaG_PIS_times=PreSpkRecTimeG{1,1}(MeanPreIsRecG);

%refine EEG spike times based on Ca  and the seizure recruitment time
EEGisTrue=and(min(DiffPreEEG_G,[],2)<CaDetectWin, EEG_spkTimes_prelim'<MeanSzRecTimeG);%final EEG spikes
EEG_PIS_times=EEG_spkTimes_prelim(EEGisTrue);%pre-ictal spike times 


%% EEG spike wave kernel and PIS refinement
%This section determines a spike kernel and is used to find missing PIS
%that might not have had high calcium response (e.g. <20% cells)

%determines a kernel by averaging the EEG around each PIS determined above 
EEGkernWin=1;% kernel full width in seconds
EEGswd_KernSet=nan(numel(EEG_PIS_times),EEGkernWin*fs_eeg);
for ii=1:numel(EEG_PIS_times)
    EEGswd_KernSet(ii,:)=EEG(floor((EEG_PIS_times(ii)-EEGkernWin/2)*fs_eeg)+1:floor((EEG_PIS_times(ii)+EEGkernWin/2)*fs_eeg));
end
EEGswd_kern=mean(EEGswd_KernSet,1);

% Convolve the kernel with eeg and then finds the local maxima of greatest 
% convolution (>5 std above baseline and at least 1 second apart) 
EEGswdConv=conv(flip(EEGswd_kern), EEG);
[~,EEGswd_spikeTimes]=findpeaks(EEGswdConv(1:round(MeanSzRecTimeG*fs_eeg)),'MinPeakHeight',5*std(EEGswdConv(1:round(MeanSzRecTimeG*fs_eeg))),'MinPeakDistance',fs_eeg);
EEGswd_spikeTimes=EEGswd_spikeTimes/fs_eeg-0.5;

%combine the new detected spikes into the calcium contstrained spikes
EEGswdPISdiff=EEG_PIS_times-EEGswd_spikeTimes';
EEGswdPISmin=min(abs(EEGswdPISdiff),[],2);
EEG_PIS_new=EEGswd_spikeTimes(EEGswdPISmin>1);
EEG_PIS_times2=sort([EEG_PIS_times,EEG_PIS_new]);
EEG_PIS_times2=EEG_PIS_times2(EEG_PIS_times2<MeanSzRecTimeG); %make sure to grab sentinal spike but not seizure

%plotting of the SWD kernel set and average
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


%% Interspike interval (ISI) and first spike determination
%This section looks for an ISI rolling average of 10s to designate 
%consistent spiking but then looks for the first spike in this train that 
%has an ISI within 15 seconds as an upper limit on the allowable ISI
ISI_RAwin=5; %rolling average window of ISI to include
RAtolerance=10; %minimum ISI length required
IndivTolerance=15; %minimum ISI following first spike

ISI=diff(EEG_PIS_times2); %Interspike interval lengths

%determine the first spike interval in the list of PIS spikes
if numel(ISI)<ISI_RAwin
    firstSpkI=find(ISI<IndivTolerance,1);%if less than spikes RA window allow a slightly larger tolerance for ISI to grab spikes to maximize spikes collected
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

firstSpkT=EEG_PIS_times2(firstSpkI); %Time of the first ISI in the recording
EEG_PIS_times2=EEG_PIS_times2(firstSpkI:end); %new EEG PIS times starting at first spike

end


%% Filtering for Individual Trace Recruitment Time Determination
% This section performs filtering to smooth for recruitment detection

lpfilt=1; %low pass filter in hz
%hpfilt=0.1;

Fc1Gdff_flt2=zeros(size(Fc1Gdff));
%Fc1Gdff_flt2_hp=zeros(size(Fc1Gdff_flt2));
for ii=1:size(Fc1Gdff_flt1,1)
    Fc1Gdff_flt2(ii,:)=lofi(Fc1Gdff(ii,:),10^6/fs_2p,lpfilt,'verbose',0);
    %Fc1Gdff_flt2_hp(ii,:)=hifi(Fc1Gdff_flt2(ii,:),10^6/fs_2p,hpfilt);
end

clear lpfilt ii


%% Individual seizure recruitment detection (slope integral)
%This section performs seizure detection on the individual traces using the
%following function. It grabs the point of a cell of greatest slope integral
%feature weighted by a gaussian around the seizure seed time.
%
% RecTimes=IndivRecTimes3(CaTraceMatrix,SeedTime,fs_2p,sigma1)
%
%Inputs:
%CaTraceMatrix: matrix of calcium traces (cells by row, time by column)
%SeedTime: seed time for an event (e.g. TSW or Seizure)
%fs_2p: sampling frequency of 2p imaging
%sigma1: Standard deviations for width of gaussian weighting


if Sz_On==1
    
RecTimeG_sz3_struct=IndivRecTimes3(Fc1Gdff_flt2,MeanSzRecTimeG,fs_2p,1); %selecting a 1s sigma for gaussian to prevent sentinal spike from being grabbed and also if the average speed of the seizure is 200um/s then it should not take more than 1.5 seconds to cross the field so a 1s sigma gives us a 3 second window at at least 33% (1.5sigma) amplification of feature occurs
RecTimeG_sz2=RecTimeG_sz3_struct.time;

% Seizure Recruitment filtering
%This section determines if the cells are recruited at the times above
%Determined if each cell's trace has at least a 20% increase in signal 
%during seizure and the seizure detected time needs to be determined to be 
%within 3 seconds of the seed time (if average speed is 100um/s and field 
%is <300um all values should be within 3s of mean giving a 6s window)

tprd=10;%time period around seizure (s)
gthresh=0.2; %green threshold typical G threshold .3

preSzMeanDffG=mean(Fc1Gdff_flt2(:,round((MeanSzRecTimeG-tprd)*fs_2p):round(MeanSzRecTimeG*fs_2p)),2);
postSzMeanDffG=mean(Fc1Gdff_flt2(:,round(MeanSzRecTimeG*fs_2p):round((MeanSzRecTimeG+tprd)*fs_2p)),2);
SzMeanDffDiffG=postSzMeanDffG-preSzMeanDffG;
isRecruitedG=and(SzMeanDffDiffG>gthresh,abs(RecTimeG_sz2-MeanSzRecTimeG)<3); 

clear tprd gthresh rthresh

else
[RecTimeG_sz2, isRecruitedG, SzMeanDffDiffG] = deal(NaN);
end


%% Individual CSD / Death recruitment detection
%This section determines the recruitment times for the TSW using the same
%function as used for the seizure detection above.

if TSW_On==1

RecTimeG_csd_struct=IndivRecTimes3(Fc1Gdff_flt2,MeanCSDRecTimeG,fs_2p,5);%average speed of 40um/s gives means it takes 7s to cross 280um field of view, therefore a sigma of 5s gives us a 10s window for gaussian weight 
RecTimeG_csd=RecTimeG_csd_struct.time;

% Recruited cell or not (CSD)
%trace needs to have at least a 20% increase in signal during CSD relative 
%to before it and the CSD detected time needs to be determined to be within 
%10 seconds of the seed time

tprd=5;%time period for Ca difference (s)
gthresh=0.2; %percentage increase needed

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
clear tprd gthresh rthresh

else
[RecTimeG_csd, isRecruitedG_csd, CSDMeanDffDiffG] = deal(NaN);
end


%% Individual Cell PIS Recruitment Time Detection
%This section determines the recruitment times of the individual cells to the
%PIS, taking all PIS including sentintal spike, and using f'-steepest slope for
%time indexing and from the top 20% of cells using the slope integral
%feature for detection and a peak height threshold of atleast 4 standard
%deviations above baseline.
%
% IndivRecTimes3pre(CaTraceMatrix,SeedTimeVec,fs_2p,topP,Dflag,StDevX)
%
%CaTraceMatrix: matrix of calcium traces (cells by row, time by column)
%SeedTimeVec: vector of seeds time for an event (e.g. pre-ictal spikes)
%fs_2p: sampling frequency of 2p imaging
%topP: percentage of top cells to include in the selection of traces close to the seed times
%Dflag: Indicates which feature to index by (1:max slope or 2: max concavity (elbow))
%StDevX: Standard deviations for peak height threshold

if PIS_On==1
if isnan(MeanSzRecTimeG)
    RecTimeG_PIS=IndivRecTimes3pre(Fc1Gdff_flt2,EEG_PIS_times2,fs_2p,20,1,4);%fed the EEG PIS times
else
    RecTimeG_PIS=IndivRecTimes3pre(Fc1Gdff_flt2(:,1:MeanSzRecTimeG*fs_2p),EEG_PIS_times2,fs_2p,20,1,4);%fed the EEG PIS times
end

%Refine Spikes with min 10% recruitment of cells in each channel
%and remove spikes where the mean is clearly before or after the EEG
%event window of 0.75 seconds (window used for spike seperation)

PISisRecG=nansum(~isnan(RecTimeG_PIS.v3),1)/size(RecTimeG_PIS.v3,1)>.1;
PISisRecG2=abs(nanmean(RecTimeG_PIS.v3,1)-EEG_PIS_times2)<.75;
PISisRec=(PISisRecG+PISisRecG2)==2;
RecTimeG_PIS.true=RecTimeG_PIS.v3(:,PISisRec);
EEG_PIS_times3=EEG_PIS_times2(PISisRec);

else
    [MeanCaG_PIS_times,PISisRec,firstSpkT,firstSpkI,RecTimeG_PIS,EEG_PIS_times2,EEG_PIS_times3]=deal(nan);
end


%% RecTime PLOTS
if plotON==1
tempG=Fc1Gdff_flt2;

%Mean PIS, Sz and CSD traces and times 
figure
suptitle('Means')
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
plot(F1_ts,Fc1Gdff_flt1_Mean,'color', [0 1 0]);
line([MeanSzRecTimeG;MeanSzRecTimeG],repmat([-4;5],1,2),'color',[0 0.5 0],'LineStyle','--')
line([MeanCSDRecTimeG;MeanCSDRecTimeG],repmat([-4;5],1,2),'color',[0 0.5 0],'LineStyle','-.')
line([EEG_PIS_times2;EEG_PIS_times2],repmat([-4;5],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)
if death==1
    line([deathTime;deathTime],[-4;5],'color',[.5 0 0],'linewidth',2)
end


% Indivdual cell Seizure recruitment 
figure
suptitle('Sz')
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
hold off
clear ii 

% Individual cell TSW recruitment
figure
suptitle('CSD')
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
hold off
clear ii


% Individual cell PIS recruitment
figure
suptitle('PIS')
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt2,1)
    plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 0]);
    line([RecTimeG_PIS.v1(ii,:);RecTimeG_PIS.v1(ii,:)],repmat([ii*10+1;ii*10+1.75],1,size(RecTimeG_PIS.v1,2)),'color',[0 0 0],'LineWidth',1)
    line([RecTimeG_PIS.v3(ii,:);RecTimeG_PIS.v3(ii,:)],repmat([ii*10+3;ii*10+3.75],1,size(RecTimeG_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)
end
line([EEG_PIS_times2;EEG_PIS_times2],repmat([0;(size(Fc1Gdff_flt2,1)+1)*10],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)
clear ii tempG

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
    %Features
    'EEG_spkTimes_prelim'
    'EEG_PIS_times2'
    'EEG_PIS_times3'
    'MeanCaG_PIS_times'
    'PISisRec'
    'RecTimeG_PIS'
    'MeanSzRecTimeG'
    'MeanCSDRecTimeG'
    'RecTimeG_sz2'
    'RecTimeG_csd'
    'isRecruitedG'
    'isRecruitedG_csd'
    };    
for kk=1:length(StructFieldNames)
    eval(strcat('szStruct.',StructFieldNames{kk},'=',StructFieldNames{kk},';'));
end

%% SAVE Structure Varaible as the SeizureID (szID)
assignin('base', szID, szStruct);
save(strcat(szID,'.mat'),szID,'-v7.3')
%modified to save in the compressed form of .mat to save disc space





