% Determine PREICTAL Recruitment Times - Individual traces

%this script take a matrix of calcium traces (cells by row, time by column)
%and a vector of seeds time for an event (e.g. pre-ictal spikes) and finds the 'elbow' of
%the trace closest to that event time for each trace. The elbow is the
%beginning of the block of time where the integral of first deriviate of
%the trace is largest (selecting for the longest sustatined and steepist slope).
%
%This script can calc the recTime by first (max slope) or second derivative
%(max concavity) use Dflag varaible to indicate which derivative to use
%('1' or '2').

function RecT=IndivRecTimes3pre(CaTraceMatrix,SeedTimeVec,fs_2p,topP,Dflag,StDevX)

%Missing Variables 
if ~exist('Dflag')
    Dflag = 1;
end

%define trace matrix
Fc1Gdff_flt2=CaTraceMatrix;
%weight trace by gaussian around seizure start time
% trace1=1:size(CaTraceMatrix,2);
% sigma1=10;%in seconds
% gauss1=exp(-(trace1-SeedTime*fs_2p).^2/(sigma1*fs_2p)^2);
% %Fc1Gdff_flt2=CaTraceMatrix.*gauss1;


%calculate the derivates of the filtered traces
Fc1Gdff_flt2_df=diff(Fc1Gdff_flt2,1,2);
Fc1Gdff_flt2_d2f=diff(Fc1Gdff_flt2,2,2);

% determine all the points that are positive
temp1G=Fc1Gdff_flt2_df>=0;

% find block bounds
temp2G=diff(temp1G,1,2);

temp3G=cell(size(temp2G,1),2);
for ii=1:size(temp3G,1)%green
    temp3G{ii,1}=find(temp2G(ii,:)==1);
    temp3G{ii,2}=find(temp2G(ii,:)==-1);
end
clear ii

%integrate over blocks
temp4G=cell(size(temp2G,1),1);
maxDFwin=cell(size(temp2G,1),1);
for ii=1:length(temp4G) %GREEN
    if temp3G{ii,1}(1)>temp3G{ii,2}(1)%correct for first slope value index being negative
        temp3G{ii,2}=temp3G{ii,2}(2:end);
    end
    
    if temp3G{ii,1}(end)>temp3G{ii,2}(end)%correct for extra final slope value being positive
        temp3G{ii,1}=temp3G{ii,1}(1:end-1);
    end
    
    for jj=1:length(temp3G{ii,1}) %find the area under the first derivative curve of each postitive region
        temp4G{ii}(jj)=trapz(Fc1Gdff_flt2_df(ii,[temp3G{ii,1}(jj):temp3G{ii,2}(jj)]));
        %Also find the max slope in this region
        [maxDFwin{ii}(jj,1),maxDFwin{ii}(jj,2)]=max(Fc1Gdff_flt2_df(ii,[temp3G{ii,1}(jj):temp3G{ii,2}(jj)]));
        %Also find the max second derivative in this region
        [maxD2Fwin{ii}(jj,1),maxD2Fwin{ii}(jj,2)]=max(Fc1Gdff_flt2_d2f(ii,[temp3G{ii,1}(jj):temp3G{ii,2}(jj)]));
    end
end
clear ii jj

% find maximum block and pull boundariy indices of that block
% find each largest value in each cell array and then use this to index into temp3 to find the index (time point) of the seizure in the origional trace
temp5G=cell(size(temp2G,1),1);%gives the block in each cell that is seizure
RecTimes=nan(size(temp2G,1),length(SeedTimeVec));

for ii=1:size(temp2G,1)
    % sorting to find top % of values
    [feats_sorted, feats_sorted_inds] = sort(temp4G{ii},2,'descend');
    topN=ceil(length(temp4G{ii})*topP/100);
    feats_sorted = feats_sorted(1:topN);
    feats_sorted_inds = feats_sorted_inds(1:topN);
    
    %recalculate block index times based upon max slope time
    if Dflag==2 %choose to use max 2nd derivative for rec time
        temp5G{ii,1}=temp3G{ii,1}+maxD2Fwin{ii}(:,2)';
    else %default is max first derivative for rec time
        temp5G{ii,1}=temp3G{ii,1}+maxDFwin{ii}(:,2)';
    end
    
    for jj=1:length(SeedTimeVec)
        [time_min2, time_min_ind2] = min(abs((temp5G{ii,1}(feats_sorted_inds))/fs_2p - SeedTimeVec(jj)));
        if Dflag==2 %choose to use max 2nd derivative for rec time
            RecTimes(ii,jj) = (temp5G{ii,1}(feats_sorted_inds(time_min_ind2))+2)/fs_2p;
        else %default is max first derivative for rec time
            RecTimes(ii,jj) = (temp5G{ii,1}(feats_sorted_inds(time_min_ind2))+1)/fs_2p;
        end
    end
    
    clear jj
    
end

clear ii


%Filter out Traces not within window for calcium
EventWin=1.25; %seconds

if size(RecTimes,2)==size(SeedTimeVec,2)
    EventDiff=RecTimes-SeedTimeVec;
else
    EventDiff=RecTimes-SeedTimeVec';
end

RecTimes(abs(EventDiff)>EventWin)=NaN;
RecT.v1=RecTimes;

%Filter Out Traces Not Recruited

valGauss=nan(size(RecTimes));
peakPIS=nan(size(RecTimes));
for ii=1:size(RecTimes,1)%loop over cells
    for jj=1:size(RecTimes,2)%loop over spikes
        if ~isnan(RecTimes(ii,jj)) 
            if RecTimes(ii,jj)+EventWin<size(CaTraceMatrix,2)/fs_2p
                tempTrace=CaTraceMatrix(ii,RecTimes(ii,jj)*fs_2p:(RecTimes(ii,jj)+EventWin)*fs_2p);%window for integration and peak finding
                valGauss(ii,jj)=trapz(tempTrace-min(tempTrace));
                peakPIS(ii,jj)=max(tempTrace);
                %meanPIS(ii,jj)=mean(CaTraceMatrix(ii,(RecTimes(ii,jj)-EventWin)*fs_2p+1:(RecTimes(ii,jj)+EventWin)*fs_2p));
                clear tempTrace
            else
                tempTrace=CaTraceMatrix(ii,RecTimes(ii,jj)*fs_2p:end);
                valGauss(ii,jj)=trapz(tempTrace-min(tempTrace));
                peakPIS(ii,jj)=max(tempTrace);
                %meanPIS(ii,jj)=mean(CaTraceMatrix(ii,(RecTimes(ii,jj)-EventWin)*fs_2p+1:end);
                clear tempTrace
            end
        end
    end
end
clear ii jj


stdThreshTraceSort=sort(CaTraceMatrix(:,(SeedTimeVec(1)-10)*fs_2p:end),2);
%CaSTDthresh=std(CaTraceMatrix(:,(SeedTimeVec(1)-10)*fs_2p:end),[],2);
CaSTDthresh=std(stdThreshTraceSort(:,1:size(stdThreshTraceSort,2)/2),[],2);
CaGaussThresh=2.5*CaSTDthresh*fs_2p/2/0.3989;%gaussian threshhold with height=std of preictal period, sigma=1/2 second in frames (2s width at 2 sigma two tail)
CaPeakThresh=mean(CaTraceMatrix(:,(SeedTimeVec(1)-10)*fs_2p:end),2)+StDevX*CaSTDthresh;
%CaMaxThresh=0.25*max(CaTraceMatrix-min(CaTraceMatrix,[],2),[],2);
%CaPeakThresh=mean(stdThreshTraceSort(:,1:size(stdThreshTraceSort,2)/2),2)+StDevX*CaSTDthresh;

RecTimes2=RecTimes;
%RecTimes3=RecTimes;

RecTimes(valGauss<repmat(CaGaussThresh,1,size(RecTimes,2)))=NaN;
RecT.v2=RecTimes;

RecTimes2(peakPIS<repmat(CaPeakThresh,1,size(RecTimes2,2)))=NaN;%ATTEMPTING TO FIX DETECTION

RecT.v3=RecTimes2;

%RecTimes3(peakPIS<repmat(CaMaxThresh,1,size(RecTimes2,2)))=NaN;

%RecT.v4=RecTimes3;


