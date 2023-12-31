%COMPANION SCRIPT FOR PLOTTING DDO/DIO PROCESSING STEPS
%This script is to be used alongside the processing while variables are
%loaded into memory. In its current form this will not work to look at 
%data from the saved structure output of this processing script. This is
%primarily used to troubleshoot.

%% Run to define cell number to use for below recordings
F1G_n=[1:size(F1G,1)];%number of green ROIs

ictalFront_start = 1/fs_2p;  %Enter manual start and end times if desiredfor ictal front segment (minutes)
ictalFront_end = F1_ts(end)/60;

%defines start and end times for eeg/2p sampling rates
ictalFront_start_eeg = ictalFront_start*fs_eeg*60;
ictalFront_start_2p = ictalFront_start*fs_2p*60;
ictalFront_end_eeg = ictalFront_end*fs_eeg*60;
ictalFront_end_2p = ictalFront_end*fs_2p*60;

t_2p_ictalFront = F1_ts(ictalFront_start_2p:ictalFront_end_2p)/60;
t_eeg_ictalFront = EEG_ts(ictalFront_start_eeg:ictalFront_end_eeg)/60;


%% Investigate traces: raw trace raster
figure
%Green Ca Raster
ax(1)= subplot(2,1,1);
imagesc(F1_ts,F1G_n,F1G(xCoorIG,:)) %note xCoorIG is used here to reorder cells by x coordiantes
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
xlabel('Times (s)')
ylabel('Neuron')
%xlim([0 60])


%% Investigate Raster of F-neu compared to background subtrated
figure
%Green Ca Raster
ax(1)= subplot(2,1,2);
imagesc(F1_ts,F1G_n,Fb1G(xCoorIG,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('Soma Calcium Transient Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')
%xlim([0 60])

%Green Ca Raster
ax(1)= subplot(2,1,1);
imagesc(F1_ts,F1G_n,Fc1G(xCoorIG,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('Bkgnd sub Calcium Transient Raster ','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')
%xlim([0 60])


%% plot filtered mean data version 1
figure
subplot(2,1,1)
plot(Fb1Gdff_Mean,'color',[0 1 0])
hold on
plot(Fneub1Gdff_Mean+2,'color',[0 0.6 0])
plot(Fc1Gdff_Mean+4,'color',[0 0.3 0])
title('dF/F Transients','fontsize',12)
ylim([-1 6])

subplot(2,1,2)
plot(Fb1Gdff_flt1_Mean,'color',[0 1 0])
hold on
plot(Fneub1Gdff_flt1_Mean+2,'color',[0 0.6 0])
plot(Fc1Gdff_flt1_Mean+4,'color',[0 0.3 0])
title('Filtered dF/F Transients','fontsize',12)
ylim([-1 6])


%% Investigate: Population mean traces by processing step
figure
subplot(2,1,1)
%plot(F1G_Mean,'color',[0 1 0])
hold on
plot(Fb1G_Mean,'color',[0 0.6 0])
plot(Fc1G_Mean,'color',[0 0.3 0])
%plot(Fneu1G_Mean,'color',[0 0 1])
plot(Fneub1G_Mean,'color',[0 0 0.5])
title('Mean Transients ','fontsize',12)

subplot(2,1,2)
plot(Fb1Gdff_Mean,'color',[0 1 0])
hold on
plot(Fneub1Gdff_Mean+1,'color',[0 0.6 0])
plot(Fc1Gdff_Mean+3,'color',[0 0.3 0])
title('dF/F Transients','fontsize',12)
ylim([-1 10])


%% Plot Rasters of Fb Fneub and dff

figure

%Ca Rasters
ax(1)= subplot(4,1,1);
imagesc(F1_ts,F1G_n,Fc1G(yCoorIG,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('Soma Calcium Transient Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')

ax(1)= subplot(4,1,3);
imagesc(F1_ts,F1G_n,Fneub1G(yCoorIG,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('Neuropil Calcium Transient Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')

ax(1)= subplot(4,1,2);
imagesc(F1_ts,F1G_n,Fc1Gdff(yCoorIG,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='dF/F';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('Soma dF/F Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')

ax(1)= subplot(4,1,4);
imagesc(F1_ts,F1G_n,Fneub1Gdff(yCoorIG,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='dF/F';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('Neuropil dF/F Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')


%% Plot Population mean traces by processing step
figure
subplot(4,1,1)
plot(F1_ts,Fb1G_Mean,'color',[0 0.6 0])
hold on
plot(F1_ts,Fc1G_Mean,'color',[0 0.3 0])
plot(F1_ts,Fneub1G_Mean-50,'color',[0 0 1])
title('Mean Transients ','fontsize',12)
xlabel('Times (s)')
ylabel('F (AU)')
hold off

subplot(4,1,2)
plot(F1_ts,Fb1Gdff_Mean,'color',[0 1 0])
title('F dF/F Transients','fontsize',12)
xlabel('Times (s)')
ylabel('dF/F')
%ylim([-1 2])

subplot(4,1,3)
plot(F1_ts,Fneub1Gdff_Mean,'color',[0 1 0])

title('Fneu dF/F Transients','fontsize',12)
xlabel('Times (s)')
ylabel('dF/F')
%ylim([-1 2])

subplot(4,1,4)
plot(F1_ts,Fc1Gdff_flt1_Mean,'color',[0 1 0])
title('Fc dF/F Transients','fontsize',12)
xlabel('Times (s)')
ylabel('dF/F')
%ylim([-1 2])


%% Plot Population Mean Traces with EEG

%ranges in seconds
fullRng=[10.5*60 15.5*60];

%Calcium trace selection
caTraceSelG=Fc1Gdff_flt1_Mean;

figure
subplot(2,2,1:2)
plot(EEG_ts,EEG,'color',[0 0 0])
%title('EEG')
%xlabel('time (s)')
%xlim(fullRng)
%ylabel('uV')
%ylim([-3000 2000])
pbaspect([5 1 1])

subplot(2,2,3:4)
plot([1/fs_2p:1/fs_2p:length(caTraceSelG)/fs_2p],caTraceSelG,'color',[0 1 0])
%title('Mean Fc by Population')
%xlabel('time (s)')
%xlim(fullRng)
%ylabel('dF/F')
%ylim([-0.4 .8])

hold off
pbaspect([3 1 1])

clear caTraceSel*
clear *Rng


%% Spectrogram of EEG

PSD_maxF=300;
PSD_ylim=find(PSD_F<PSD_maxF,1,'last');
figure
subplot(2,1,1)
pcolor(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_F(1:PSD_ylim),10*log10(PSD_P(1:PSD_ylim,1:numel(EEGstart))))
shading flat
colormap(jet)
%ylim([0.5,50])
colorbar
subplot(2,1,2)
plot(EEG_ts,EEG.*10+30000,'color',[0.7 0.7 0.7]);
%plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_totP,'color',[0 0 0])
hold on
%plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_theta,'color',[1 0 0])
%plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_gammaL*5,'color',[0 0 1])
plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_100plus,'color',[1 0 0])
plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_100less,'color',[0 0 1])
plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,Ratio100.*10000,'color',[0 1 0])
plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,RA_Ratio100.*10000,'color',[0 1 0.5])
if death>0
    line([deathTime;deathTime],[0;100000],'color',[.5 0 0],'linewidth',2)
end
ylim([0 10^5])   

%% Examine mean traces relative to EEG for ictal front period

figure
plot(t_eeg_ictalFront,EEG(ictalFront_start_eeg:ictalFront_end_eeg)./1000+3,'color',[0.7 0.7 0.7]);
hold on
plot(t_2p_ictalFront,Fc1Gdff_flt1_Mean(ictalFront_start_2p:ictalFront_end_2p),'color', [0 1 0]);
%title('Mean Fc Transients')
%xlabel('time (min)')
%ylabel('dF/F')
ylim([-4 5])
xlim([ictalFront_start ictalFront_end])
yticks([])
pbaspect([1 .5 1])
hold off

set(gca,'Visible','off')


%% examine individual traces of filtered data for ictal front period
figure
plot(t_eeg_ictalFront,EEG(ictalFront_start_eeg:ictalFront_end_eeg)./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt1,1)
    plot(t_2p_ictalFront,Fc1Gdff_flt1(ii,ictalFront_start_2p:ictalFront_end_2p)-4+ii*6,'color', [0.25 1 0.25]);
end
%ylim([-6.5 5])
%xlim([6 7])
hold off
clear ii jj


%% Investigate traces: filtered trace raster

figure
%Ca Raster
ax(1)= subplot(2,1,1);
imagesc(F1_ts/60,size(Fc1Gdff_flt1,1),Fc1Gdff_flt1(yCoorIG,:)) %note xCoorIG is used here to reorder cells by x coordiantes
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
%c1.Label.String='dF/F';
clear c1
caxis([-4 10])
pbaspect([3 1 1])
%title('Filtered Fc Transient Raster','fontsize',12)
%xlabel('Times (s)')
%ylabel('Neuron')
xlim([ictalFront_start ictalFront_end])


%% EEG spike wave kernel

figure
for ii=1:size(EEGswd_KernSet,1)
    subplot(ceil(sqrt(size(EEGswd_KernSet,1))),ceil(sqrt(size(EEGswd_KernSet,1))),ii)
    plot(EEGswd_KernSet(ii,:))
    title(ii)
end
clear ii

figure
plot(EEGswd_kern)


%% Method to do detection based upon relative theta to gamma power

figure
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
line([EEG_spkTimes_prelim;EEG_spkTimes_prelim],repmat([3;3.75],1,length(EEG_spkTimes_prelim)),'color',[0 1 0],'linewidth',1)
%line([EEG_spkTimes_prelim1;EEG_spkTimes_prelim1],repmat([2.75;3.5],1,length(EEG_spkTimes_prelim1)),'color',[1 0 0],'linewidth',1)
line([EEG_PIS_times2;EEG_PIS_times2],repmat([2;2.75],1,length(EEG_PIS_times2)),'color',[0 0 1],'linewidth',1)

%% plot mean calcium and EEG
figure
plot(t_eeg_ictalFront,EEG(ictalFront_start_eeg:ictalFront_end_eeg)./1000+3,'color',[0.7 0.7 0.7]);
hold on
%plot(t_2p_ictalFront,Fc1Gdff_flt1_Mean(ictalFront_start_2p:ictalFront_end_2p),'color', [0 1 0]);
plot(t_2p_ictalFront,FGdff_MeanSm(ictalFront_start_2p:ictalFront_end_2p),'color', [0 1 0]);
line([MeanSzRecTimeG/60;MeanSzRecTimeG/60],repmat([-4;5],1,2),'color',[0 0.5 0],'LineStyle','--')
line([MeanCSDRecTimeG/60;MeanCSDRecTimeG/60],repmat([-4;5],1,2),'color',[0 0.5 0],'LineStyle','-.')
%title('Mean Fc Transients')
%xlabel('time (min)')
%ylabel('dF/F')
ylim([-4 5])
xlim([ictalFront_start ictalFront_end])
yticks([])
pbaspect([1 .5 1])
hold off

set(gca,'Visible','off')


%% plot the preictal recruitment times
figure
subplot(2,1,1)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
%line([EEG_spkTimes_prelim;EEG_spkTimes_prelim],repmat([2.75;3.5],1,length(EEG_spkTimes_prelim)),'color',[0 1 1],'linewidth',1)
%line([EEG_spkTimes_prelim_old;EEG_spkTimes_prelim_old],repmat([2.75;3.5],1,length(EEG_spkTimes_prelim_old)),'color',[1 1 0.5],'linewidth',1)
plot(F1_ts,Fc1Gdff_flt1_Mean,'color', [0.8 1 0.8]);
plot(F1_ts,FGdff_MeanSm,'color', [0 1 0]);


for ii=1:length(PreSpkRecTimeG{1,1})
        subplot(2,1,1)
        line([PreSpkRecTimeG{1,1}(ii);PreSpkRecTimeG{1,1}(ii)],[-.25;.5],'color',[0 0 1],'LineWidth',1)
%         subplot(2,1,2)
%         scatter(ii,PreSpkRecTimeG{1,3}(ii),'MarkerEdgeColor', [0 0 0])
end
subplot(2,1,1)
    line([MeanCaG_PIS_times;MeanCaG_PIS_times],repmat([-.25;.5],1,length(MeanCaG_PIS_times)),'color',[1 0 0],'LineWidth',1)
    line([EEG_PIS_times;EEG_PIS_times],repmat([3.75;4.5],1,length(EEG_PIS_times)),'color',[1 0 1],'linewidth',1)
    line([EEG_PIS_times2;EEG_PIS_times2],repmat([2.75;3.5],1,length(EEG_PIS_times2)),'color',[1 0 0],'linewidth',1)
    line([EEGswd_spikeTimes;EEGswd_spikeTimes],repmat([4.75;5.5],1,length(EEGswd_spikeTimes)),'color',[0 0 1],'linewidth',1)
    
subplot(2,1,2)
plot([-999/fs_eeg:1/fs_eeg:0,EEG_ts,EEG_ts(end)+1/fs_eeg:1/fs_eeg:EEG_ts(end)+999/fs_eeg],EEGswdConv/10^7)
hold on
line([0,EEG_ts(end)],[5*std(EEGswdConv(1:MeanSzRecTimeG*fs_eeg)/10^7),5*std(EEGswdConv(1:MeanSzRecTimeG*fs_eeg)/10^7)],'color',[0 0 0])
line([EEGswd_spikeTimes;EEGswd_spikeTimes],repmat([1;2],1,length(EEGswd_spikeTimes)),'color',[1 0 0],'linewidth',.5)


%% Spectral power relative to Ca signal

PSD_maxF=200;
PSD_ylim=find(PSD_F<PSD_maxF,1,'last');
figure
subplot(2,1,1)
pcolor(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_F(1:PSD_ylim),10*log10(PSD_P(1:PSD_ylim,1:numel(EEGstart))))
shading flat
colormap(jet)
%ylim([0.5,50])
%colorbar

subplot(2,1,2)
plot(EEG_ts,EEG.*10+30000,'color',[0.7 0.7 0.7]);
%plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_totP,'color',[0 0 0])
hold on
%plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_theta,'color',[1 0 0])
%plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_gammaL*5,'color',[0 0 1])
plot(F1_ts,Fc1Gdff_flt1_Mean.*10000,'color', [0 1 0]);
% plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_100plus,'color',[1 0 0])
% plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_100less,'color',[0 0 1])
% plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,Ratio100.*10000,'color',[0 1 0])
line([MeanSzRecTimeG;MeanSzRecTimeG],[0;50000],'color',[0 .7 0],'LineWidth',.5)
if death==1
    line([deathTime;deathTime],[0;50000],'color',[.5 0 0],'linewidth',2)
end


%% to write figure with sz recruitment times inserted
tempG=Fc1Gdff_flt2;

figure
plot(t_eeg_ictalFront,EEG(ictalFront_start_eeg:ictalFront_end_eeg)./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt2,1)
    if isRecruitedG(ii)==1 %incldG(ii)==1
        plot(t_2p_ictalFront,tempG(ii,ictalFront_start_2p:ictalFront_end_2p)+ii*10,'color', [0 1 0]);
    else
        plot(t_2p_ictalFront,tempG(ii,ictalFront_start_2p:ictalFront_end_2p)+ii*10,'color', [0 1 1]);
    end
    line([RecTimeG_sz2(ii)/60;RecTimeG_sz2(ii)/60],repmat([ii*10;ii*10+1],1,2),'color',[0 0 1],'LineWidth',2)
    line([MeanSzRecTimeG/60;MeanSzRecTimeG/60],repmat([0;2000],1,2),'color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',.5)
end
hold off
clear ii tempG


%% pre ictal individual ploting
tempG=Fc1Gdff_flt2;

figure
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt2,1)
    plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 0]);
    line([RecTimeG_PIS.v1(ii,:);RecTimeG_PIS.v1(ii,:)],repmat([ii*10+1;ii*10+1.75],1,size(RecTimeG_PIS.v1,2)),'color',[0 0 0],'LineWidth',1)
    line([RecTimeG_PIS.v2(ii,:);RecTimeG_PIS.v2(ii,:)],repmat([ii*10+2;ii*10+2.75],1,size(RecTimeG_PIS.v2,2)),'color',[1 0 0],'LineWidth',1)
    line([RecTimeG_PIS.v3(ii,:);RecTimeG_PIS.v3(ii,:)],repmat([ii*10+3;ii*10+3.75],1,size(RecTimeG_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)

end
line([EEG_PIS_times2;EEG_PIS_times2],repmat([0;size(Fc1Gdff_flt2,1)*10],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)



