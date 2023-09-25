%COMPANION SCRIPT FOR PLOTTING DDO/DIO PROCESSING STEPS

F1G_n=[1:size(F1G,1)];%number of green ROIs
F1R_n=[1:size(F1R,1)];%number of red ROIs

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
title('Yellow','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')
%xlim([0 60])

%Red Ca Raster
ax(1)= subplot(2,1,2);
imagesc(F1_ts,F1R_n,F1R(xCoorIR,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-150 250])
%pbaspect([3 1 1])
title('Red','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')
%xlim([0 60])


%% Investigate Raster of F-neu compared to background subtrated
figure
%Green Ca Raster
ax(1)= subplot(2,2,2);
imagesc(F1_ts,F1G_n,Fb1G(xCoorIG,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('Non-vGAT (jYCaMP1s) Soma Calcium Transient Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')
%xlim([0 60])

%Red Ca Raster
ax(1)= subplot(2,2,4);
imagesc(F1_ts,F1R_n,Fb1R(xCoorIR,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-150 250])
%pbaspect([3 1 1])
title('vGAT (jRGECO1a) Soma Calcium Transient Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')
%xlim([0 60])

%Green Ca Raster
ax(1)= subplot(2,2,1);
imagesc(F1_ts,F1G_n,Fc1G(xCoorIG,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('Non-vGAT (jYCaMP1s) Bkgnd sub Calcium Transient Raster ','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')
%xlim([0 60])

%Red Ca Raster
ax(1)= subplot(2,2,3);
imagesc(F1_ts,F1R_n,Fc1R(xCoorIR,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-150 250])
%pbaspect([3 1 1])
title('vGAT (jRGECO1a) Bkgnd sub Calcium Transient Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')
%xlim([0 60])


%% plot filtered mean data version 1
figure
subplot(2,2,1)
plot(Fb1Gdff_Mean,'color',[0 1 0])
hold on
plot(Fneub1Gdff_Mean+2,'color',[0 0.6 0])
plot(Fc1Gdff_Mean+4,'color',[0 0.3 0])
title('Green dF/F Transients','fontsize',12)
ylim([-1 6])

subplot(2,2,2)
plot(Fb1Gdff_flt1_Mean,'color',[0 1 0])
hold on
plot(Fneub1Gdff_flt1_Mean+2,'color',[0 0.6 0])
plot(Fc1Gdff_flt1_Mean+4,'color',[0 0.3 0])
title('Filtered Green dF/F Transients','fontsize',12)
ylim([-1 6])

subplot(2,2,3)
plot(Fb1Rdff_Mean,'color',[1 0 0])
hold on
plot(Fneub1Rdff_Mean+2,'color',[0.6 0 0])
plot(Fc1Rdff_Mean+4,'color',[0.3 0 0])
title('Red dF/F Transients','fontsize',12)
ylim([-1 6])

subplot(2,2,4)
plot(Fb1Rdff_flt1_Mean,'color',[1 0 0 ])
hold on
plot(Fneub1Rdff_flt1_Mean+2,'color',[0.6 0 0])
plot(Fc1Rdff_flt1_Mean+4,'color',[0.3 0 0])
title('Filtered Red dF/F Transients','fontsize',12)
ylim([-1 6])


%% Investigate: Population mean traces by processing step
figure
subplot(3,1,1)
%plot(F1G_Mean,'color',[0 1 0])
hold on
plot(Fb1G_Mean,'color',[0 0.6 0])
plot(Fc1G_Mean,'color',[0 0.3 0])
%plot(Fneu1G_Mean,'color',[0 0 1])
plot(Fneub1G_Mean,'color',[0 0 0.5])
title('Non-vGAT (jYCaMP1s) Mean Transients ','fontsize',12)


subplot(3,1,2)
%plot(F1R_Mean,'color',[1 0 0])
hold on
plot(Fb1R_Mean,'color',[0.6 0 0])
plot(Fc1R_Mean,'color',[0.3 0 0])
%plot(Fneu1R_Mean,'color',[0 0 1])
plot(Fneub1R_Mean,'color',[0 0 0.5])
title('vGAT (jRGECO1a) Mean Transients','fontsize',12)

subplot(3,1,3)
plot(Fb1Gdff_Mean,'color',[0 1 0])
hold on
plot(Fb1Rdff_Mean,'color',[1 0 0])
plot(Fneub1Gdff_Mean+1,'color',[0 0.6 0])
plot(Fneub1Rdff_Mean+1,'color',[0.6 0 0])
plot(Fc1Gdff_Mean+3,'color',[0 0.3 0])
plot(Fc1Rdff_Mean+3,'color',[0.3 0 0])
title('dF/F Transients','fontsize',12)
ylim([-1 10])


% %% Investigate: Plot of dff (by mean)
% figure
% plot(1/fs_2p:1/fs_2p:size(Fc1Gdff_Mean,2)/fs_2p,Fc1Gdff_Mean,'color',[0 1 0],'linewidth',1)
% hold on
% plot(1/fs_2p:1/fs_2p:size(Fc1Rdff_Mean,2)/fs_2p,Fc1Rdff_Mean,'color',[1 0 1],'linewidth',1)
% plot(1/fs_2p:1/fs_2p:size(Fc1Gdff_Mean,2)/fs_2p,Fc1Gdff_Mean+Fc1Gdff_std,'color',[0.5 1 0.5],'linewidth',0.5)
% plot(1/fs_2p:1/fs_2p:size(Fc1Gdff_Mean,2)/fs_2p,Fc1Gdff_Mean-Fc1Gdff_std,'color',[0.5 1 0.5],'linewidth',0.5)
% plot(1/fs_2p:1/fs_2p:size(Fc1Rdff_Mean,2)/fs_2p,Fc1Rdff_Mean+Fc1Rdff_std,'color',[1 0.5 1],'linewidth',0.5)
% plot(1/fs_2p:1/fs_2p:size(Fc1Rdff_Mean,2)/fs_2p,Fc1Rdff_Mean-Fc1Rdff_std,'color',[1 0.5 1],'linewidth',0.5)
% %xlim([])
% %ylim([])
% pbaspect([1.5 1 1])
% title('Mean Calcium Transients')
% xlabel('Time (s)')
% ylabel('dF/F')
% legend('Non-VGAT (jYCaMP1s)','VGAT (jRGECO1a)')
% hold off
% 
% %set(gca, 'box', 'off','color',[0 0 0])%plot background black


%% Plot Rasters of Fb Fneub and dff

figure

%Green Ca Rasters
ax(1)= subplot(4,2,1);
imagesc(F1_ts,F1G_n,Fc1G(yCoorIG,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('Non-vGAT (jYCaMP1s) Soma Calcium Transient Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')

ax(1)= subplot(4,2,3);
imagesc(F1_ts,F1G_n,Fneub1G(yCoorIG,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('Non-vGAT (jYCaMP1s) Neuropil Calcium Transient Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')

ax(1)= subplot(4,2,2);
imagesc(F1_ts,F1G_n,Fc1Gdff(yCoorIG,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='dF/F';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('Non-vGAT (jYCaMP1s) Soma dF/F Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')

ax(1)= subplot(4,2,4);
imagesc(F1_ts,F1G_n,Fneub1Gdff(yCoorIG,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='dF/F';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('Non-vGAT (jYCaMP1s) Neuropil dF/F Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')


%Red Ca Rasters
ax(1)= subplot(4,2,5);
imagesc(F1_ts,F1R_n,Fc1R(yCoorIR,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('vGAT (jRGECO1a) Soma Calcium Transient Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')

ax(1)= subplot(4,2,7);
imagesc(F1_ts,F1R_n,Fneub1R(yCoorIR,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='F (AU)';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('vGAT (jRGECO1a) Neuropil Calcium Transient Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')

ax(1)= subplot(4,2,6);
imagesc(F1_ts,F1R_n,Fc1Rdff(yCoorIR,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='dF/F';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('vGAT (jRGECO1a) Soma dF/F Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')

ax(1)= subplot(4,2,8);
imagesc(F1_ts,F1R_n,Fneub1Rdff(yCoorIR,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
c1.Label.String='dF/F';
clear c1
%caxis([-250 250])
%pbaspect([3 1 1])
title('vGAT (jRGECO1a) Neuropil dF/F Raster','fontsize',12)
xlabel('Times (s)')
ylabel('Neuron')


%% Plot Population mean traces by processing step
figure
subplot(3,2,1)
plot(F1_ts,Fb1G_Mean,'color',[0 0.6 0])
hold on
plot(F1_ts,Fc1G_Mean,'color',[0 0.3 0])
plot(F1_ts,Fneub1G_Mean-50,'color',[0 0 1])
title('Non-vGAT (jYCaMP1s) Mean Transients ','fontsize',12)
xlabel('Times (s)')
ylabel('F (AU)')
hold off

subplot(3,2,3)
plot(F1_ts,Fb1R_Mean,'color',[0.6 0 0.6])
hold on
plot(F1_ts,Fc1R_Mean,'color',[0.3 0 0.3])
plot(F1_ts,Fneub1R_Mean,'color',[0 0 1])
title('vGAT (jRGECO1a) Mean Transients','fontsize',12)
xlabel('Times (s)')
ylabel('F (AU)')
hold off

subplot(3,2,2)
plot(F1_ts,Fb1Gdff_Mean,'color',[0 1 0])
hold on
plot(F1_ts,Fb1Rdff_Mean,'color',[1 0 1])
title('F dF/F Transients','fontsize',12)
xlabel('Times (s)')
ylabel('dF/F')
%ylim([-1 2])

subplot(3,2,4)
plot(F1_ts,Fneub1Gdff_Mean,'color',[0 1 0])
hold on
plot(F1_ts,Fneub1Rdff_Mean,'color',[1 0 1])
title('Fneu dF/F Transients','fontsize',12)
xlabel('Times (s)')
ylabel('dF/F')
%ylim([-1 2])

subplot(3,2,5:6)
plot(F1_ts,Fc1Gdff_flt1_Mean,'color',[0 1 0])
hold on
plot(F1_ts,Fc1Rdff_flt1_Mean,'color',[1 0 1])
title('Fc dF/F Transients','fontsize',12)
xlabel('Times (s)')
ylabel('dF/F')
ylim([-1 2])


%% Plot Population Mean Traces with EEG

%ranges in seconds
fullRng=[10.5*60 15.5*60];

%Calcium trace selection
caTraceSelG=Fc1Gdff_flt1_Mean;
caTraceSelR=Fc1Rdff_flt1_Mean;

figure
subplot(2,2,1:2)
plot(EEG_ts,EEG,'color',[0 0 0])
%title('LFP')
%xlabel('time (s)')
%xlim(fullRng)
%ylabel('uV')
%ylim([-3000 2000])
pbaspect([5 1 1])

subplot(2,2,3:4)
plot([1/fs_2p:1/fs_2p:length(caTraceSelG)/fs_2p],caTraceSelG,'color',[0 1 0])
hold on
plot([1/fs_2p:1/fs_2p:length(caTraceSelR)/fs_2p],caTraceSelR,'color',[1 0 1])
%title('Mean Fc by Population')
%xlabel('time (s)')
%xlim(fullRng)
%ylabel('dF/F')
%ylim([-0.4 .8])
%legend({'non-VGAT','VGAT'},'Location','northwest')
%legend('boxoff')
hold off
pbaspect([3 1 1])

clear caTraceSel*
clear *Rng



%% CODE FOR MORE ASTHETIC VERSIONS OF ABOVE PLOTS 
% %% Plot population Mean Traces with pre ictal phases
% 
% %ranges in seconds
% fullRng=[300 900];
% earlyRng=[475 500];
% lateRng=[720 745];
% 
% %Calcium trace selection
% caTraceSelG=Fc1Gdff_flt1_Mean;
% caTraceSelR=Fc1Rdff_flt1_Mean;
% 
% figure
% subplot(4,2,1:2)
% plot(EEG_ts,EEG,'color',[0 0 0])
% title('LFP')
% %xlabel('time (s)')
% xlim(fullRng)
% ylabel('uV')
% ylim([-2500 2500])
% 
% subplot(4,2,3:4)
% plot([1/fs_2p:1/fs_2p:length(caTraceSelG)/fs_2p],caTraceSelG,'color',[0 1 0])
% hold on
% plot([1/fs_2p:1/fs_2p:length(caTraceSelR)/fs_2p],caTraceSelR,'color',[1 0 1])
% title('Mean Fc by Population')
% xlabel('time (s)')
% xlim(fullRng)
% ylabel('dF/F')
% %ylim([-0.4 .8])
% legend({'non-VGAT','VGAT'},'Location','northwest')
% legend('boxoff')
% hold off
% 
% subplot(4,2,5)
% plot(EEG_ts,EEG,'color',[0 0 0])
% title('Early Spiking')
% %xlabel('time (s)')
% ylabel('uV')
% xlim(earlyRng)
% ylim([-2500 2500])
% set(gca, 'box', 'on')
% 
% subplot(4,2,7)
% plot([1/fs_2p:1/fs_2p:length(caTraceSelG)/fs_2p],caTraceSelG,'color',[0 1 0])
% hold on
% plot([1/fs_2p:1/fs_2p:length(caTraceSelR)/fs_2p],caTraceSelR,'color',[1 0 1])
% %title('Early Spiking')
% xlabel('time (s)')
% ylabel('dF/F')
% xlim(earlyRng)
% ylim([-.5 1.5])
% hold off
% 
% subplot(4,2,6)
% plot(EEG_ts,EEG,'color',[0 0 0])
% title('Late Spiking')
% %xlabel('time (s)')
% ylabel('uV')
% xlim(lateRng)
% ylim([-2500 2500])
% set(gca, 'box', 'on')
% 
% subplot(4,2,8)
% plot([1/fs_2p:1/fs_2p:length(caTraceSelG)/fs_2p],caTraceSelG,'color',[0 1 0])
% hold on
% plot([1/fs_2p:1/fs_2p:length(caTraceSelR)/fs_2p],caTraceSelR,'color',[1 0 1])
% %title('Late Spiking')
% xlabel('time (s)')
% ylabel('dF/F')
% xlim(lateRng)
% ylim([-0.5 1.5])
% hold off
% 
% clear caTraceSel*
% clear *Rng
% 
% 

% %% Plot Rasters of Fb Fneub and dff
% 
% xLimRange=[300 800];
% 
% figure
% 
% %Green Ca Rasters
% 
% ax(1)= subplot(2,1,1);
% imagesc(F1_ts,F1G_n,Fc1Gdff_flt1(yCoorIG,:))
% set(gca,'Ydir','normal')
% shading flat
% colormap(ax(1),jet)
% c1=colorbar;
% %c1.Label.String='dF/F';
% clear c1
% caxis([-4 10])
% %pbaspect([3 1 1])
% %title('Non-VGAT (jYCaMP1s) Somatic Caclium Transients','fontsize',12)
% %xlabel('Times (s)')
% %ylabel('Neuron')
% xlim(xLimRange)
% 
% ax(1)= subplot(2,1,2);
% imagesc(F1_ts,F1R_n,Fc1Rdff_flt1(yCoorIR,:))
% set(gca,'Ydir','normal')
% shading flat
% colormap(ax(1),jet)
% c1=colorbar;
% %c1.Label.String='dF/F';
% clear c1
% caxis([-4 10])
% %pbaspect([3 1 1])
% %title('VGAT (jRGECO1a) Somatic Calcium Transients','fontsize',12)
% %xlabel('Times (s)')
% %ylabel('Neuron')
% xlim(xLimRange)
% 
% clear  xLimRange














%% Manually Entered Values for Plotting

ictalFront_start = 1/fs_2p;     %s51ptz1 3.5; s69sz1:7.0; s71sz1 3.3; s75ptz1:12.5; s75ptz2:5.5; s77sz1e1:2.5; s77sz1e2:5.8; s77sz2:4.7; s77sz3:1.4 ;  %manual start and end times for ictal front segment (minutes)
ictalFront_end = F1_ts(end)/60;      %s51ptz1 5.0; s69sz1:8.67; s71sz1 4.8; s75ptz1:13.5; s75ptz2:6.5; s77sz1e1:3.1; s77sz1e2:6.8; s77sz2:5.4; s77sz3:2.4 ;

    %defines start and end times for eeg/2p sampling rates
ictalFront_start_eeg = ictalFront_start*fs_eeg*60;
ictalFront_start_2p = ictalFront_start*fs_2p*60;
ictalFront_end_eeg = ictalFront_end*fs_eeg*60;
ictalFront_end_2p = ictalFront_end*fs_2p*60;

t_2p_ictalFront = F1_ts(ictalFront_start_2p:ictalFront_end_2p)/60;
t_eeg_ictalFront = EEG_ts(ictalFront_start_eeg:ictalFront_end_eeg)/60;


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

%% examine mean traces relative to EEG for ictal front period

figure
plot(t_eeg_ictalFront,EEG(ictalFront_start_eeg:ictalFront_end_eeg)./1000+3,'color',[0.7 0.7 0.7]);
hold on
plot(t_2p_ictalFront,Fc1Gdff_flt1_Mean(ictalFront_start_2p:ictalFront_end_2p),'color', [0 1 0]);
plot(t_2p_ictalFront,Fc1Rdff_flt1_Mean(ictalFront_start_2p:ictalFront_end_2p)-1,'color',[1 0 1]);
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
subplot(1,2,1)
plot(t_eeg_ictalFront,EEG(ictalFront_start_eeg:ictalFront_end_eeg)./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt1,1)
    plot(t_2p_ictalFront,Fc1Gdff_flt1(ii,ictalFront_start_2p:ictalFront_end_2p)-4+ii*6,'color', [0.25 1 0.25]);
end

subplot(1,2,2)
plot(t_eeg_ictalFront,EEG(ictalFront_start_eeg:ictalFront_end_eeg)./1000+3,'color',[0.7 0.7 0.7]);
hold on
for jj=1:size(Fc1Rdff_flt1,1)
    plot(t_2p_ictalFront,Fc1Rdff_flt1(jj,ictalFront_start_2p:ictalFront_end_2p)-4+jj*6,'color',[1 0.25 1]);
end
%ylim([-6.5 5])
%xlim([6 7])
hold off
clear ii jj


%% Investigate traces: filtered trace raster

figure
%Green Ca Raster
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
%title('Non-vGAT (jYCaMP1s) Filtered Fc Transient Raster','fontsize',12)
%xlabel('Times (s)')
%ylabel('Neuron')
xlim([ictalFront_start ictalFront_end])

%Red Ca Raster
ax(1)= subplot(2,1,2);
imagesc(F1_ts/60,size(Fc1Rdff_flt1,1),Fc1Rdff_flt1(yCoorIR,:))
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
%c1.Label.String='dF/F';
clear c1
caxis([-4 10])
pbaspect([3 1 1])
%title('vGAT (jRGECO1a) Filtered Fc Transient Raster','fontsize',12)
%xlabel('Times (s)')
%ylabel('Neuron')
xlim([ictalFront_start ictalFront_end])


%%
figure
%Green Ca Raster
[~,xCoorIG]=sort(aID(kk).xCoorG);
ax(1)= subplot(2,1,1);
imagesc(F1_ts,size(Fc1Gdff_flt1,1),Fc1Gdff_flt1(yCoorIG,:)) %note xCoorIG is used here to reorder cells by x coordiantes
set(gca,'Ydir','normal')
shading flat
colormap(ax(1),jet)
c1=colorbar;
%c1.Label.String='dF/F';
clear c1
caxis([-4 10])
pbaspect([3 1 1])
%title('Non-vGAT (jYCaMP1s) Filtered Fc Transient Raster','fontsize',12)
%xlabel('Times (s)')
%ylabel('Neuron')
xlim([ictalFront_start ictalFront_end])

%Red Ca Raster
[~,xCoorIR]=sort(aID(kk).xCoorR);
ax(2)= subplot(2,1,2);
imagesc(F1_ts,size(Fc1Rdff_flt1,1),Fc1Rdff_flt1(xCoorIR,:)) %note xCoorIG is used here to reorder cells by x coordiantes
set(gca,'Ydir','normal')
shading flat
colormap(ax(2),jet)
c1=colorbar;
%c1.Label.String='dF/F';
clear c1
caxis([-4 10])
pbaspect([3 1 1])
%title('Non-vGAT (jYCaMP1s) Filtered Fc Transient Raster','fontsize',12)
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

% figure
% plot([-999/fs_eeg:1/fs_eeg:0,t_eeg,t_eeg(end)+1/fs_eeg:1/fs_eeg:t_eeg(end)+999/fs_eeg],EEGswdConv)
% hold on
% line([0,t_eeg(end)],[5*std(EEGswdConv(1:MeanSzRecTimeG*fs_eeg)),5*std(EEGswdConv(1:MeanSzRecTimeG*fs_eeg))],'color',[1 0 0])
% line([EEGswd_spikeTimes;EEGswd_spikeTimes],repmat([0;1*10^8],1,length(EEGswd_spikeTimes)),'color',[0 0 1],'linewidth',1)
% 


%% Method to do detection based upon relative theta to gamma power

figure
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
line([EEG_spkTimes_prelim;EEG_spkTimes_prelim],repmat([3;3.75],1,length(EEG_spkTimes_prelim)),'color',[0 1 0],'linewidth',1)
%line([EEG_spkTimes_prelim1;EEG_spkTimes_prelim1],repmat([2.75;3.5],1,length(EEG_spkTimes_prelim1)),'color',[1 0 0],'linewidth',1)
line([EEG_PIS_times2;EEG_PIS_times2],repmat([2;2.75],1,length(EEG_PIS_times2)),'color',[0 0 1],'linewidth',1)

%%
figure
plot(t_eeg_ictalFront,EEG(ictalFront_start_eeg:ictalFront_end_eeg)./1000+3,'color',[0.7 0.7 0.7]);
hold on
%plot(t_2p_ictalFront,Fc1Gdff_flt1_Mean(ictalFront_start_2p:ictalFront_end_2p),'color', [0 1 0]);
%plot(t_2p_ictalFront,Fc1Rdff_flt1_Mean(ictalFront_start_2p:ictalFront_end_2p)-1,'color',[1 0 1]);
plot(t_2p_ictalFront,FGdff_MeanSm(ictalFront_start_2p:ictalFront_end_2p),'color', [0 1 0]);
plot(t_2p_ictalFront,FRdff_MeanSm(ictalFront_start_2p:ictalFront_end_2p)-1,'color',[1 0 1]);
line([MeanSzRecTimeG/60;MeanSzRecTimeG/60],repmat([-4;5],1,2),'color',[0 0.5 0],'LineStyle','--')
line([MeanSzRecTimeR/60;MeanSzRecTimeR/60],repmat([-4;5],1,2),'color',[0.5 0 0.5],'LineStyle','--')
line([MeanCSDRecTimeG/60;MeanCSDRecTimeG/60],repmat([-4;5],1,2),'color',[0 0.5 0],'LineStyle','-.')
line([MeanCSDRecTimeR/60;MeanCSDRecTimeR/60],repmat([-4;5],1,2),'color',[0.5 0 0.5],'LineStyle','-.')
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
plot(F1_ts,Fc1Rdff_flt1_Mean-1,'color',[1 0.8 1]);
plot(F1_ts,FGdff_MeanSm,'color', [0 1 0]);
plot(F1_ts,FRdff_MeanSm-1,'color',[1 0 1]);

% subplot(2,1,2)
% plot(PreSpkRecTimeG{1,3},'color', [0 1 0])
% hold on
% plot(PreSpkRecTimeR{1,3},'color',[1 0 1])
for ii=1:length(PreSpkRecTimeG{1,1})
        subplot(2,1,1)
        line([PreSpkRecTimeG{1,1}(ii);PreSpkRecTimeG{1,1}(ii)],[-.25;.5],'color',[0 0 1],'LineWidth',1)
%         subplot(2,1,2)
%         scatter(ii,PreSpkRecTimeG{1,3}(ii),'MarkerEdgeColor', [0 0 0])
    if ii<=length(PreSpkRecTimeR{1,1})
        subplot(2,1,1)
        line([PreSpkRecTimeR{1,1}(ii);PreSpkRecTimeR{1,1}(ii)],[-1.25;-.5],'color',[0 0 1],'LineWidth',1)
%         subplot(2,1,2)
%         scatter(ii,PreSpkRecTimeR{1,3}(ii),'MarkerEdgeColor', [0 0 0])
    end
    %pause(5)
end
%FPevent=48;%first non-preictal spike insert
% subplot(2,1,1)
%     line([PreSpkRecTimeG{1,1}(FPevent:end);PreSpkRecTimeG{1,1}(FPevent:end)],repmat([-.25;.5],1,length(PreSpkRecTimeG{1,1}(FPevent:end))),'color',[1 0 0],'LineWidth',1)
%     line([PreSpkRecTimeR{1,1}(FPevent:end);PreSpkRecTimeR{1,1}(FPevent:end)],repmat([-1.25;-.5],1,length(PreSpkRecTimeR{1,1}(FPevent:end))),'color',[1 0 0],'LineWidth',1)
%     title('S69PTZ1 top 15%; cuttoff '+string(FPevent-1)+' spikes mean')
subplot(2,1,1)
    line([MeanCaG_PIS_times;MeanCaG_PIS_times],repmat([-.25;.5],1,length(MeanCaG_PIS_times)),'color',[1 0 0],'LineWidth',1)
    line([MeanCaR_PIS_times;MeanCaR_PIS_times],repmat([-1.25;-.5],1,length(MeanCaR_PIS_times)),'color',[1 0 0],'LineWidth',1)
    line([EEG_PIS_times;EEG_PIS_times],repmat([3.75;4.5],1,length(EEG_PIS_times)),'color',[1 0 1],'linewidth',1)
    line([EEG_PIS_times2;EEG_PIS_times2],repmat([2.75;3.5],1,length(EEG_PIS_times2)),'color',[1 0 0],'linewidth',1)
    line([EEGswd_spikeTimes;EEGswd_spikeTimes],repmat([4.75;5.5],1,length(EEGswd_spikeTimes)),'color',[0 0 1],'linewidth',1)
    %title('CaIntegral 1.5 gauss R|G S69PTZ1 top 15%; EEG cuttoff delta through alpha; Ca detect win:'+string(CaDetectWin)+' EEG interval win:'+string(winW)+' EEG spk threshold:'+string(S_th))

subplot(2,1,2)
plot([-999/fs_eeg:1/fs_eeg:0,EEG_ts,EEG_ts(end)+1/fs_eeg:1/fs_eeg:EEG_ts(end)+999/fs_eeg],EEGswdConv/10^7)
hold on
line([0,EEG_ts(end)],[5*std(EEGswdConv(1:MeanSzRecTimeG*fs_eeg)/10^7),5*std(EEGswdConv(1:MeanSzRecTimeG*fs_eeg)/10^7)],'color',[0 0 0])
line([EEGswd_spikeTimes;EEGswd_spikeTimes],repmat([1;2],1,length(EEGswd_spikeTimes)),'color',[1 0 0],'linewidth',.5)


%%

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
plot(F1_ts,Fc1Rdff_flt1_Mean.*10000,'color', [1 0 1]);
% plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_100plus,'color',[1 0 0])
% plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,PSD_100less,'color',[0 0 1])
% plot(PSDwin/2:PSDwin/2:PSDwin*numel(EEGstart)/2,Ratio100.*10000,'color',[0 1 0])
line([MeanSzRecTimeG;MeanSzRecTimeG],[0;50000],'color',[0 .7 0],'LineWidth',.5)
line([MeanSzRecTimeR;MeanSzRecTimeR],[0;50000],'color',[.7 0 .7],'LineWidth',.5)
if death==1
    line([deathTime;deathTime],[0;50000],'color',[.5 0 0],'linewidth',2)
end

    



%%
figure
scatter(RecTimeGsortReG2(isRecruitedG),1:length(RecTimeGsortReG2(isRecruitedG)),10,[0 1 0],'filled')
hold on
scatter(RecTimeRsortReR2(isRecruitedR),1:length(RecTimeRsortReR2(isRecruitedR)),10,[1 0 1],'filled')
title('Ordered with Respect to Recruitment')
xlabel('Times (s)')
ylabel('Neurons')


%% to write figure with sz recruitment times inserted
tempG=Fc1Gdff_flt2; tempR=Fc1Rdff_flt2;

figure
subplot(1,2,1)
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

subplot(1,2,2)
plot(t_eeg_ictalFront,EEG(ictalFront_start_eeg:ictalFront_end_eeg)./1000+3,'color',[0.7 0.7 0.7]);
hold on
for jj=1:size(Fc1Rdff_flt2,1)
    if isRecruitedR(jj)==1
        plot(t_2p_ictalFront,tempR(jj,ictalFront_start_2p:ictalFront_end_2p)+jj*6,'color',[1 0.25 1]);
    else
        plot(t_2p_ictalFront,tempR(jj,ictalFront_start_2p:ictalFront_end_2p)+jj*6,'color',[1 0.75 1]);
    end
    line([RecTimeR_sz2(jj)/60;RecTimeR_sz2(jj)/60],repmat([jj*6;jj*6+1],1,2),'color',[0 0 1],'LineWidth',2)
    line([MeanSzRecTimeR/60;MeanSzRecTimeR/60],repmat([0;750],1,2),'color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',.5)
end
%ylim([-6.5 5])
%xlim([6 7])
hold off
clear ii jj

clear tempG tempR


%% pre ictal indiv ploting
tempG=Fc1Gdff_flt2; tempR=Fc1Rdff_flt2;
%tempG=Fc1Gdff_flt1; tempR=Fc1Rdff_flt1;
figure
subplot(1,2,1)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Gdff_flt2,1)
    plot(F1_ts,tempG(ii,:)+ii*10,'color', [0 1 0]);
    line([RecTimeG_PIS.v1(ii,:);RecTimeG_PIS.v1(ii,:)],repmat([ii*10+1;ii*10+1.75],1,size(RecTimeG_PIS.v1,2)),'color',[0 0 0],'LineWidth',1)
    line([RecTimeG_PIS.v2(ii,:);RecTimeG_PIS.v2(ii,:)],repmat([ii*10+2;ii*10+2.75],1,size(RecTimeG_PIS.v2,2)),'color',[1 0 0],'LineWidth',1)
    line([RecTimeG_PIS.v3(ii,:);RecTimeG_PIS.v3(ii,:)],repmat([ii*10+3;ii*10+3.75],1,size(RecTimeG_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)

end
line([EEG_PIS_times2;EEG_PIS_times2],repmat([0;size(Fc1Gdff_flt2,1)*10],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)

subplot(1,2,2)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=1:size(Fc1Rdff_flt2,1)
    plot(F1_ts,tempR(ii,:)+ii*10,'color', [1 0 1]);
    line([RecTimeR_PIS.v1(ii,:);RecTimeR_PIS.v1(ii,:)],repmat([ii*10+1;ii*10+1.75],1,size(RecTimeR_PIS.v1,2)),'color',[0 0 0],'LineWidth',1)
    line([RecTimeR_PIS.v2(ii,:);RecTimeR_PIS.v2(ii,:)],repmat([ii*10+2;ii*10+2.75],1,size(RecTimeR_PIS.v2,2)),'color',[1 0 0],'LineWidth',1)
    line([RecTimeR_PIS.v3(ii,:);RecTimeR_PIS.v3(ii,:)],repmat([ii*10+3;ii*10+3.75],1,size(RecTimeR_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)

end
line([EEG_PIS_times2;EEG_PIS_times2],repmat([0;size(Fc1Rdff_flt2,1)*10],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)

%% Z-scored individual trace plotting
tempG=Fc1Gdff_flt2; tempR=Fc1Rdff_flt2;
tempGmean1=mean(tempG(:,(EEG_PIS_times2(1)-10)*fs_2p:MeanSzRecTimeG*fs_2p),2);
tempRmean1=mean(tempR(:,(EEG_PIS_times2(1)-10)*fs_2p:MeanSzRecTimeR*fs_2p),2);
tempGstd1=std(tempG(:,(EEG_PIS_times2(1)-10)*fs_2p:MeanSzRecTimeG*fs_2p),[],2);
tempRstd1=std(tempR(:,(EEG_PIS_times2(1)-10)*fs_2p:MeanSzRecTimeR*fs_2p),[],2);
tempGstd=sort(tempG(:,(EEG_PIS_times2(1)-10)*fs_2p:MeanSzRecTimeG*fs_2p),2);
tempRstd=sort(tempR(:,(EEG_PIS_times2(1)-10)*fs_2p:MeanSzRecTimeR*fs_2p),2);
tempGmean=mean(tempGstd(:,1:size(tempGstd,2)/2),2);
tempRmean=mean(tempRstd(:,1:size(tempRstd,2)/2),2);
tempGstd=std(tempGstd(:,1:size(tempGstd,2)/2),[],2);
tempRstd=std(tempRstd(:,1:size(tempRstd,2)/2),[],2);

tempGz=(tempG-tempGmean)./tempGstd;
tempRz=(tempR-tempRmean)./tempRstd;

%tempG=Fc1Gdff_flt1; tempR=Fc1Rdff_flt1;
figure
subplot(1,2,1)
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
cc=1;%counter
for ii=[149,153]%list of ROIs
    plot(F1_ts,tempG(ii,:)+cc*10-2,'color', [0 .5 0]);
    plot(F1_ts,tempGz(ii,:)+cc*10,'color', [0 1 0]);
    line([RecTimeG_PIS.v1(ii,:);RecTimeG_PIS.v1(ii,:)],repmat([cc*10+1;cc*10+1.75],1,size(RecTimeG_PIS.v1,2)),'color',[0 0 0],'LineWidth',1)
    line([RecTimeG_PIS.v2(ii,:);RecTimeG_PIS.v2(ii,:)],repmat([cc*10+2;cc*10+2.75],1,size(RecTimeG_PIS.v2,2)),'color',[1 0 0],'LineWidth',1)
    line([RecTimeG_PIS.v3(ii,:);RecTimeG_PIS.v3(ii,:)],repmat([cc*10+3;cc*10+3.75],1,size(RecTimeG_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)
    line([0;EEG_ts(end)],[cc*10+tempGmean(ii)+4*tempGstd(ii)-2;cc*10+tempGmean(ii)+4*tempGstd(ii)-2],'color',[0.8 0.8 0.8],'LineWidth',.5)
    cc=cc+1;
end
line([EEG_PIS_times2;EEG_PIS_times2],repmat([0;cc*10],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)
clear ii cc


subplot(1,2,2)
cc=1;%counter
plot(EEG_ts,EEG./1000+3,'color',[0.7 0.7 0.7]);
hold on
for ii=[21,36]%list of ROIs
    plot(F1_ts,tempR(ii,:)+cc*10-2,'color', [.5 0 .5]);
    plot(F1_ts,tempRz(ii,:)+cc*10,'color', [1 0 1]);
    line([RecTimeR_PIS.v1(ii,:);RecTimeR_PIS.v1(ii,:)],repmat([cc*10+1;cc*10+1.75],1,size(RecTimeR_PIS.v1,2)),'color',[0 0 0],'LineWidth',1)
    line([RecTimeR_PIS.v2(ii,:);RecTimeR_PIS.v2(ii,:)],repmat([cc*10+2;cc*10+2.75],1,size(RecTimeR_PIS.v2,2)),'color',[1 0 0],'LineWidth',1)
    line([RecTimeR_PIS.v3(ii,:);RecTimeR_PIS.v3(ii,:)],repmat([cc*10+3;cc*10+3.75],1,size(RecTimeR_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)
    line([0;EEG_ts(end)],[cc*10+tempRmean(ii)+4*tempRstd(ii)-2;cc*10+tempRmean(ii)+4*tempRstd(ii)-2],'color',[0.8 0.8 0.8],'LineWidth',.5)
    cc=cc+1;
end
line([EEG_PIS_times2;EEG_PIS_times2],repmat([0;cc*10],1,size(EEG_PIS_times2,2)),'color',[0.8 0.8 0.8],'LineWidth',.5)
clear ii cc

%clear tempG* tempR*

