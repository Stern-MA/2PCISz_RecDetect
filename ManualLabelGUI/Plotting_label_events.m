%This is a companion script to the labeling GUI that allows you to plot
%your data against your labels.

%% Plot Mean Trace selected times PIS (GREEN)

recTimes=cell2mat(label_struct.G_preictal_mean.labels);
recTimes=recTimes(:,1);
figure
plot(t_2p,mean(data.Fc1Gdff_flt2),'color',[0 1 0])
hold on
line([recTimes recTimes]',repmat([.5;1.25],1,length(recTimes)),'color',[0 0 0])
line([data.EEG_PIS_times3;data.EEG_PIS_times3],repmat([1.5;2.25],1,length(data.EEG_PIS_times3)),'color',[1 0 0])

%% Plot Indiv Traces selected Times PIS (GREEN)

figure
hold on
cc=1;
for ii=1:length(label_struct.G_preictal.labels)
    %will look to see if a cell labeled and if so does it have at least one
    %label (i,.e is not a cell that was decided to not be recruited in
    %labeling)
    if iscell(label_struct.G_preictal.labels{ii}) && any(~cellfun(@isempty,label_struct.G_preictal.labels{ii}))
        recTimes=cell2mat(label_struct.G_preictal.labels{ii});
        recTimes=recTimes(:,1);
        plot(t_2p,data.Fc1Gdff_flt2(ii,:)+cc*10,'color',[0 1 0])
        line([recTimes recTimes]',repmat([cc*10+1;cc*10+1.75],1,length(recTimes)),'color',[0 0 0],'LineWidth',1)
        line([data.RecTimeG_PIS.v3(ii,:);data.RecTimeG_PIS.v3(ii,:)],repmat([cc*10+2;cc*10+2.75],1,size(data.RecTimeG_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)
        cc=cc+1;
    elseif ~isempty(label_struct.G_preictal.labels{ii}) %will plot the labeled unrecruited cell against the algorithim event detection
        plot(t_2p,data.Fc1Gdff_flt2(ii,:)+cc*10,'color',[0 1 0])
        line([data.RecTimeG_PIS.v3(ii,:);data.RecTimeG_PIS.v3(ii,:)],repmat([cc*10+2;cc*10+2.75],1,size(data.RecTimeG_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)
        cc=cc+1;
    end
end

%% Plot Mean Trace selected times Sz (GREEN)

recTimes=label_struct.G_seizure_mean.labels(1,1);
figure
plot(t_2p,mean(data.Fc1Gdff_flt2),'color',[0 1 0])
hold on
line([recTimes recTimes]',repmat([.5;1.25],1,length(recTimes)),'color',[0 0 0])
line([data.MeanSzRecTimeG;data.MeanSzRecTimeG],[1.5;2.25],'color',[1 0 0])
xlim([data.MeanSzRecTimeG-10 data.MeanSzRecTimeG+10])


%% Plot Indiv Traces selected Times Seizure (GREEN)
figure
hold on
cc=1;
for ii=1:length(label_struct.G_seizure.labels)
    if ~isempty(label_struct.G_seizure.labels{ii})
        recTimes=label_struct.G_seizure.labels{ii}(1,1);
        if data.isRecruitedG(ii)==1
            plot(t_2p,data.Fc1Gdff_flt2(ii,:)+cc*10,'color',[0 1 0])
        else
            plot(t_2p,data.Fc1Gdff_flt2(ii,:)+cc*10,'color',[0.6 1 0.6])
        end
        line([recTimes recTimes]',repmat([cc*10+1;cc*10+1.75],1,2),'color',[0 0 0],'LineWidth',1)
        line([data.RecTimeG_sz2(ii);data.RecTimeG_sz2(ii)],repmat([cc*10+2;cc*10+2.75],1,2),'color',[0 0 1],'LineWidth',1)
        cc=cc+1;
    end
end
xlim([mean(data.RecTimeG_sz2)-10 mean(data.RecTimeG_sz2)+10])


%% Plot Mean Trace selected times CSD (GREEN)

recTimes=label_struct.G_csd_mean.labels(1,1);
figure
plot(t_2p,mean(data.Fc1Gdff_flt2),'color',[0 1 0])
hold on
line([recTimes recTimes]',repmat([.5;1.25],1,length(recTimes)),'color',[0 0 1])
line([data.MeanCSDRecTimeG;data.MeanCSDRecTimeG],[1.5;2.25],'color',[1 0 0])
xlim([data.MeanCSDRecTimeG-10 data.MeanCSDRecTimeG+10])


%% Plot Indiv Traces selected Times CSD (GREEN)
figure
hold on
cc=1;
for ii=1:length(label_struct.G_csd.labels)
    if ~isempty(label_struct.G_csd.labels{ii})
        recTimes=label_struct.G_csd.labels{ii}(1,1);
        if data.isRecruitedG_csd(ii)==1
            plot(t_2p,data.Fc1Gdff_flt2(ii,:)+cc*10,'color',[0 1 0])
        else
            plot(t_2p,data.Fc1Gdff_flt2(ii,:)+cc*10,'color',[0.6 1 0.6])
        end
        line([recTimes recTimes]',repmat([cc*10+1;cc*10+1.75],1,2),'color',[0 0 0],'LineWidth',1)
        line([data.RecTimeG_csd(ii);data.RecTimeG_csd(ii)],repmat([cc*10+2;cc*10+2.75],1,2),'color',[0 0 1],'LineWidth',1)
        cc=cc+1;
    end
end
xlim([mean(data.RecTimeG_csd)-10 mean(data.RecTimeG_csd)+10])


%% Plot Mean Trace selected times death (GREEN)

recTimes=label_struct.G_death_mean.labels(1,1);
figure
plot(t_2p,mean(data.Fc1Gdff_flt2),'color',[0 1 0])
hold on
line([recTimes recTimes]',repmat([.5;1.25],1,length(recTimes)),'color',[0 0 1])
line([data.MeanCSDRecTimeG;data.MeanCSDRecTimeG],[1.5;2.25],'color',[1 0 0])
xlim([data.MeanCSDRecTimeG-10 data.MeanCSDRecTimeG+10])


%% Plot Indiv Traces selected Times death (GREEN)
figure
hold on
cc=1;
for ii=1:length(label_struct.G_death.labels)
    if ~isempty(label_struct.G_death.labels{ii})
        recTimes=label_struct.G_death.labels{ii}(1,1);
        if data.isRecruitedG_csd(ii)==1
            plot(t_2p,data.Fc1Gdff_flt2(ii,:)+cc*10,'color',[0 1 0])
        else
            plot(t_2p,data.Fc1Gdff_flt2(ii,:)+cc*10,'color',[0.6 1 0.6])
        end
        line([recTimes recTimes]',repmat([cc*10+1;cc*10+1.75],1,2),'color',[0 0 0],'LineWidth',1)
        line([data.RecTimeG_csd(ii);data.RecTimeG_csd(ii)],repmat([cc*10+2;cc*10+2.75],1,2),'color',[0 0 1],'LineWidth',1)
        cc=cc+1;
    end
end
xlim([mean(data.RecTimeG_csd)-10 mean(data.RecTimeG_csd)+10])






%% Plot Mean Trace selected times PIS (RED)

recTimes=cell2mat(label_struct.R_preictal_mean.labels);
recTimes=recTimes(:,1);
figure
plot(t_2p,mean(data.Fc1Rdff_flt2),'color',[1 0 1])
hold on
line([recTimes recTimes]',repmat([.5;1.25],1,length(recTimes)),'color',[0 0 0])
line([data.EEG_PIS_times3;data.EEG_PIS_times3],repmat([1.5;2.25],1,length(data.EEG_PIS_times3)),'color',[1 0 0])

%% Plot Indiv Traces selected Times PIS (RED)
figure
hold on
cc=1;
for ii=1:length(label_struct.R_preictal.labels)
    if iscell(label_struct.R_preictal.labels{ii}) && any(~cellfun(@isempty,label_struct.R_preictal.labels{ii}))
        recTimes=cell2mat(label_struct.R_preictal.labels{ii});
        recTimes=recTimes(:,1);
        plot(t_2p,data.Fc1Rdff_flt2(ii,:)+cc*10,'color',[1 0 1])
        line([recTimes recTimes]',repmat([cc*10+1;cc*10+1.75],1,length(recTimes)),'color',[0 0 0],'LineWidth',1)
        line([data.RecTimeR_PIS.v3(ii,:);data.RecTimeR_PIS.v3(ii,:)],repmat([cc*10+2;cc*10+2.75],1,size(data.RecTimeR_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)
        cc=cc+1;
    elseif ~isempty(label_struct.R_preictal.labels{ii})
        plot(t_2p,data.Fc1Rdff_flt2(ii,:)+cc*10,'color',[1 .5 1])
        line([data.RecTimeR_PIS.v3(ii,:);data.RecTimeR_PIS.v3(ii,:)],repmat([cc*10+2;cc*10+2.75],1,size(data.RecTimeR_PIS.v3,2)),'color',[0 0 1],'LineWidth',1)
        cc=cc+1;
    end
end


%% Plot Mean Trace selected times Sz (RED)

recTimes=label_struct.R_seizure_mean.labels(1,1);
figure
plot(t_2p,mean(data.Fc1Rdff_flt2))
hold on
line([recTimes recTimes]',repmat([.5;1.25],1,length(recTimes)),'color',[1 0 1])
line([data.MeanSzRecTimeR;data.MeanSzRecTimeR],[1.5;2.25],'color',[1 0 0])
xlim([data.MeanSzRecTimeR-10 data.MeanSzRecTimeR+10])


%% Plot Indiv Traces selected Times Seizure (RED)
figure
hold on
cc=1;
for ii=1:length(label_struct.R_seizure.labels)
    if ~isempty(label_struct.R_seizure.labels{ii})
        recTimes=label_struct.R_seizure.labels{ii}(1,1);
        if data.isRecruitedR(ii)==1
            plot(t_2p,data.Fc1Rdff_flt2(ii,:)+cc*10,'color',[1 0 1])
        else
            plot(t_2p,data.Fc1Rdff_flt2(ii,:)+cc*10,'color',[1 0.6 1])
        end
        line([recTimes recTimes]',repmat([cc*10+1;cc*10+1.75],1,2),'color',[0 0 0],'LineWidth',1)
        line([data.RecTimeR_sz2(ii);data.RecTimeR_sz2(ii)],repmat([cc*10+2;cc*10+2.75],1,2),'color',[0 0 1],'LineWidth',1)
        cc=cc+1;
    end
end
xlim([mean(data.RecTimeR_sz2)-10 mean(data.RecTimeR_sz2)+10])

%% Plot Mean Trace selected times CSD (RED)

recTimes=label_struct.R_csd_mean.labels(1,1);
figure
plot(t_2p,mean(data.Fc1Rdff_flt2))
hold on
line([recTimes recTimes]',repmat([.5;1.25],1,length(recTimes)),'color',[1 0 1])
line([data.MeanCSDRecTimeR;data.MeanCSDRecTimeR],[1.5;2.25],'color',[1 0 0])
xlim([data.MeanCSDRecTimeR-10 data.MeanCSDRecTimeR+10])


%% Plot Indiv Traces selected Times CSD (RED)
figure
hold on
cc=1;
for ii=1:length(label_struct.R_csd.labels)
    if ~isempty(label_struct.R_csd.labels{ii})
        recTimes=label_struct.R_csd.labels{ii}(1,1);
        if data.isRecruitedR_csd(ii)==1
            plot(t_2p,data.Fc1Rdff_flt2(ii,:)+cc*10,'color',[1 0 1])
        else
            plot(t_2p,data.Fc1Rdff_flt2(ii,:)+cc*10,'color',[1 0.6 1])
        end
        line([recTimes recTimes]',repmat([cc*10+1;cc*10+1.75],1,2),'color',[0 0 0],'LineWidth',1)
        line([data.RecTimeR_csd(ii);data.RecTimeR_csd(ii)],repmat([cc*10+2;cc*10+2.75],1,2),'color',[0 0 1],'LineWidth',1)
        cc=cc+1;
    end
end
xlim([mean(data.RecTimeR_csd)-10 mean(data.RecTimeR_csd)+10])


%% Plot Mean Trace selected times death (RED)

recTimes=label_struct.R_death_mean.labels(1,1);
figure
plot(t_2p,mean(data.Fc1Rdff_flt2))
hold on
line([recTimes recTimes]',repmat([.5;1.25],1,length(recTimes)),'color',[1 0 1])
line([data.MeanCSDRecTimeR;data.MeanCSDRecTimeR],[1.5;2.25],'color',[1 0 0])
xlim([data.MeanCSDRecTimeR-10 data.MeanCSDRecTimeR+10])


%% Plot Indiv Traces selected Times death (RED)
figure
hold on
cc=1;
for ii=1:length(label_struct.R_death.labels)
    if ~isempty(label_struct.R_death.labels{ii})
        recTimes=label_struct.R_death.labels{ii}(1,1);
        if data.isRecruitedR_csd(ii)==1
            plot(t_2p,data.Fc1Rdff_flt2(ii,:)+cc*10,'color',[1 0 1])
        else
            plot(t_2p,data.Fc1Rdff_flt2(ii,:)+cc*10,'color',[1 0.6 1])
        end
        line([recTimes recTimes]',repmat([cc*10+1;cc*10+1.75],1,2),'color',[0 0 0],'LineWidth',1)
        line([data.RecTimeR_csd(ii);data.RecTimeR_csd(ii)],repmat([cc*10+2;cc*10+2.75],1,2),'color',[0 0 1],'LineWidth',1)
        cc=cc+1;
    end
end
xlim([mean(data.RecTimeR_csd)-10 mean(data.RecTimeR_csd)+10])
