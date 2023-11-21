% GUI for manual labeling of seizure data
%
%Input is the output of the SzRecTime_processing.m script
%
%The purpose of this GUI is to enable allow for manual labeling of your 
%data to compare to algorithm's performance.


% Note that this script is not designed to rewrite over previously scored
% data, but rather to append new scores to the end. The user will need to
% manual delete a structure if they want to correct for someting previosu
% scored incorrectly.

%label_events: labeling script to save manual annotations marking
% timing of preictal spikes, seizure start times, etc

% Version 231116 Eric R. Cole and Matthew A. Stern, Emory University
% Contact: matthew.a.stern@emory.edu or matt@matthewastern.com
%
% If using this code please cite our work.


%% Load data files:

file = 'S69ptz1'; %file name located in directory
szdata = load(sprintf('%s.mat',file));
subfiles = fieldnames(szdata);

data = szdata.(subfiles{1});

i1 = imread('important_image_1.jpg');
i2 = imread('important_image_2.jpeg');

%% Load Working Structure
%change the variable to true if this is the first timing running the script
%to score a given recording. If a structure already exists from a
%previously partially scored recording that you wish to continue with here,
%change this varaible to false.
first_time = false; %if first run through of new recording

if first_time
    label_struct=struct();  
    save(sprintf('%s_labeled.mat',file),'label_struct');
else
    file_to_load = sprintf('%s_labeled.mat',file);
    load(file_to_load);
    save(sprintf('%s_labeled_copy.mat',file),'label_struct');
end
    
%% Set plotting settings:

w_size = 5; %seconds to view at a time
pct_overlap = .1; %fraction of overlap between consecutive windows
y_pct = [0 0.9]; %sets y limit range (percentiles of mean trace)

filt_type = '1Hz'; %'1Hz' (default), '5Hz', or 'None'

show_labels = true; %display labels on time series after entering

%file_type = 'new'; %'new' to start new file, 'continue' to load existing

w_overlap = pct_overlap*w_size;

%Input settings:
q = 'quit'; d = 'done'; r = 'redo'; n = 'next'; j = 'jump'; nc = 'next cell';% To be implimented in future: b = 'back'

input_s = sprintf('Select a command:\nq = %s, d = %s, r = %s, n = %s, nc = %s, j = %s\n\n',q,d,r,n,nc,j);

%% Plot overview of data
%use this to determine start and end times of the event periods to enter in
%the below sections

channel = 'G';

figure

if strcmp(channel, 'G')
    data_2p = data.Fc1Gdff_flt1;
elseif strcmp(channel, 'R')
    data_2p = data.Fc1Rdff_flt1;
end

t_2p = (1:size(data_2p,2))/data.fs_2p; 
t_eeg = (1:size(data.EEG,2))/data.fs_eeg; 
mean_dff = mean(data_2p);

subplot(2,1,1)
plot(t_eeg,data.EEG)

subplot(2,1,2)
plot(t_2p,mean_dff)
xlabel('Time (s)')


%% label mean trace preictal spikes:
%Run this section to score the mean trace preictal spikes input the start 
%and end times for this region below. Make sureto run the section after to 
%save the data.

channel = 'G';%select channel for analysis
start_time = 200;  %start time in seconds
end_time = 500;  %end time in seconds

file_to_load = sprintf('%s_labeled.mat',file);
load(file_to_load);
save(sprintf('%s_labeled_copy.mat',file),'label_struct');

if strcmp(channel, 'G')
    data_2p = data.Fc1Gdff_flt2;
elseif strcmp(channel, 'R')
    data_2p = data.Fc1Rdff_flt2;
end

t_2p = (1:size(data_2p,2))/data.fs_2p; 
t_eeg = (1:size(data.EEG,2))/data.fs_eeg;
mean_dff = mean(data_2p);
ylims = [prctile(mean_dff,100*y_pct(1)), prctile(mean_dff,100*y_pct(2))];

if isfield(label_struct,sprintf('%s_preictal_mean',channel))
    mean_labeled = label_struct.(sprintf('%s_preictal_mean',channel));
    labels = mean_labeled.labels;
    window_inds = mean_labeled.window_inds;
    window_times = mean_labeled.window_times;
    window_ind = mean_labeled.window_ind;
else
    labels = {};
    window_inds = {};
    window_times = {};
	window_ind = 1;
end

norm_gs = @(x) (x-min(x))/(max(x)-min(x));
eeg_norm = norm_gs(data.EEG);
max_dff = max(mean_dff) ;
range_dff = max(mean_dff) - min(mean_dff);


figure('Position', [200 200 1000 800])

subplot(4,1,1)
title(sprintf('%s; preictal spikes; mean trace; %s channel',file,channel));

subplot(4,1,2)
plot(t_eeg, data.EEG/1000 + max_dff+1,'k')
%plot(t_eeg, (eeg_norm*range_dff)*2 + max_dff - 2,'k')
hold on
plot(t_2p,mean_dff,'color',[0 0 .7])
ylim([-.5 max(mean_dff)+3])
if start_time-100<=0
    xlim([0 end_time+100])
else
    xlim([start_time-100 end_time+100])
end

subplot(4,1,[3 4])
plot(t_2p,mean_dff)
prev_labels_to_plot = cell2mat(labels);
hold on
if show_labels && ~(isempty(prev_labels_to_plot))
    lb_plot = scatter(prev_labels_to_plot(:,1),prev_labels_to_plot(:,2),200,'k*');
end
hold on

loop_flag = true;
while(loop_flag)
    
    
    subplot(4,1,[3 4])
    hold on
    t_start = start_time + (w_size)*(window_ind-1);
    t_end = t_start + w_size;
    xlim([t_start-w_overlap,t_end+w_overlap])
    ylim(ylims + [-0.5 1])
    pbaspect([1 1.25 1])
    
    b1 = xline(t_start, '--');
    b2 = xline(t_end, '--');
    
    subplot(4,1,2)
    hold on
    m1 = xline(t_start);
    m2 = xline(t_end);
    hold off 
    
%     t1 = round(t_start*fs_2p);
%     t2 = round(t_end*fs_2p);
    
    
    temp_labels = ginput();
    
    if show_labels && ~(isempty(temp_labels))
        lb_plot = scatter(temp_labels(:,1),temp_labels(:,2),200,'k*');
    end
    
    hold off
    
    action = input(input_s);
    if isempty(action)
        action = 'next';
    end
    
    switch action
        case 'quit'
            window_inds = [window_inds; {window_ind}];
            loop_flag = false;
        case 'done'
            window_inds = [window_inds; {window_ind}];
            loop_flag = false;
        case 'redo'
            delete(lb_plot)
            temp_labels = [];
            
            delete(m1); delete(m2);
            continue
        case 'next'
            window_inds = [window_inds; {window_ind}];
            window_ind = window_ind + 1;
        case 'jump'
            window_inds = [window_inds; {window_ind}];
            new_t = input('Enter a time (s) to jump to: ');
            window_ind = floor((new_t-start_time)/w_size);
        otherwise
            window_inds = [window_inds; {window_ind}];
            window_ind = window_ind + 1;
    end
    labels = [labels; {temp_labels}];
    window_times = [window_times; {[t_start, t_end]}];
    
    delete(m1); delete(m2);
end

%% Save labels for mean trace, preictal:
%Run this section to save the mean trace preictal output from above. Note
%that if you do not run this section after quiting the above GUI that data
%will not be saved.

%Quality Control
%remove labels times that fall outside of the window within which they were
%selected
labelsNew=cell(size(labels));
for ii=1:length(labels)
    if ~isempty(labels{ii})
        labelsNew{ii}=labels{ii}(and(labels{ii}(:,1)<=window_times{ii}(2),labels{ii}(:,1)>=window_times{ii}(1)),:);
    else
        labelsNew{ii}=[];
    end
end

%Saving structure
mean_labeled = struct();
mean_labeled.labels = labelsNew;
mean_labeled.window_inds = window_inds;
mean_labeled.window_times = window_times;
mean_labeled.start_time = start_time;
mean_labeled.w_size = w_size;
mean_labeled.window_ind = window_ind;

if strcmp(channel, 'G')
    label_struct.G_preictal_mean = mean_labeled;
elseif strcmp(channel, 'R')
    label_struct.R_preictal_mean = mean_labeled;
end

save(sprintf('%s_labeled.mat',file),'label_struct');
disp('saved preictal mean to struct')

%% Label preictal period for individual cells:
%Run this section to score preictal spikes on individual traces. Make sure
%to run the section after to save the data.

channel = 'G';
start_time = 200;  %start time in seconds
end_time = 500;  %start time in seconds

fixcellorder=0;
if strcmp(channel, 'G')
%s69-G
    fixedcellorder=[8,12,33,45,32,28,43,26,60,3,52,22,37,27,30,39,2,16,21,17,4,36,67,38,58,13,53,44,64,1,47,19,50,29,46,40,34,57,20,42,62,18,7,14,56,5,15,61,65,24,11,35,6,49,23,63,66,9,51,31,41,55,10,54,25,59,48];
elseif strcmp(channel, 'R')
%s69-R
    fixedcellorder=[18,11,10,38,31,40,15,34,20,24,32,26,5,22,19,8,39,30,17,23,12,36,3,9,16,2,25,27,28,29,7,6,13,21,35,1,33,4,14,37];
end

file_to_load = sprintf('%s_labeled.mat',file);
load(file_to_load);
save(sprintf('%s_labeled_copy.mat',file),'label_struct');

if strcmp(channel, 'G')
    data_2p = data.Fc1Gdff_flt2;
elseif strcmp(channel, 'R')
    data_2p = data.Fc1Rdff_flt2;
end

t_2p = (1:size(data_2p,2))/data.fs_2p; 
t_eeg = (1:size(data.EEG,2))/data.fs_eeg;
mean_dff = mean(data_2p);
ylims = [prctile(mean_dff,100*y_pct(1)), prctile(mean_dff,100*y_pct(2))];

n_cells = size(data_2p,1);

if isfield(label_struct,sprintf('%s_preictal',channel))
    mean_labeled = label_struct.(sprintf('%s_preictal',channel));
    labels = mean_labeled.labels;
    window_inds = mean_labeled.window_inds;
    window_times = mean_labeled.window_times;
    window_ind = mean_labeled.window_ind;
    cell_ind = mean_labeled.cell_ind;
    cell_order=mean_labeled.cell_order;
else
    window_ind = 1;
    cell_ind = 1;
    
    if fixcellorder==1
        cell_order=fixedCellOrder;
    else
        cell_order = randperm(n_cells);   
    end
    
    labels = cell(n_cells,1);
    window_inds = cell(n_cells,1);
    window_times = cell(n_cells,1);
end

norm_gs = @(x) (x-min(x))/(max(x)-min(x));
eeg_norm = norm_gs(data.EEG);
max_dff = max(mean_dff);
range_dff = max(mean_dff) - min(mean_dff);

figure('Position', [200 200 1000 800])
subplot(5,1,1)
mouse1 = image(i1,'xdata',[0 .15],'ydata',[0.2 0.7]);
hold on; 
ylim([0 1])
hold on;
xlim([0 1])
image(i2,'xdata',[0.9 1],'ydata',[0 0.9])
yticks([0.5])
yticklabels({''})
xticklabels([])
m1_pos = 0;

subplot(5,1,2)
plot(t_eeg, data.EEG/1000 + max_dff+1,'k')
%plot(t_eeg, (eeg_norm*range_dff)*2 + max_dff - 2,'k')
hold on
plot(t_2p,mean_dff,'color',[0 0 .7])
ylim([-.5 max(mean_dff)+3])
if start_time-100<=0
    xlim([0 end_time+100])
else
    xlim([start_time-100 end_time+100])
end

for kk = cell_ind:n_cells
    loop_flag = true;
    subplot(5,1,1)
    title(sprintf('%s; cell %d; %d/%d done (%2.1f pct.); preictal; %s channel',file,cell_order(kk),kk-1,length(cell_order),100*kk/length(cell_order),channel));
    
    progress = 100*kk/length(cell_order);
    
    subplot(5,1,3)
    plot(t_2p,data_2p(cell_order(kk),:))
    if start_time-100<=0
        xlim([0 end_time+100])
    else
        xlim([start_time-100 end_time+100])
    end
    ylim([-.5 max([max(mean_dff) max(data_2p(cell_order(kk),:))])])
    
    subplot(5,1,[4 5])
    hold off
    plot(t_2p,data_2p(cell_order(kk),:));
%     if load_file
%         prev_labels_to_plot = cell2mat(labels);
%         hold on
%         if show_labels && ~(isempty(prev_labels_to_plot))
%             lb_plot = scatter(temp_labels(:,1),temp_labels(:,2),200,'k*');
%         end
%     end
    hold on
    while(loop_flag)
        
        subplot(5,1,[4 5])
        hold on
        t_start = start_time + (w_size)*(window_ind-1);
        t_end = t_start + w_size;
        xlim([t_start-w_overlap,t_end+w_overlap])
        ylim(ylims + [-0.5 3])
        pbaspect([1 2 1])
        
        b1 = xline(t_start, '--');
        b2 = xline(t_end, '--');
        
        subplot(5,1,3)
        
        mm1 = xline(t_start);
        mm2 = xline(t_end);
        
        subplot(5,1,2)
        
        m1 = xline(t_start);
        m2 = xline(t_end);
        
        temp_labels = ginput;
        
        if show_labels && ~(isempty(temp_labels))
            subplot(5,1,[4 5])
            lb_plot = scatter(temp_labels(:,1),temp_labels(:,2),200,'k*');
        end
        
        hold off
        
        action = input(input_s);
        if isempty(action)
            action = 'next';
        end
        
        switch action
            case 'quit'
                window_inds{cell_order(kk)} = [window_inds{cell_order(kk)}; {window_ind}];
                loop_flag = false;
                quit_flag = true;
            case 'done'
                window_inds{cell_order(kk)} = [window_inds{cell_order(kk)}; {window_ind}];
                loop_flag = false;
                quit_flag = true;
            case 'redo'
                delete(lb_plot)
                temp_labels = [];
                quit_flag = false;
                delete(m1); delete(m2); delete(mm1); delete(mm2);
                continue
            case 'next'
                window_inds{cell_order(kk)} = [window_inds{cell_order(kk)}; {window_ind}];
                window_ind = window_ind + 1;
                quit_flag = false;
            case 'next cell'
                window_inds{cell_order(kk)} = [window_inds{cell_order(kk)}; {window_ind}];
                quit_flag = false;
                delete(m1); delete(m2); delete(mm1); delete(mm2);
                break
            case 'jump'
                window_inds{cell_order(kk)} = [window_inds{cell_order(kk)}; {window_ind}];
                new_t = input('Enter a time (s) to jump to: ');
                window_ind = floor((new_t-start_time)/w_size);
                quit_flag = false;
            otherwise
                window_inds{cell_order(kk)} = [window_inds{cell_order(kk)}; {window_ind}];
                window_ind = window_ind + 1;
                quit_flag = false;
        end
        labels{cell_order(kk)} = [labels{cell_order(kk)}; {temp_labels}];
        window_times{cell_order(kk)} = [window_times{cell_order(kk)}; {[t_start, t_end]}];
        
        delete(m1); delete(m2); delete(mm1); delete(mm2);
    end
    if quit_flag
       break 
    end
    cell_ind = cell_ind + 1;
    window_ind = 1;
    
    if m1_pos < 0.9
        delete(mouse1)
        
        m1_pos = progress/100;
        
        subplot(5,1,1)
        hold on
        mouse1 = image(i1,'xdata',[0 .15]+m1_pos,'ydata',[0.2 0.7]);
        
    end
end

%% Save labels for preictal events per cell:
%Run this section to save the individual trace preictal output from above. 
%Note that if you do not run this section after quiting the above GUI that 
%data will not be saved.

%Quality Control
%remove labels times that fall outside of the window within which they were
%selected
labelsNew=cell(size(labels));
for ii=1:length(labels)
    if ~isempty(labels{ii})
        labelsNew{ii}=cell(size(labels{ii}));
        for jj=1:length(labels{ii})
            if ~isempty(labels{ii}{jj})
                labelsNew{ii}{jj}=labels{ii}{jj}(and(labels{ii}{jj}(:,1)<=window_times{ii}{jj}(2),labels{ii}{jj}(:,1)>=window_times{ii}{jj}(1)),:);
            else
                labelsNew{ii}{jj}=[];
            end
        end
    end
end


cells_labeled = struct();
cells_labeled.labels = labelsNew;
cells_labeled.window_inds = window_inds;
cells_labeled.window_times = window_times;
cells_labeled.start_time = start_time;
cells_labeled.w_size = w_size;
cells_labeled.window_ind = window_ind;
cells_labeled.cell_ind = cell_ind; 
cells_labeled.cell_order = cell_order;

if strcmp(channel, 'G')
    label_struct.G_preictal = cells_labeled;
elseif strcmp(channel, 'R')
    label_struct.R_preictal = cells_labeled;
end

save(sprintf('%s_labeled.mat',file),'label_struct');

disp('saved preictal indiv to struct')

%% Label mean trace for single event (seizure start time, CSD, wave of death)
%Run this section to score single event times on the mean trace. Make sure
%to run the section after to save the data.

channel = 'G';
mode = 'death'; %designate the event type: 'seizure','csd', or 'death'
start_time = 420; %enter start time of window
end_time = 500; %enter end time of window

figure('Position',[200 200 1000 800])

if strcmp(channel, 'G')
    data_2p = data.Fc1Gdff_flt2;
elseif strcmp(channel, 'R')
    data_2p = data.Fc1Rdff_flt2;
end

t_2p = (1:size(data_2p,2))/data.fs_2p; 
t_eeg = (1:size(data.EEG,2))/data.fs_eeg; 
mean_dff = mean(data_2p);

xlims = [start_time, end_time];

subplot(2,1,1)
plot(t_eeg,data.EEG)
xlim(xlims)

subplot(2,1,2)
plot(t_2p,mean_dff)
hold on
xlabel('Time (s)')
xlim(xlims)
ylim([-.5 max(mean_dff)+1])

temp_labels = ginput();

scatter(temp_labels(1),temp_labels(2),200,'k*');

%% Save mean label for single event
%Run this section to save the mean trace event output from above. 
%Note that if you do not run this section after quiting the above GUI that 
%data will not be saved.

event_labeled = struct();
event_labeled.labels = temp_labels;
event_labeled.start_time = start_time;
event_labeled.end_time = end_time;

label_struct.(sprintf('%s_%s_mean',channel,mode)) = event_labeled;

save(sprintf('%s_labeled.mat',file),'label_struct');

disp(sprintf('saved %s mean to struct',mode))


%% Label individual cells for single event (seizure start time, CSD, wave of death)
%Run this section to score single event times on individual traces. Make 
%sure to run the section after to save the data.

channel = 'G';
mode = 'seizure'; %designate the event type: 'seizure','csd', or 'death'
start_time = 420; %enter start time of window
end_time = 520; %enter ends time of window

file_to_load = sprintf('%s_labeled.mat',file);
load(file_to_load);
save(sprintf('%s_labeled_copy.mat',file),'label_struct');

if strcmp(channel, 'G')
    data_2p = data.Fc1Gdff_flt2;
elseif strcmp(channel, 'R')
    data_2p = data.Fc1Rdff_flt2;
end

t_2p = (1:size(data_2p,2))/data.fs_2p; 
t_eeg = (1:size(data.EEG,2))/data.fs_eeg;
mean_dff = mean(data_2p);
ylims = [prctile(mean_dff,100*y_pct(1)), prctile(mean_dff,100*y_pct(2))];

if isfield(label_struct,sprintf('%s_%s',channel,mode))
    event_labeled = label_struct.(sprintf('%s_%s',channel,mode));
    labels = event_labeled.labels;
    cell_ind = event_labeled.cell_ind;
    cell_order = event_labeled.cell_order;
else
    cell_ind = 1;
    cell_order = label_struct.(sprintf('%s_preictal',channel)).cell_order;
    labels = cell(size(cell_order));
end

n_cells = size(data_2p,1);

figure('Position',[200 200 1000 800])
subplot(4,1,1)
mouse1 = image(i1,'xdata',[0 .15],'ydata',[0.2 0.7]);
hold on; 
ylim([0 1])
hold on; 
xlim([0 1])
image(i2,'xdata',[0.9 1],'ydata',[0 0.9])
yticks([0.5])
yticklabels({''})
xticklabels([])

m1_pos = 0;

xlims = [start_time, end_time];

subplot(4,1,2)
plot(t_2p,mean_dff)
hold on
xlabel('Time (s)')
xlim(xlims)
ylim([-.5 max(mean_dff)+1])

loop_flag = true;
while(loop_flag)
    subplot(4,1,1)
    title(sprintf('%s; cell %d; %d/%d done (%2.1f pct.); preictal; %s channel',file,cell_order(cell_ind),cell_ind,length(cell_order),100*cell_ind/length(cell_order),channel));
    progress = 100*cell_ind/length(cell_order);
    
    subplot(4,1,[3 4])
    plot(t_2p,data_2p(cell_order(cell_ind),:));
    xlim(xlims)
    ylim([-.5 2*max(mean_dff)])
    hold on
    
%     t1 = round(t_start*fs_2p);
%     t2 = round(t_end*fs_2p);
    
    temp_labels = ginput();
    
    if show_labels && ~(isempty(temp_labels))
        lb_plot = scatter(temp_labels(:,1),temp_labels(:,2),200,'k*');
    end
    
    hold off
    
    action = input(input_s);
    if isempty(action)
        action = 'next';
    end
        
    switch action
        case 'quit'
            labels{cell_order(cell_ind)} = temp_labels;
            loop_flag = false;
        case 'done'
            labels{cell_order(cell_ind)} = temp_labels;
            loop_flag = false;
        case 'redo'
            delete(lb_plot)
            temp_labels = [];
            continue
        case 'next'
            labels{cell_order(cell_ind)} = temp_labels;
            cell_ind = cell_ind + 1;
        case 'jump'
            labels{cell_order(cell_ind)} = temp_labels;
            new_t = input('Enter a cell number to jump to: ');
            cell_ind = find(cell_order == round(new_t),1);
        otherwise
            labels{cell_order(cell_ind)} = temp_labels;
            cell_ind = cell_ind + 1;
    end
    
    
    if m1_pos < 0.9
        delete(mouse1)
        
        m1_pos = progress/100;
        
        subplot(4,1,1)
        hold on
        mouse1 = image(i1,'xdata',[0 .15]+m1_pos,'ydata',[0.2 0.7]);
        
    end
end


%% Save indiv. cell labels for big event
%Run this section to save the individual trace event output from above. 
%Note that if you do not run this section after quiting the above GUI that 
%data will not be saved.

for ii=1:cell_ind-1
    if isempty(labels{cell_order(ii)})
        labels{cell_order(ii)}=NaN;
    end
end

event_labeled = struct();
event_labeled.labels = labels;
event_labeled.start_time = start_time;
event_labeled.end_time = end_time;
event_labeled.cell_ind = cell_ind;
event_labeled.cell_order = cell_order;

label_struct.(sprintf('%s_%s',channel,mode)) = event_labeled;

save(sprintf('%s_labeled.mat',file),'label_struct');
disp(sprintf('saved %s indiv to struct',mode))

