%%%%% Get EMG, process and statistics
addpath /home/mikkel/PD_motor/global_scripts
[dirs, sub_info, lh_subs] = PD_proj_setup('rebound');

cd(dirs.megDir);

subs = dir(dirs.megDir);                                %Find subjects in folder
subs = {subs([subs.isdir]).name};                       %Make list
subs = subs(~(strcmp('.',subs)|strcmp('..',subs)));     %Remove dots

% Remove bad subject
badsubs = {'0393'}; % Too few trials/bad data
subs(strcmp(badsubs,subs)) = [];

PD_subs = intersect(sub_info.PD,subs);
ctrl_subs = intersect(sub_info.ctrl,subs);

%% Load data (get trials from ERF files, extract EMG from raw)
PD_emg1 = cell(1,length(PD_subs));
PD_emg2 = cell(1,length(PD_subs));
ctrl_emg1 = cell(1,length(ctrl_subs));
ctrl_emg2 = cell(1,length(ctrl_subs));

PD_emg1_freq = cell(1,length(PD_subs));
PD_emg2_freq = cell(1,length(PD_subs));
ctrl_emg1_freq = cell(1,length(ctrl_subs));
ctrl_emg2_freq = cell(1,length(ctrl_subs));

for ii = 1:length(subs)
    subID = subs{ii};
    disp(subID);
    sub_dir = [dirs.megDir,'/',subID];
    
    %%% SESSION 1 %%%

    data_file = [sub_dir,'/',subID,'_1-ica_raw.fif'];
    load([sub_dir,'/',subID,'_1-epochs.mat']);
    disp('loaded epochs file 1');
    
    % Get trial structure
    trl_def = [data.sampleinfo,-1500*ones(length(data.trialinfo),1),data.trialinfo];
    
    clear data
    
    disp(['Reading data from ',data_file]);
    
    cfg = [];
    if any(strcmp(lh_subs,subID))
        cfg.channel     = 'EMG005';
    else
        cfg.channel     = 'EMG004';
    end 
    
    cfg.dataset     = data_file;
    cfg.trl         = trl_def;
    cfg.hpfilter    = 'no';
    cfg.lpfilter    = 'no';
    cfg.dftfilter   = 'yes';
    cfg.hpfreq      = 1;
%     cfg.lpfreq      = 400;
    cfg.dftfreq     = [50 100 150];
    cfg.demean      = 'no';         
    cfg.rectify     = 'no';

    emg_data_raw = ft_preprocessing(cfg);
    
    % adjust timeline and baselinecorrect
    cfg = [];
    cfg.offset = -25;

    emg_data_raw = ft_redefinetrial(cfg, emg_data_raw);
    
    cfg = [];
    cfg.latency = [-1.250 2.500];
    
    emg_data_slct = ft_selectdata(cfg,emg_data_raw);
    
    cfg = [];
    cfg.rectify     = 'yes';    
    emg_data_rect = ft_preprocessing(cfg,emg_data_slct);
    
    cfg = [];
    cfg.demean      = 'yes';    
    cfg.baselinewindow = [-inf 0];
    emg_data_rect = ft_preprocessing(cfg,emg_data_rect);
    
    if any(strcmp(subID,PD_subs))
        idx = find(strcmp(subID,PD_subs));
        PD_emg1{idx} = ft_timelockanalysis([],emg_data_rect);
        PD_emg1{idx}.label = {'EMG004'};
    elseif any(strcmp(subID,ctrl_subs))
        idx = find(strcmp(subID,ctrl_subs));
        ctrl_emg1{idx} = ft_timelockanalysis([],emg_data_rect);
    end
           
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'hanning';
    cfg.output      = 'pow';
    cfg.foilim     = [1 250];
    cfg.keeptrials = 'no';
    cfg.channel    = 'all';
    cfg.pad        = 'nextpow2'; 
    
    if any(strcmp(subID,PD_subs))
        idx = find(strcmp(subID,PD_subs));
        PD_emg1_freq{idx} = ft_freqanalysis(cfg,emg_data_slct);
        PD_emg1_freq{idx}.label = {'EMG004'};

    elseif any(strcmp(subID,ctrl_subs))
        idx = find(strcmp(subID,ctrl_subs));
        ctrl_emg1_freq{idx} = ft_freqanalysis(cfg,emg_data_slct);
    end
    
    clear data* emg*
      
    %%% SESSION 2 %%%
    data_file = [sub_dir,'/',subID,'_2-ica_raw.fif'];
    load([sub_dir,'/',subID,'_2-epochs.mat']);
    disp('loaded epochs file 2');
    
    % Get trial structure
    trl_def = [data.sampleinfo,-1500*ones(length(data.trialinfo),1),data.trialinfo];
    
    clear data
    
    disp(['Reading data from ',data_file]);
    
    cfg = [];
    if any(strcmp(lh_subs,subID))
        cfg.channel     = 'EMG005';
    else
        cfg.channel     = 'EMG004';
    end 
    
    cfg.dataset     = data_file;
    cfg.trl         = trl_def;
    cfg.hpfilter    = 'no';
    cfg.lpfilter    = 'no';
    cfg.dftfilter   = 'yes';
    cfg.hpfreq      = 1;
%     cfg.lpfreq      = 400;
    cfg.dftfreq     = [50 100 150];
    cfg.demean      = 'no';         
    cfg.rectify     = 'no';

    emg_data_raw = ft_preprocessing(cfg);
    
    % adjust timeline and baselinecorrect
    cfg = [];
    cfg.offset = -25;

    emg_data_raw = ft_redefinetrial(cfg, emg_data_raw);
    
    cfg = [];
    cfg.latency = [-1.250 2.500];
    
    emg_data_slct = ft_selectdata(cfg,emg_data_raw);
    
    cfg = [];
    cfg.rectify     = 'yes';    
    emg_data_rect = ft_preprocessing(cfg,emg_data_slct);
    cfg = [];
    cfg.demean      = 'yes';    
    cfg.baselinewindow = [-inf 0];
    emg_data_rect = ft_preprocessing(cfg,emg_data_rect);
    
    if any(strcmp(subID,PD_subs))
        idx = find(strcmp(subID,PD_subs));
        PD_emg2{idx} = ft_timelockanalysis([],emg_data_rect);
        PD_emg2{idx}.label = {'EMG004'};

    elseif any(strcmp(subID,ctrl_subs))
        idx = find(strcmp(subID,ctrl_subs));
        ctrl_emg2{idx} = ft_timelockanalysis([],emg_data_rect);
    end
    
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'hanning';
    cfg.output      = 'pow';
    cfg.foilim     = [1 250];
%     cfg.tapsmofrq  = .5;
    cfg.keeptrials = 'no';
    cfg.channel    = 'all';
    cfg.pad        = 'nextpow2'; 
    
    if any(strcmp(subID,PD_subs))
        idx = find(strcmp(subID,PD_subs));
        PD_emg2_freq{idx} = ft_freqanalysis(cfg,emg_data_slct);
        PD_emg2_freq{idx}.label = {'EMG004'};

    elseif any(strcmp(subID,ctrl_subs))
        idx = find(strcmp(subID,ctrl_subs));
        ctrl_emg2_freq{idx} = ft_freqanalysis(cfg,emg_data_slct);
    end
        
end
disp('Done processing. Saving...');
save([dirs.output,'/emg_raw.mat'],'PD_emg1','PD_emg2','ctrl_emg1','ctrl_emg2','-v7.3');
save([dirs.output,'/emg_freq.mat'],'PD_emg1_freq','PD_emg2_freq','ctrl_emg1_freq','ctrl_emg2_freq','-v7.3');
disp('DONE');

%% (re-)load
load([dirs.output,'/emg_raw.mat']);
load([dirs.output,'/emg_freq.mat']);
disp('done')

%% Remove sub 0376 who have muscle movements
% ctrl_emg1e = ctrl_emg1;
% ctrl_emg2e = ctrl_emg2;
% ctrl_emg1e{12} = [];
% ctrl_emg2e{12} = [];
% ctrl_emg1e = ctrl_emg1e(~cellfun('isempty',ctrl_emg1e));
% ctrl_emg2e = ctrl_emg2e(~cellfun('isempty',ctrl_emg2e));
% 
% ctrl_emg1_freqE = ctrl_emg1_freq;
% ctrl_emg2_freqE = ctrl_emg2_freq;
% ctrl_emg1_freqE{12} = [];
% ctrl_emg2_freqE{12} = [];
% ctrl_emg1_freqE = ctrl_emg1_freqE(~cellfun('isempty',ctrl_emg1_freqE));
% ctrl_emg2_freqE = ctrl_emg2_freqE(~cellfun('isempty',ctrl_emg2_freqE));

%% Grand average
PD1_avg = ft_timelockgrandaverage([],PD_emg1{:});
ctrl1_avg = ft_timelockgrandaverage([],ctrl_emg1{:});
PD2_avg = ft_timelockgrandaverage([],PD_emg2{:});
ctrl2_avg = ft_timelockgrandaverage([],ctrl_emg2{:});

PD1_avgFreq = ft_freqgrandaverage([],PD_emg1_freq{:});
ctrl1_avgFreq = ft_freqgrandaverage([],ctrl_emg1_freq{:});
PD2_avgFreq = ft_freqgrandaverage([],PD_emg2_freq{:});
ctrl2_avgFreq = ft_freqgrandaverage([],ctrl_emg2_freq{:});

figure; hold on
plot(PD1_avg.time,PD1_avg.avg,'b')
plot(PD2_avg.time,PD2_avg.avg,'b--')
plot(ctrl1_avg.time,ctrl1_avg.avg,'r')
plot(ctrl2_avg.time,ctrl2_avg.avg,'r--')
axis([-1 2.5 -1e-3 1e-3])

figure; hold on
plot(PD1_avgFreq.freq,PD1_avgFreq.powspctrm,'b')
plot(PD2_avgFreq.freq,PD2_avgFreq.powspctrm,'b--')
plot(ctrl1_avgFreq.freq,ctrl1_avgFreq.powspctrm,'r')
plot(ctrl2_avgFreq.freq,ctrl2_avgFreq.powspctrm,'r--')

figure;
for i=1:length(ctrl_emg1)
%     if i == 10; continue; end;
    subplot(length(ctrl_emg1)/2,2,i); hold on
    plot(ctrl_emg1{i}.time,ctrl_emg1{i}.avg);
    plot(ctrl_emg2{i}.time,ctrl_emg2{i}.avg,'r');
    title(ctrl_subs(i))
end; hold off

figure;
for i=1:length(PD_emg1)
%     if i == 10; continue; end;
    subplot(round(length(PD_emg1)/2),2,i); hold on
    plot(PD_emg1{i}.time,PD_emg1{i}.avg);
    plot(PD_emg2{i}.time,PD_emg2{i}.avg,'r');
    title(PD_subs(i))
end; hold off

figure;
for i=1:length(ctrl_emg1_freq)
    plot(ctrl_emg1_freq{i}.freq,ctrl_emg1_freq{i}.powspctrm); hold on
end; hold off

%% Statistics
% PD1 vs Ctrl1
cfg = [];
cfg.method              = 'montecarlo';
cfg.neighbours          = [];
cfg.channel             = 'EMG004';

cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';

cfg.design = [ones(1,length(PD_emg1_freq)) 2*ones(1,length(ctrl_emg1_freq))];
cfg.design = [cfg.design; 1:length(cfg.design)];

emg_stat_1v1 = ft_freqstatistics(cfg, PD_emg1_freq{:}, ctrl_emg1_freq{:});
emg_stat_1v1raw = ft_timelockstatistics(cfg, PD_emg1{:}, ctrl_emg1{:});
disp('done');

% PD2 vs Ctrl2
cfg.design = [ones(1,length(PD_emg2_freq)) 2*ones(1,length(ctrl_emg2_freq))];
cfg.design = [cfg.design; 1:length(cfg.design)];

emg_stat_2v2 = ft_freqstatistics(cfg, PD_emg2_freq{:}, ctrl_emg2_freq{:});
emg_stat_2v2raw = ft_timelockstatistics(cfg, PD_emg2{:}, ctrl_emg2{:});

disp('done');

% PD1 vs PD2
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';

cfg.design  = [ones(1,length(PD_emg1_freq)) 2*ones(1,length(PD_emg2_freq))];
cfg.design  = [cfg.design; [1:length(PD_emg1_freq),1:length(PD_emg2_freq)]];
cfg.ivar    = 1;
cfg.uvar    = 2;

emg_stat_1v2PD = ft_freqstatistics(cfg, PD_emg1_freq{:}, PD_emg2_freq{:});
emg_stat_1v2PDraw = ft_timelockstatistics(cfg, PD_emg1{:}, PD_emg2{:});
disp('done');

% Ctrl1 vs Ctrl2
cfg.design  = [ones(1,length(ctrl_emg1_freq)) 2*ones(1,length(ctrl_emg2_freq(:)))];
cfg.design  = [cfg.design; [1:length(ctrl_emg1_freq),1:length(ctrl_emg2_freq(:))]];

emg_stat_1v2ctrl = ft_freqstatistics(cfg, ctrl_emg1_freq{:}, ctrl_emg2_freq{:});
emg_stat_1v2ctrlraw = ft_timelockstatistics(cfg, ctrl_emg1{:}, ctrl_emg2{:});
disp('done');

%% Save
save(fullfile(dirs.output,'EMG_stats.mat'), ...
    'emg_stat_1v2ctrl','emg_stat_1v2ctrlraw', ...
    'emg_stat_1v2PD', 'emg_stat_1v2PDraw', ...
    'emg_stat_2v2','emg_stat_2v2raw', ...
    'emg_stat_1v1','emg_stat_1v1raw','-v7.3');
disp('Done')
    





















