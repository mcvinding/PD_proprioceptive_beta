%%%%% GROUP SUMMARY TFR FOR PD-PROJ: BETA PART 
addpath /home/mikkel/PD_motor/global_scripts
[dirs, sub_info, lh_subs] = PD_proj_setup('rebound');
dirs.output = '/home/mikkel/PD_motor/rebound/groupAnalysis';

cd(dirs.megDir);

subs = dir(dirs.megDir);                                %Find subjects in folder
subs = {subs([subs.isdir]).name};                       %Make list
subs = subs(~(strcmp('.',subs)|strcmp('..',subs)));     %Remove dots

badsubs = {'0393'};  % Too few trials/bad data
subs(cellfun(@(x) any(strcmp(x,badsubs)),subs)) = [];

PD_subs = intersect(sub_info.PD,subs);
ctrl_subs = intersect(sub_info.ctrl,subs);

%% Define bands
theta_band      = [4 7];
betamu_band     = [8 30];
lowgamma_band   = [35 45];
highgamma_band  = [55 100];
gamma_band      = [35 100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Log/transform and extract single channel for Beta & Mu analysis

% Find peak
time = [0.05 .110];

data_PD1 = cell(1,length(PD_subs));
data_PD2 = cell(1,length(PD_subs));
data_ctrl1 = cell(1,length(ctrl_subs));
data_ctrl2 = cell(1,length(ctrl_subs));

PD_beta1 = cell(1,length(PD_subs));
PD_beta2 = cell(1,length(PD_subs));
ctrl_beta1 = cell(1,length(ctrl_subs));
ctrl_beta2 = cell(1,length(ctrl_subs));

data_PD1log = PD_beta1;
data_PD2log = PD_beta2;
data_ctrl1log = ctrl_beta1;
data_ctrl2log = ctrl_beta2;

PD_beta1bs = PD_beta1;
PD_beta2bs = PD_beta2;
ctrl_beta1bs = ctrl_beta1;
ctrl_beta2bs = ctrl_beta2;

channam_PD = cell(1,length(PD_subs));
channam_ctrl = cell(1,length(ctrl_subs));


for ii = 1:length(PD_subs)
    subID = PD_subs{ii};
    disp(subID);
    sub_dir = [dirs.megDir,'/',subID];
    load([sub_dir,'/',subID,'_1-evoked.mat']);
    timelocked = ft_combineplanar([],timelocked);
    
    [channam_PD{ii},~,~] = find_ERF_peak(timelocked,time);
    
    % Session 1
    load([sub_dir,'/',subID,'_1-tfr.mat']);
    
    data_PD1{ii} = tfr;
    
    cfg           = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'log10';
    
    data_PD1log{ii} = ft_math(cfg,data_PD1{ii});
    
    % Select channel
    cfg = [];
    cfg.channel = channam_PD{ii};
    cfg.frequency = betamu_band;
    
    PD_beta1{ii} = ft_selectdata(cfg, data_PD1log{ii});
    
    PD_beta1{ii}.label = {'peak_channel'};
    
    cfg = [];
    cfg.baseline        = [-inf -0.2];
    cfg.baselinetype    = 'absolute';
    cfg.parameter       = 'powspctrm';
    
    PD_beta1bs{ii} = ft_freqbaseline(cfg, PD_beta1{ii});
        
    clear tfr
    
    % Session 2
    load([sub_dir,'/',subID,'_2-tfr.mat']);
    
    data_PD2{ii} = tfr;
      
    cfg           = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'log10';
    
    data_PD2log{ii} = ft_math(cfg, data_PD2{ii});
       
    cfg = [];
    cfg.channel = channam_PD{ii};
    cfg.frequency = betamu_band;
    
    PD_beta2{ii} = ft_selectdata(cfg, data_PD2log{ii});
    
    PD_beta2{ii}.label = {'peak_channel'};

    cfg = [];
    cfg.baseline        = [-inf -0.2];
    cfg.baselinetype    = 'absolute';
    cfg.parameter       = 'powspctrm';
    
    PD_beta2bs{ii} = ft_freqbaseline(cfg,  PD_beta2{ii});
        
    clear tfr
end

for ii = 1:length(ctrl_subs)
    subID = ctrl_subs{ii};
    disp(subID);
    sub_dir = [dirs.megDir,'/',subID];
    load([sub_dir,'/',subID,'_1-evoked.mat']);
    timelocked = ft_combineplanar([],timelocked);
    
    [channam_ctrl{ii},~,~] = find_ERF_peak(timelocked,time);
    
    % Session 1
    load([sub_dir,'/',subID,'_1-tfr.mat']);
      
    data_ctrl1{ii} = tfr;

    cfg           = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'log10';
    
    data_ctrl1log{ii} = ft_math(cfg,data_ctrl1{ii});
    
    cfg = [];
    cfg.channel = channam_ctrl{ii};
    cfg.frequency = betamu_band;
%     cfg.latency = [-.5 2.5];
    
    ctrl_beta1{ii} = ft_selectdata(cfg, data_ctrl1log{ii});
    
    ctrl_beta1{ii}.label = {'peak_channel'};
    
    cfg = [];
    cfg.baseline        = [-inf -0.2];
    cfg.baselinetype    = 'absolute';
    cfg.parameter       = 'powspctrm';
        
    ctrl_beta1bs{ii} = ft_freqbaseline(cfg, ctrl_beta1{ii});

    % Session 2
    load([sub_dir,'/',subID,'_2-tfr.mat']);
    
    data_ctrl2{ii} = tfr;

    cfg           = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'log10';
    
    data_ctrl2log{ii} = ft_math(cfg, data_ctrl2{ii});
    
    cfg = [];
    cfg.channel = channam_ctrl{ii};
    cfg.frequency = betamu_band;
%     cfg.latency = [-.5 2.5];
    
    ctrl_beta2{ii} = ft_selectdata(cfg, data_ctrl2log{ii});
    
    ctrl_beta2{ii}.label = {'peak_channel'};
    
    cfg = [];
    cfg.baseline        = [-inf -0.2];
    cfg.baselinetype    = 'absolute';
    cfg.parameter       = 'powspctrm';
        
    ctrl_beta2bs{ii} = ft_freqbaseline(cfg, ctrl_beta2{ii});
    
end

% ctrl_betaAvg = ft_freqgrandaverage([],ctrl_beta1bs{:});
% PD_betaAvg = ft_freqgrandaverage([],PD_beta1bs{:});


save([dirs.output,'/betaChan1log.mat'],'ctrl_beta1','PD_beta1','ctrl_beta2','PD_beta2','-v7.3');
disp('saved 1 of 5');
save([dirs.output,'/betaChan1bs.mat'],'ctrl_beta1bs','PD_beta1bs','ctrl_beta2bs','PD_beta2bs','-v7.3');
disp('saved 2 of 5');
save([dirs.output,'/all_TFRlog.mat'],'data_ctrl1log','data_ctrl2log','data_PD1log','data_PD2log','-v7.3');
disp('saved 3 of 5');
save([dirs.output,'/all_TFR.mat'],'data_ctrl1','data_ctrl2','data_PD1','data_PD2','-v7.3');
disp('saved 4 of 5');
save([dirs.output,'/peak_chans.mat'],'channam_ctrl','channam_PD','-v7.3');
disp('done')

cd(dirs.output);
clear data* ctrl_beta* PD_beta*

%% Low-freq: Flip if necessary and baseline correct
load([dirs.output,'/all_TFRlog.mat']);
disp('done')

%%
data_PD1bs = cell(1,length(PD_subs));
data_PD2bs = cell(1,length(PD_subs));
data_ctrl1bs = cell(1,length(ctrl_subs));
data_ctrl2bs = cell(1,length(ctrl_subs));

%Baseline options
cfg = [];
cfg.baseline        = [-inf -0.2];
cfg.baselinetype    = 'absolute';
cfg.parameter       = 'powspctrm';

for ii = 1:length(PD_subs)
    subID = PD_subs{ii};
    disp(subID);
%     sub_dir = [dirs.megDir,'/',subID];
    
    %Session 1+2
    if any(strcmp(lh_subs,subID))
        disp('Flip')
        data_PD1log{ii} = flip_sens_neuromag(data_PD1log{ii});
        data_PD2log{ii} = flip_sens_neuromag(data_PD2log{ii});
    end    
    data_PD1bs{ii} = ft_freqbaseline(cfg,data_PD1log{ii});
    data_PD2bs{ii} = ft_freqbaseline(cfg,data_PD2log{ii});
    
end

for ii = 1:length(ctrl_subs)
    subID = ctrl_subs{ii};
    disp(subID);
    
    %Session 1+2
    if any(strcmp(lh_subs,subID))
        disp('Flip')
        data_ctrl1log{ii} = flip_sens_neuromag(data_crtl1log{ii});
        data_ctrl2log{ii} = flip_sens_neuromag(data_ctrl2log{ii});
    end    
    data_ctrl1bs{ii} = ft_freqbaseline(cfg,data_ctrl1log{ii});
    data_ctrl2bs{ii} = ft_freqbaseline(cfg,data_ctrl2log{ii});
    
end

save([dirs.output,'/all_TFRbs.mat'],'data_ctrl1bs','data_ctrl2bs','data_PD1bs','data_PD2bs','-v7.3');
disp('done')

%% Average
% cfg = [];
% cfg.keepindividual = 'yes';
% 
% avgTFR.PD1 = ft_freqgrandaverage(cfg, data_PD1{:});
% avgTFR.PD2 = ft_freqgrandaverage(cfg, data_PD2{:});
% avgTFR.ctrl1 = ft_freqgrandaverage(cfg, data_ctrl1{:});
% avgTFR.ctrl2 = ft_freqgrandaverage(cfg, data_ctrl2{:});
% 
% % Save data
% save([dirs.output,'/avgTFR.mat'],'avgTFR','-v7.3');
% disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% High-freq: Combine into one structure, flip if necessary and baseline correct
data_PD1_hf = cell(1,length(PD_subs));
data_PD2_hf = cell(1,length(PD_subs));
data_ctrl1_hf = cell(1,length(ctrl_subs));
data_ctrl2_hf = cell(1,length(ctrl_subs));

data_PD1_hf_log = data_PD1_hf;
data_PD2_hf_log = data_PD2_hf;
data_ctrl1_hf_log = data_ctrl1_hf;
data_ctrl2_hf_log = data_ctrl1_hf;

PD_hfsingle1 = data_PD1_hf;
PD_hfsingle2 = data_PD2_hf;
ctrl_hfsingle1 = data_ctrl1_hf;
ctrl_hfsingle2 = data_ctrl2_hf;

PD_hfsingle1bs = data_PD1_hf;
PD_hfsingle2bs = data_PD2_hf;
ctrl_hfsingle1bs = data_ctrl1_hf;
ctrl_hfsingle2bs = data_ctrl2_hf;

data_PD1_hf_bs = data_PD1_hf;
data_PD2_hf_bs = data_PD2_hf;
data_ctrl1_hf_bs = data_ctrl1_hf;
data_ctrl2_hf_bs = data_ctrl2_hf;

for ii = 1:length(PD_subs)
    subID = PD_subs{ii};
    disp(subID);
    sub_dir = [dirs.megDir,'/',subID];
    
    % Session 1
    load([sub_dir,'/',subID,'_1-hfreq-tfr.mat']);
    
    data_PD1_hf{ii} = tfr_hf;
    
    cfg           = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'log10';
    
    data_PD1_hf_log{ii} = ft_math(cfg,data_PD1_hf{ii});
    
    % Select channel
    cfg = [];
    cfg.channel = channam_PD{ii};
%     cfg.frequency = betamu_band;
    
    PD_hfsingle1{ii} = ft_selectdata(cfg, data_PD1_hf_log{ii});
    
    PD_hfsingle1{ii}.label = {'peak_channel'};
    
    cfg = [];
    cfg.baseline        = [-inf -0.05];
    cfg.baselinetype    = 'absolute';
    cfg.parameter       = 'powspctrm';
    
    PD_hfsingle1bs{ii} = ft_freqbaseline(cfg, PD_hfsingle1{ii});
       
    if any(strcmp(lh_subs,subID))
        disp('Flip')
        data_PD1_hf_log{ii} = flip_sens_neuromag(data_PD1_hf_log{ii});
    end
    
%     cfg = [];
%     cfg.baseline        = [-inf -0.05];
%     cfg.baselinetype    = 'absolute';
%     cfg.parameter       = 'powspctrm';    
    data_PD1_hf_bs{ii} = ft_freqbaseline(cfg,data_PD1_hf_log{ii});
    
    % Session 2
    load([sub_dir,'/',subID,'_2-hfreq-tfr.mat']);
    
    data_PD2_hf{ii} = tfr_hf;
    
    cfg           = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'log10';
    
    data_PD2_hf_log{ii} = ft_math(cfg,data_PD2_hf{ii});
    
    % Select channel
    cfg = [];
    cfg.channel = channam_PD{ii};
%     cfg.frequency = betamu_band;
    
    PD_hfsingle2{ii} = ft_selectdata(cfg, data_PD2_hf_log{ii});
    
    PD_hfsingle2{ii}.label = {'peak_channel'};
    
    cfg = [];
    cfg.baseline        = [-inf -0.05];
    cfg.baselinetype    = 'absolute';
    cfg.parameter       = 'powspctrm';
    
    PD_hfsingle2bs{ii} = ft_freqbaseline(cfg, PD_hfsingle2{ii});
       
    if any(strcmp(lh_subs,subID))
        disp('Flip')
        data_PD2_hf_log{ii} = flip_sens_neuromag(data_PD2_hf_log{ii});
    end
    
%     cfg = [];
%     cfg.baseline        = [-inf -0.05];
%     cfg.baselinetype    = 'absolute';
%     cfg.parameter       = 'powspctrm';    
    data_PD2_hf_bs{ii} = ft_freqbaseline(cfg,data_PD2_hf_log{ii});
end

for ii = 1:length(ctrl_subs)
    subID = ctrl_subs{ii};
    disp(subID);
    sub_dir = [dirs.megDir,'/',subID];
    
    % Session 1
    load([sub_dir,'/',subID,'_1-hfreq-tfr.mat']);
    
    data_ctrl1_hf{ii} = tfr_hf;
    
    cfg           = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'log10';
    
    data_ctrl1_hf_log{ii} = ft_math(cfg,data_ctrl1_hf{ii});
    
    % Select channel
    cfg = [];
    cfg.channel = channam_ctrl{ii};
    
    ctrl_hfsingle1{ii} = ft_selectdata(cfg, data_ctrl1_hf_log{ii});
    
    ctrl_hfsingle1{ii}.label = {'peak_channel'};
    
    cfg = [];
    cfg.baseline        = [-inf -0.05];
    cfg.baselinetype    = 'absolute';
    cfg.parameter       = 'powspctrm';
    
    ctrl_hfsingle1bs{ii} = ft_freqbaseline(cfg, ctrl_hfsingle1{ii});
             
    data_ctrl1_hf_bs{ii} = ft_freqbaseline(cfg,data_ctrl1_hf_log{ii});
    
    % Session 2
    load([sub_dir,'/',subID,'_2-hfreq-tfr.mat']);
    
    data_ctrl2_hf{ii} = tfr_hf;
    
    cfg           = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'log10';
    
    data_ctrl2_hf_log{ii} = ft_math(cfg,data_ctrl2_hf{ii});
    
    % Select channel
    cfg = [];
    cfg.channel = channam_ctrl{ii};
%     cfg.frequency = betamu_band;
    
    ctrl_hfsingle2{ii} = ft_selectdata(cfg, data_ctrl2_hf_log{ii});
    
    ctrl_hfsingle2{ii}.label = {'peak_channel'};
    
    cfg = [];
    cfg.baseline        = [-inf -0.05];
    cfg.baselinetype    = 'absolute';
    cfg.parameter       = 'powspctrm';
    
    ctrl_hfsingle2bs{ii} = ft_freqbaseline(cfg, ctrl_hfsingle2{ii});
      
    data_ctrl2_hf_bs{ii} = ft_freqbaseline(cfg, data_ctrl2_hf_log{ii});
end

save([dirs.output,'/all_hfreq_data.mat'],'data_ctrl1_hf','data_ctrl2_hf','data_PD1_hf','data_PD2_hf','-v7.3');
disp('saved 1 of 5')
save([dirs.output,'/all_hfreq_log.mat'],'data_ctrl1_hf_log','data_ctrl2_hf_log','data_PD1_hf_log','data_PD2_hf_log','-v7.3');
disp('saved 2 of 5')
save([dirs.output,'/all_hfreq_bs.mat'],'data_ctrl1_hf_bs','data_ctrl2_hf_bs','data_PD1_hf_bs','data_PD2_hf_bs','-v7.3');
disp('saved 3 of 5')
save([dirs.output,'/singlechan_hfreq_log.mat'],'ctrl_hfsingle1','ctrl_hfsingle2','PD_hfsingle2','PD_hfsingle2','-v7.3');
disp('saved 4 of 5')
save([dirs.output,'/singlechan_hfreq_bs.mat'],'ctrl_hfsingle1bs','ctrl_hfsingle2bs','PD_hfsingle2bs','PD_hfsingle2bs','-v7.3');

disp('done');

%% Average
% cfg = [];
% cfg.keepindividual = 'yes';
% 
% avgTFR_hfreq.PD1 = ft_freqgrandaverage(cfg, data_PD1_hf_bs{:});
% avgTFR_hfreq.PD2 = ft_freqgrandaverage(cfg, data_PD2_hf_bs{:});
% avgTFR_hfreq.ctrl1 = ft_freqgrandaverage(cfg, data_ctrl1_hf_bs{:});
% avgTFR_hfreq.ctrl2 = ft_freqgrandaverage(cfg, data_ctrl2_hf_bs{:});
% 
% save([dirs.output,'/avgTFR_hfreq.mat'],'avgTFR_hfreq','-v7.3');
% disp('done');



exit