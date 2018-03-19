%%%%% STATISTICS %%%%%%%
addpath /home/mikkel/PD_motor/global_scripts
[dirs, sub_info, lh_subs] = PD_proj_setup('rebound');

cd(dirs.megDir);

subs = dir(dirs.megDir);                                %Find subjects in folder
subs = {subs([subs.isdir]).name};                       %Make list
subs = subs(~(strcmp('.',subs)|strcmp('..',subs)));     %Remove dots

badsubs = {'0393'}; % Too few trials/bad data
subs(cellfun(@(x) any(strcmp(x,badsubs)),subs)) = [];

PD_subs = intersect(sub_info.PD,subs);
ctrl_subs = intersect(sub_info.ctrl,subs);

%% Define bands
theta_band      = [4 7];
betamu_band     = [8 30];
lowgamma_band   = [31 49];
highgamma_band  = [51 100];
gamma_band      = [31 100];

%% ########################################################################
%       Beta/Mu statistics
% #########################################################################
%% Load
cd(dirs.output);
load([dirs.output,'/betaChan1bs.mat']);
disp('data loaded');

%% Stats
% 1st session
cfg = [];
cfg.method              = 'montecarlo';
cfg.neighbours          = [];
cfg.channel             = 'peak_channel';

cfg.design = [ones(1,length(PD_beta1bs)) 2*ones(1,length(ctrl_beta1bs))];
cfg.design = [cfg.design; 1:length(cfg.design)];

cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic = 'maxsum';

cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

stat_beta1 = ft_freqstatistics(cfg, PD_beta1bs{:}, ctrl_beta1bs{:});
disp('done');

% 2nd session
cfg = [];
cfg.method              = 'montecarlo';
cfg.neighbours          = [];
cfg.channel             = 'peak_channel';

cfg.design = [ones(1,length(PD_beta2bs)) 2*ones(1,length(ctrl_beta2bs))];
cfg.design = [cfg.design; 1:length(cfg.design)];

cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic = 'maxsum';

cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

stat_beta2 = ft_freqstatistics(cfg, PD_beta2bs{:}, ctrl_beta2bs{:});
disp('done');

%% Interaction with session in beta/mu band
%Get difference
load([dirs.output,'/betaChan1log.mat']);
disp('loaded data')

PD_betaDiff = cell(1,length(PD_subs));
ctrl_betaDiff = cell(1,length(ctrl_subs));
PD_betaDiffbs = PD_betaDiff;
ctrl_betaDiffbs = ctrl_betaDiff;

cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = '(x1-x2)';
% cfg.operation = 'x1-x2';
for ii = 1:length(PD_beta1)
    PD_betaDiff{ii} = ft_math(cfg, PD_beta2{ii},PD_beta1{ii});
end

for ii = 1:length(ctrl_subs)
    ctrl_betaDiff{ii} = ft_math(cfg,ctrl_beta2{ii},ctrl_beta1{ii});
end

GA_diff1 = ft_freqgrandaverage([],ctrl_betaDiff{:});
GA_diff2 = ft_freqgrandaverage([],PD_betaDiff{:});

cfg                 = [];
cfg.baseline        = [-inf -0.2];
cfg.baselinetype    = 'absolute';
cfg.parameter       = 'powspctrm';

for ii = 1:length(PD_beta1)
    PD_betaDiffbs{ii} = ft_freqbaseline(cfg, PD_betaDiff{ii});
end

for ii = 1:length(ctrl_subs)
    ctrl_betaDiffbs{ii} = ft_freqbaseline(cfg,ctrl_betaDiff{ii});
end

GA_diff1bs = ft_freqgrandaverage([],ctrl_betaDiffbs{:});
GA_diff2bs = ft_freqgrandaverage([],PD_betaDiffbs{:});

% Run stats
cfg = [];
cfg.method              = 'montecarlo';
cfg.neighbours          = [];
cfg.channel             = 'peak_channel';

cfg.design = [ones(1,length(PD_betaDiff)) 2*ones(1,length(ctrl_betaDiff))];
cfg.design = [cfg.design; 1:length(cfg.design)];

cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic = 'maxsum';

cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

stat_betaIntraction = ft_freqstatistics(cfg, PD_betaDiffbs{:}, ctrl_betaDiffbs{:});

%% Correlation with UPDRS score (PD only)

% Scores from behavious assesment (manually entered here)
scores1 = [38 25 31 10 29 22 27 29 28 21 41 61];
scores2 = [21  7 16  5 19  6  9 14 15 12 31 39];
subs = [1:12, 1:12];
reg_design = [scores1,scores2; subs];

% Difference
score_diff = scores1-scores2;
subs = [1:12];
regDif_design = [score_diff;subs];

% Difference based UPDRS
cfg = [];
cfg.method              = 'montecarlo';
cfg.neighbours          = [];
cfg.channel             = 'peak_channel';

cfg.design              =  reg_design;
cfg.statistic           = 'ft_statfun_depsamplesregrT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';

cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2;
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

stats_UPDRS = ft_freqstatistics(cfg, PD_beta1bs{:}, PD_beta2bs{:});

%% NEW 2018-01-10: Collapse across session (due to n.s. interaction)
% 1st session
PD_pool = cell(1,length(PD_beta1bs));
for i = 1:length(PD_beta1bs)
    PD_pool{i} = ft_freqgrandaverage([],PD_beta1bs{i},PD_beta2bs{i});
end
ctrl_pool = cell(1,length(ctrl_beta1bs));
for i = 1:length(ctrl_beta1bs)
    ctrl_pool{i} = ft_freqgrandaverage([],ctrl_beta1bs{i},ctrl_beta2bs{i});
end

cfg = [];
cfg.method              = 'montecarlo';
cfg.neighbours          = [];
cfg.channel             = 'peak_channel';

design = [ones(1,length(PD_beta1bs)) 2*ones(1,length(ctrl_beta1bs))];
design = [design; 1:length(design)];
cfg.design = design;

cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic = 'maxsum';

cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

stat_grpMain = ft_freqstatistics(cfg, PD_pool{:}, ctrl_pool{:});
disp('done');

%% Save
save([dirs.output,'/beta_stats.mat'],'stat_betaIntraction','stat_beta1','stat_beta2','stats_UPDRS','stat_grpMain','-v7.3');
disp('done');

clear ctrl_beta* PD_beta* stat*

%% ########################################################################
% Gamma
% #########################################################################
%% Load data
cd(dirs.output);
% load('singlechan_hfreq_bs.mat')
load('all_hfreq_bs')
disp('done');

%% Prepare neighbours structure
cfg = [];
cfg.layout  = 'neuromag306cmb.lay';
cfg.method  = 'distance';
neighbours = ft_prepare_neighbours(cfg, data_PD1_hf_bs{1});


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% High gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Difference between 1st sessions
cfg = [];
% cfg.latency             = [-.5 2.5];
cfg.frequency           = highgamma_band;
cfg.avgoverchan         = 'no';         % Analyse all channels
cfg.avgoverfreq         = 'yes';         % Collapse frequency bins

cfg.method              = 'montecarlo';
cfg.channel             = 'MEGGRAD';

cfg.design = [ones(1,length(data_PD1_hf_bs)) 2*ones(1,length(data_ctrl1_hf_bs))];
cfg.design = [cfg.design; 1:length(cfg.design)];

cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.neighbours          = neighbours;
cfg.clusteralpha        = 0.05;
cfg.minnbchan           = 2;
cfg.clusterstatistic    = 'maxsum';

cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

stat_hiGamma_avg = ft_freqstatistics(cfg, data_PD1_hf_bs{:}, data_ctrl1_hf_bs{:});

%% save
save([dirs.output,'/stat_hiGamma_avg.mat'],'stat_hiGamma_avg','-v7.3');
disp('done')

%% ########################################################################
% Low gamma
% #########################################################################
% Difference between 1st sessions
cfg = [];
cfg.frequency           = lowgamma_band;
cfg.avgoverchan         = 'no';             % Analyse all channels
cfg.avgoverfreq         = 'yes';            % Collapse frequency bins

cfg.method              = 'montecarlo';
cfg.channel             = 'MEGGRAD';

cfg.design = [ones(1,length(data_PD1_hf_bs)) 2*ones(1,length(data_ctrl1_hf_bs))];
cfg.design = [cfg.design; 1:length(cfg.design)];

cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.neighbours          = neighbours;
cfg.clusteralpha        = 0.05;
cfg.minnbchan           = 2;
cfg.clusterstatistic    = 'maxsum';

cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

stat_lowGamma_avg = ft_freqstatistics(cfg, data_PD1_hf_bs{:}, data_ctrl1_hf_bs{:});

%% save
save([dirs.output,'/stat_lowGamma_avg.mat'],'stat_lowGamma_avg','-v7.3');
disp('done')

%% Clear some stuff
clear stat* data*

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAMMA effect of meds. (interaction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get difference
load([dirs.output,'/all_hfreq_log.mat']);
disp('loaded data')

PD_hfDiff = cell(1,length(PD_subs));
ctrl_hfDiff = cell(1,length(ctrl_subs));
PD_hfDiffbs = PD_hfDiff;
ctrl_hfDiffbs = ctrl_hfDiff;

cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = '(x2-x1)';
for ii = 1:length(data_PD1_hf_log)
    PD_hfDiff{ii} = ft_math(cfg, data_PD1_hf_log{ii},data_PD2_hf_log{ii});
end

for ii = 1:length(data_ctrl1_hf_log)
    ctrl_hfDiff{ii} = ft_math(cfg,data_ctrl1_hf_log{ii},data_ctrl2_hf_log{ii});
end

cfg                 = [];
cfg.baseline        = [-inf -0.05];
cfg.baselinetype    = 'absolute';
cfg.parameter       = 'powspctrm';

for ii = 1:length(PD_hfDiff)
    PD_hfDiffbs{ii} = ft_freqbaseline(cfg, PD_hfDiff{ii});
end

for ii = 1:length(ctrl_hfDiff)
    ctrl_hfDiffbs{ii} = ft_freqbaseline(cfg,ctrl_hfDiff{ii});
end

GA_diff1 = ft_freqgrandaverage([],ctrl_hfDiffbs{:});
GA_diff2 = ft_freqgrandaverage([],PD_hfDiffbs{:});

disp('done')

%% Low gamma
cfg = [];
cfg.frequency           = lowgamma_band;
cfg.avgoverchan         = 'no';             % Analyse all channels
cfg.avgoverfreq         = 'yes';            % Collapse frequency bins

cfg.method              = 'montecarlo';
cfg.channel             = 'MEGGRAD';

cfg.design = [ones(1,length(PD_hfDiffbs)) 2*ones(1,length(ctrl_hfDiffbs))];
cfg.design = [cfg.design; 1:length(cfg.design)];

cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.neighbours          = neighbours;
cfg.clusteralpha        = 0.05;
cfg.minnbchan           = 2;
cfg.clusterstatistic    = 'maxsum';

cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

stat_lowGammaInteraction = ft_freqstatistics(cfg, PD_hfDiffbs{:}, ctrl_hfDiffbs{:});

%% Save
save([dirs.output,'/stat_lowGammaInteraction.mat'],'stat_lowGammaInteraction','-v7.3');
disp('done')

%% High gamma
cfg = [];
cfg.frequency           = highgamma_band;
cfg.avgoverchan         = 'no';             % Analyse all channels
cfg.avgoverfreq         = 'yes';            % Collapse frequency bins

cfg.method              = 'montecarlo';
cfg.channel             = 'MEGGRAD';

cfg.design = [ones(1,length(PD_hfDiffbs)) 2*ones(1,length(ctrl_hfDiffbs))];
cfg.design = [cfg.design; 1:length(cfg.design)];

cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.neighbours          = neighbours;
cfg.clusteralpha        = 0.05;
cfg.minnbchan           = 2;
cfg.clusterstatistic    = 'maxsum';

cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

stat_hiGammaInteraction = ft_freqstatistics(cfg, PD_hfDiffbs{:}, ctrl_hfDiffbs{:});

%% Save
save([dirs.output,'/stat_hiGammaInteraction.mat'],'stat_hiGammaInteraction','-v7.3');
disp('done')

%% ########################################################################
% Theta
% #########################################################################
%% Load data and prepare
cd(dirs.output);
load('all_TFRlog.mat')
disp('data loaded')

%% Prepare neighbours structure
cfg = [];
cfg.layout  = 'neuromag306cmb.lay';
cfg.method  = 'distance';
neighbours = ft_prepare_neighbours(cfg, data_PD1log{1});

%% Prepare data
cfg = [];
cfg.baseline        = [-inf -0.2];
cfg.baselinetype    = 'absolute';
cfg.parameter       = 'powspctrm';

data_PD1_bs = cell(1,length(PD_subs));
data_PD2_bs = cell(1,length(PD_subs));
data_ctrl1_bs = cell(1,length(PD_subs));
data_ctrl2_bs = cell(1,length(PD_subs));


for ii = 1:length(PD_subs)
    data_PD1_bs{ii} = ft_freqbaseline(cfg,data_PD1log{ii});
    data_PD2_bs{ii} = ft_freqbaseline(cfg,data_PD2log{ii});
end

for ii = 1:length(ctrl_subs)
    data_ctrl1_bs{ii} = ft_freqbaseline(cfg, data_ctrl1log{ii});
    data_ctrl2_bs{ii} = ft_freqbaseline(cfg, data_ctrl2log{ii});
end

disp('done')

%% Difference between 1st sessions
cfg = [];
cfg.frequency           = theta_band;
cfg.avgoverchan         = 'no';             % Analyse all channels
cfg.avgoverfreq         = 'yes';            % collapse frequency bins

cfg.method              = 'montecarlo';
cfg.channel             = 'MEGGRAD';

cfg.design = [ones(1,length(data_PD1_bs)) 2*ones(1,length(data_ctrl1_bs))];
cfg.design = [cfg.design; 1:length(cfg.design)];

cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.neighbours          = neighbours;
cfg.clusteralpha        = 0.05;
cfg.minnbchan           = 2;
cfg.clusterstatistic = 'maxsum';

cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

stat_theta = ft_freqstatistics(cfg, data_PD1_bs{:}, data_ctrl1_bs{:});
disp('done')

%% Save
save([dirs.output,'/stat_theta_all.mat'],'stat_theta','-v7.3');
disp('DONE');

%% Interaction
PD_thetaDiff = cell(1,length(PD_subs));
ctrl_thetaDiff = cell(1,length(ctrl_subs));
PD_thetaDiffbs = PD_thetaDiff;
ctrl_thetaDiffbs = ctrl_thetaDiff;

cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'x1-x2';
for ii = 1:length(data_PD1log)
    PD_thetaDiff{ii} = ft_math(cfg, data_PD2log{ii}, data_PD1log{ii});
end

for ii = 1:length(data_ctrl1log)
    ctrl_thetaDiff{ii} = ft_math(cfg,data_ctrl2log{ii}, data_ctrl1log{ii});
end

cfg                 = [];
cfg.baseline        = [-inf -0.2];
cfg.baselinetype    = 'absolute';
cfg.parameter       = 'powspctrm';

for ii = 1:length(PD_thetaDiffbs)
    PD_thetaDiffbs{ii} = ft_freqbaseline(cfg, PD_thetaDiff{ii});
end

for ii = 1:length(ctrl_thetaDiff)
    ctrl_thetaDiffbs{ii} = ft_freqbaseline(cfg, ctrl_thetaDiff{ii});
end

disp('done')

%% Stat
cfg = [];
cfg.frequency           = theta_band;
cfg.avgoverchan         = 'no';         % Analyse all channels
cfg.avgoverfreq         = 'yes';         % collapse frequency bins

cfg.method              = 'montecarlo';
cfg.channel             = 'MEGGRAD';

cfg.design = [ones(1,length(PD_thetaDiffbs)) 2*ones(1,length(ctrl_thetaDiffbs))];
cfg.design = [cfg.design; 1:length(cfg.design)];

cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.neighbours          = neighbours;
cfg.clusteralpha        = 0.05;
cfg.minnbchan           = 2;
cfg.clusterstatistic = 'maxsum';

cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

stat_theta_interaction = ft_freqstatistics(cfg, PD_thetaDiffbs{:}, ctrl_thetaDiffbs{:});
disp('done')

%% save
save([dirs.output,'/stat_theta_interaction.mat'],'stat_theta_interaction','-v7.3');
disp('DONE');

%% Load previous data
% Beta
load([dirs.output,'/beta_stats.mat']);
disp('done');

% High gamma
load([dirs.output,'/stat_hiGamma_avg.mat']);
load([dirs.output,'/stat_hiGammaInteraction.mat']);
disp('done')

% Low gamma
load([dirs.output,'/stat_lowGamma_avg.mat']);
load([dirs.output,'/stat_lowGammaInteraction.mat']);
disp('done')

% Theta
load([dirs.output,'/stat_theta_all.mat']);
load([dirs.output,'/stat_theta_interaction.mat']);
disp('DONE');


