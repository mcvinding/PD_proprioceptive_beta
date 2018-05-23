%%%%%%%%%%%%% Make plots for publication for PD-proj: beta part %%%%%%%%%%%
addpath /home/mikkel/PD_motor/global_scripts
[dirs, sub_info, lh_subs] = PD_proj_setup('rebound');
addpath /home/mikkel/matlab/export_fig/
addpath /home/mikkel/matlab/align_Ylabels/

cd(dirs.megDir);
dirs.export = '/home/mikkel/PD_motor/rebound/export/publication';

subs = dir(dirs.megDir);                                %Find subjects in folder
subs = {subs([subs.isdir]).name};                       %Make list
subs = subs(~(strcmp('.',subs)|strcmp('..',subs)));     %Remove dots

badsubs = {'0393'}; % Too few trials/bad data
subs(strcmp(badsubs,subs)) = [];

PD_subs = intersect(sub_info.PD,subs);
ctrl_subs = intersect(sub_info.ctrl,subs);

%% ################ TFR FIGURE ##################
%% load data
cd(dirs.output);
load([dirs.output,'/betaChan1bs.mat']);
load([dirs.output,'/acc_ERP2.mat']);
load([dirs.output,'/emg_raw.mat']);
load([dirs.output,'/beta_stats.mat'])
cd(dirs.export);

disp('data loaded');

%% Grand average
avgTFR.PD_beta1 = ft_freqgrandaverage([],PD_beta1bs{:});
avgTFR.PD_beta2 = ft_freqgrandaverage([],PD_beta2bs{:});
avgTFR.ctrl_beta1 = ft_freqgrandaverage([],ctrl_beta1bs{:});
avgTFR.ctrl_beta2 = ft_freqgrandaverage([],ctrl_beta2bs{:});

acc_PD1e = acc_PD1(~cellfun('isempty',acc_PD1)); %Remove empty field
% accPow_PD2e = accPow_PD2(~cellfun('isempty',acc_PD1)); %Remove empty field
accAvg.PD_beta1 = ft_timelockgrandaverage([],acc_PD1e{:});
accAvg.PD_beta2 = ft_timelockgrandaverage([],acc_PD2{:});
accAvg.ctrl_beta1 = ft_timelockgrandaverage([],acc_ctrl1{:});
accAvg.ctrl_beta2 = ft_timelockgrandaverage([],acc_ctrl2{:});

accAvg.PD_beta1.avg = zscore(accAvg.PD_beta1.avg);
accAvg.PD_beta2.avg = zscore(accAvg.PD_beta2.avg);
accAvg.ctrl_beta1.avg = zscore(accAvg.ctrl_beta1.avg);
accAvg.ctrl_beta2.avg = zscore(accAvg.ctrl_beta2.avg);

emg_avg.PD_beta1 = ft_timelockgrandaverage([],PD_emg1{:});
emg_avg.PD_beta2 = ft_timelockgrandaverage([],PD_emg2{:});
emg_avg.ctrl_beta1 = ft_timelockgrandaverage([],ctrl_emg1{:});
emg_avg.ctrl_beta2 = ft_timelockgrandaverage([],ctrl_emg2{:});

conds = fields(avgTFR);

for ii = 1:length(conds)
    emg_avg.(conds{ii}).avg = emg_avg.(conds{ii}).avg*1000; % Change unit mV
end

%% Make seperate plots per group x session
conds = fieldnames(avgTFR);
% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu'))) 

cfg = [];
cfg.zlim        = [-0.30 0.30];
cfg.xlim        = [-.5 2.5];
cfg.colorbar    = 'no';
cfg.colormap    = 'jet';
cfg.interactive = 'no';

cfgTFR = cfg;
cfgACC = cfg;
cfgACC.ylim     = [-1 25];
cfgACC.fontsize = 5;
cfgACC.graphcolor = [0 0 0];
cfgEMG = cfg;
cfgEMG.ylim     = [-1 1];
cfgEMG.fontsize = 5;
cfgEMG.graphcolor = [0 0 0];

titles = {
    'PD OFF medication (session 1)'
    'PD ON medication (session 2)'
    'Healthy controls (session 1)'
    'Healthy controls (session 2)'};

cd(dirs.export);
for ii = 1:length(conds)
    fig = figure;
    subplot(20,1,1:13); ft_singleplotTFR(cfgTFR, avgTFR.(conds{ii}));
    ax1 = gca();
    ft_plot_line([0 0],ax1.YLim, 'color','k','linestyle','--');
    box on
    title(titles{ii},'fontsize',12);
    set(gca, 'LineWidth', 2);

    y1 = ylabel('Freq. (Hz)'); xlabel('Time (s)'); 
%     colormap(flipud(brewermap(64,'RdBu'))) 

    subplot(20,1,16:17); ft_singleplotER(cfgACC, accAvg.(conds{ii}));
    y2 = ylabel({'Acceleration';'(z-score)'},'fontsize',6);
    set(gca,'xtick',[],'xticklabel','','XColor','none','LineWidth', 1 )
    set(y2, 'HorizontalAlignment','center','VerticalAlignment','bottom')
    title(' ');
    
    subplot(20,1,19:20); ft_singleplotER(cfgEMG, emg_avg.(conds{ii}));
    y3 = ylabel({'EMG';'(mV)'},'fontsize',6); 
    set(gca,'LineWidth', 1)
    set(y3,'HorizontalAlignment','center','VerticalAlignment','bottom'); %    y2.Position)
    title(' ');
    
%     set(findobj(gcf,'type','axes'), 'LineWidth', 1);
    set(fig,'PaperPosition', [0 0 4 2], 'color','w');
    
%     savefig(['TFR_',conds{ii},'.fig']);
%     print('-dpng','-r500',['TFR_',conds{ii},'.png'])
%     print('-painters','-dpdf',['TFR_',conds{ii},'.pdf'])

    export_fig(['TFR_',conds{ii},'exp.png'], '-r500', '-p0.05', '-CMYK')
    export_fig(['TFR_',conds{ii},'exp.pdf'], '-r500', '-p0.05', '-CMYK', '-pdf')

end
disp('DONE');


%% Extract beta traces an plot
cd(dirs.export);
timeDim = avgTFR.PD_beta1.time;

cfg = [];
cfg.avgoverchan = 'no';
cfg.frequency = [12 25];
cfg.avgoverfreq = 'yes';
cfg.avgoverchan = 'yes';

BdataPD1 = ft_selectdata(cfg, avgTFR.PD_beta1);
BtracePD1 = squeeze(BdataPD1.powspctrm);
BdataPD2 = ft_selectdata(cfg, avgTFR.PD_beta2);
BtracePD2 = squeeze(BdataPD2.powspctrm);

BdataCtrl1 = ft_selectdata(cfg,avgTFR.ctrl_beta1);
BtraceCtrl1 = squeeze(BdataCtrl1.powspctrm);
BdataCtrl2 = ft_selectdata(cfg,avgTFR.ctrl_beta2);
BtraceCtrl2 = squeeze(BdataCtrl2.powspctrm);

% Get cluster
clustStart = stat_beta1.time(find(sum(squeeze(stat_beta1.mask)), 1,'first'));
clustEnd = stat_beta1.time(find(sum(squeeze(stat_beta1.mask)), 1,'last'));

% make plot
fig = figure('rend','painters','pos',[10 10 1000 600]); hold on
set(fig,'PaperPosition', [0 0 4 2], 'color','w');

patch([clustStart clustEnd clustEnd clustStart],[-0.3 -0.3 0.2 0.2],[.1,.1,.1],'FaceAlpha',0.2,'EdgeColor','none')
p1 = plot(timeDim,squeeze(BtracePD1),'b-','LineWidth',2);
p2 = plot(timeDim,squeeze(BtracePD2),'b--','LineWidth',2);
p3 = plot(timeDim,squeeze(BtraceCtrl1),'r-','LineWidth',2);
p4 = plot(timeDim,squeeze(BtraceCtrl2),'r--','LineWidth',2);
axis([-.5, 2.5, -0.3, 0.2])
line([0 0],[-0.3, 0.2],'color',[0.5 .5 .5],'LineStyle','--','LineWidth',2)
xlabel('Time (s)','fontsize',12);
ylabel('Relative change','fontsize',12)
title('Beta-band spectral evolution','fontsize',16);
set(gca, 'LineWidth', 2,'fontweight','bold');
legend([p1,p2,p3,p4],{'PD OFF','PD ON','HC Session1','HC Session2'},'Location','SouthEast')
legend BOXOFF
export_fig('beta_evo.png', '-r500', '-p0.05', '-CMYK', '-png', '-transparent')
% % close


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define bands
theta_band      = [4 7];
betamu_band     = [8 30];
lowgamma_band   = [31 49];
highgamma_band  = [51 100];
gamma_band      = [31 100];

%% High freq
% Load data
cd(dirs.output);
load('all_hfreq_bs')
load('stat_lowGamma_avg.mat');
load('stat_hiGamma_avg');
load('all_TFRlog.mat');
load('stat_theta_all.mat')
disp('done');

%% Grand average
cfg = [];
cfg.foilim  = lowgamma_band;

avgLowGam.PD1 = ft_freqgrandaverage(cfg, data_PD1_hf_bs{:});
avgLowGam.PD2 = ft_freqgrandaverage(cfg, data_PD2_hf_bs{:});
avgLowGam.ctrl1 = ft_freqgrandaverage(cfg, data_ctrl1_hf_bs{:});
avgLowGam.ctrl2 = ft_freqgrandaverage(cfg, data_ctrl2_hf_bs{:});

avgLowGam.PD1.mask = stat_lowGamma_avg.mask;
avgLowGam.PD2.mask = stat_lowGamma_avg.mask;
avgLowGam.ctrl1.mask = stat_lowGamma_avg.mask;
avgLowGam.ctrl2.mask = stat_lowGamma_avg.mask;

cfg = [];
cfg.foilim  = highgamma_band;

avgHiGam.PD1 = ft_freqgrandaverage(cfg, data_PD1_hf_bs{:});
avgHiGam.PD2 = ft_freqgrandaverage(cfg, data_PD2_hf_bs{:});
avgHiGam.ctrl1 = ft_freqgrandaverage(cfg, data_ctrl1_hf_bs{:});
avgHiGam.ctrl2 = ft_freqgrandaverage(cfg, data_ctrl2_hf_bs{:});

avgHiGam.PD1.mask = stat_hiGamma_avg.mask;
avgHiGam.PD2.mask = stat_hiGamma_avg.mask;
avgHiGam.ctrl1.mask = stat_hiGamma_avg.mask;
avgHiGam.ctrl2.mask = stat_hiGamma_avg.mask;

% Prepare theta band
cfg = [];
cfg.baseline        = [-inf -0.2];
cfg.baselinetype    = 'absolute';
cfg.parameter       = 'powspctrm';

dataT_PD1_bs = cell(1,length(PD_subs));
dataT_PD2_bs = cell(1,length(PD_subs));
dataT_ctrl1_bs = cell(1,length(PD_subs));
dataT_ctrl2_bs = cell(1,length(PD_subs));

for ii = 1:length(PD_subs)
    dataT_PD1_bs{ii} = ft_freqbaseline(cfg,data_PD1log{ii});
    dataT_PD2_bs{ii} = ft_freqbaseline(cfg,data_PD2log{ii});
end

for ii = 1:length(ctrl_subs)
    dataT_ctrl1_bs{ii} = ft_freqbaseline(cfg, data_ctrl1log{ii});
    dataT_ctrl2_bs{ii} = ft_freqbaseline(cfg, data_ctrl2log{ii});
end

cfg = [];
cfg.foilim  = theta_band;

avgTheta.PD1 = ft_freqgrandaverage(cfg, dataT_PD1_bs{:});
avgTheta.PD2 = ft_freqgrandaverage(cfg, dataT_PD2_bs{:});
avgTheta.ctrl1 = ft_freqgrandaverage(cfg, dataT_ctrl1_bs{:});
avgTheta.ctrl2 = ft_freqgrandaverage(cfg, dataT_ctrl2_bs{:});

avgTheta.PD1.mask = stat_theta.mask;
avgTheta.PD2.mask = stat_theta.mask;
avgTheta.ctrl1.mask = stat_theta.mask;
avgTheta.ctrl2.mask = stat_theta.mask;

disp('done')

%% Plot: Multiplots
cfg = [];
cfg.layout          = 'neuromag306cmb.lay';
cfg.showlabels      = 'no';
cfg.interactive     = 'yes';
cfg.showcomment     = 'no';
cfg.maskparameter   = 'mask';  % use the thresholded probability to mask the data
cfg.maskstyle       = 'box';
% cfg.maskalpha       = .1;
cfg.xlim            = [-.5 2.5];
cfg.ylim            = [-0.15 0.15];
% cfg.linestyle       = {'-',':','-',':'};
% cfg.graphcolor      = {'blue','blue','red','red'};
figure('rend','painters','pos',[10 10 800 600]);
ft_multiplotER(cfg, avgLowGam.PD1,avgLowGam.PD2,avgLowGam.ctrl1,avgLowGam.ctrl2);
legend('PD OFF','PD ON','Ctrl. Session1','Ctrl. Session2')

% export_fig(['sens_lowGamma.png'], '-r500', '-p0.05', '-CMYK', '-png', '-transparent')

figure('rend','painters','pos',[10 10 800 600]);
ft_multiplotER(cfg, avgHiGam.PD1,avgHiGam.PD2,avgHiGam.ctrl1,avgHiGam.ctrl2);
legend('PD OFF','PD ON','Ctrl. Session1','Ctrl. Session2')

% export_fig(['sens_hiGamma.png'], '-r500', '-p0.05', '-CMYK', '-png', '-transparent')

figure('rend','painters','pos',[10 10 800 600]);
ft_multiplotER(cfg, avgTheta.PD1,avgTheta.PD2,avgTheta.ctrl1,avgTheta.ctrl2);
legend('PD OFF','PD ON','Ctrl. Session1','Ctrl. Session2')

% export_fig(['sens_theta.png'], '-r500', '-p0.05', '-CMYK', '-png', '-transparent')

%% Representative single channels 
% Low gamma
cfg = [];
cfg.channel         = 'MEG0632+0633';
cfg.latency         = [-.5 2.5];
cfg.avgoverfreq     = 'yes';

temp = ft_selectdata(cfg, avgLowGam.PD1);
LGtracePD1 = squeeze(temp.powspctrm);
temp = ft_selectdata(cfg, avgLowGam.PD2);
LGtracePD2 = squeeze(temp.powspctrm);
temp = ft_selectdata(cfg, avgLowGam.ctrl1);
LGtraceCtrlD1 = squeeze(temp.powspctrm);
temp = ft_selectdata(cfg, avgLowGam.ctrl2);
LGtracectrl2 = squeeze(temp.powspctrm);

timeDim = temp.time;

% Get cluster
clustStart = stat_lowGamma_avg.time(find(sum(squeeze(stat_lowGamma_avg.mask)), 1,'first'));
clustEnd = stat_lowGamma_avg.time(find(sum(squeeze(stat_lowGamma_avg.mask)), 1,'last'));

% Make plot
fig = figure('rend','painters','pos',[10 10 1000 600]); hold on
set(fig,'PaperPosition', [0 0 4 2], 'color','w');

patch([clustStart clustEnd clustEnd clustStart],[-0.1 -0.1 0.1 0.1],[.1,.1,.1],'FaceAlpha',0.2,'EdgeColor','none')
h1 = plot(timeDim,LGtracePD1,'b-','LineWidth',2);
h2 = plot(timeDim,LGtracePD2,'b--','LineWidth',2);
h3 = plot(timeDim,LGtraceCtrlD1,'r-','LineWidth',2);
h4 = plot(timeDim,LGtracectrl2,'r--','LineWidth',2);
axis([-.5, 2.5, -0.1, 0.1])
line([0 0],[-0.3, 0.2],'color',[0.5 .5 .5],'LineStyle','--','LineWidth',2)
xlabel('Time (s)','fontsize',12);
ylabel('Relative change','fontsize',12)
title('Low gamma spectral evolution','fontsize',16);
set(gca, 'LineWidth', 2,'fontweight','bold');
legend([h1,h2,h3,h4], {'PD OFF','PD ON','HC Session1','HC Session2'},'Location','SouthEast')
legend BOXOFF
export_fig('lowGamma_evo.png', '-r500', '-p0.05', '-CMYK', '-png', '-transparent')
% % close

% Theta
cfg = [];
cfg.channel         = 'MEG1122+1123';
cfg.latency         = [-.5 2.5];
cfg.avgoverfreq     = 'yes';

temp = ft_selectdata(cfg, avgTheta.PD1);
THtracePD1 = squeeze(temp.powspctrm);
temp = ft_selectdata(cfg, avgTheta.PD2);
THtracePD2 = squeeze(temp.powspctrm);
temp = ft_selectdata(cfg, avgTheta.ctrl1);
THtraceCtrl1 = squeeze(temp.powspctrm);
temp = ft_selectdata(cfg, avgTheta.ctrl2);
THtracectrl2 = squeeze(temp.powspctrm);

timeDim = temp.time;

% Get cluster
clustStart = stat_lowGamma_avg.time(find(sum(squeeze(stat_theta.mask)), 1,'first'));
clustEnd = stat_lowGamma_avg.time(find(sum(squeeze(stat_theta.mask)), 1,'last'));

% Make plot
fig = figure('rend','painters','pos',[10 10 1000 600]); hold on
set(fig,'PaperPosition', [0 0 4 2], 'color','w');

patch([clustStart clustEnd clustEnd clustStart],[-0.1 -0.1 0.1 0.1],[.1,.1,.1],'FaceAlpha',0.2,'EdgeColor','none')
h1 = plot(timeDim,THtracePD1,'b-','LineWidth',2);
h2 = plot(timeDim,THtracePD2,'b--','LineWidth',2);
h3 = plot(timeDim,THtraceCtrl1,'r-','LineWidth',2);
h4 = plot(timeDim,THtracectrl2,'r--','LineWidth',2);
axis([-.5, 2.5, -0.1, 0.1])
line([0 0],[-0.3, 0.2],'color',[0.5 .5 .5],'LineStyle','--','LineWidth',2)
xlabel('Time (s)','fontsize',12);
ylabel('Relative change','fontsize',12)
title('Theta spectral evolution','fontsize',16);
set(gca, 'LineWidth', 2,'fontweight','bold');
legend([h1,h2,h3,h4], {'PD OFF','PD ON','HC Session1','HC Session2'},'Location','SouthEast')
legend BOXOFF
export_fig('theta_evo.png', '-r500', '-p0.05', '-CMYK', '-png', '-transparent')
% close

%% Plot: Topo plots
timestep = 0.5;		% timestep between time windows for each subplot (in seconds)
sampling_rate = 100;	% Time-axis has a temporal resolution of 100 Hz
sample_count = length(stat_theta.time);
					% number of temporal samples in the statistics object
j = [0:timestep:2.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples
m = m(stat_theta.time(m)>=0);

% pos = stat_theta.posclusterslabelmat == 1;
% neg = stat_theta.negclusterslabelmat == 1;

conds = {'PD1','PD2','ctrl1','ctrl2'};

figure('rend','painters','pos',[10 10 1200 600])
for c = 1:4
    for k = 1:5;
        disp(k+(c-1)*5)
        subplot(4,5,k+(c-1)*5);
        cfg = [];   
        cfg.xlim=[j(k) j(k+1)];   % time interval of the subplot
        cfg.zlim = [-0.15 .15];
        cfg.comment        = 'no';   
        cfg.commentpos     = 'title';   
        cfg.layout         = 'neuromag306cmb.lay';
        cfg.interactive 	= 'no';
        cfg.colormap       = 'jet';
        cfg.marker         = 'off';
        ft_topoplotER(cfg, avgTheta.(conds{c}));   
    end
end

cd(dirs.export);
export_fig(['Topo_theta.png'], '-r500', '-p0.05', '-CMYK', '-png', '-transparent')
export_fig(['Topo_theta.pdf'], '-r500', '-p0.05', '-CMYK', '-pdf')

figure('rend','painters','pos',[10 10 1200 600])
for c = 1:4
    for k = 1:5;
        disp(k+(c-1)*5)
        subplot(4,5,k+(c-1)*5);
        cfg = [];   
        cfg.xlim=[j(k) j(k+1)];   % time interval of the subplot
        cfg.zlim = [-0.15 .15];
        cfg.comment        = 'no';   
        cfg.commentpos     = 'title';   
        cfg.layout         = 'neuromag306cmb.lay';
        cfg.interactive 	= 'no';
        cfg.colormap       = 'jet';
        cfg.marker         = 'off';
        ft_topoplotER(cfg, avgLowGam.(conds{c}));   
    end
end

cd(dirs.export);
export_fig(['Topo_loGamma.png'], '-r500', '-p0.05', '-CMYK', '-png','-transparent')
export_fig(['Topo_loGamma.pdf'], '-r500', '-p0.05', '-CMYK', '-pdf')

figure('rend','painters','pos',[10 10 1200 600])
for c = 1:4
    for k = 1:5;
        disp(k+(c-1)*5)
        subplot(4,5,k+(c-1)*5);
        cfg = [];   
        cfg.xlim=[j(k) j(k+1)];   % time interval of the subplot
        cfg.zlim = [-0.15 .15];
        cfg.comment        = 'no';   
        cfg.commentpos     = 'title';   
        cfg.layout         = 'neuromag306cmb.lay';
        cfg.interactive 	= 'no';
        cfg.colormap       = 'jet';
        cfg.marker         = 'off';
        ft_topoplotER(cfg, avgHiGam.(conds{c}));   
    end
end

cd(dirs.export);
export_fig(['Topo_hiGamma.png'], '-r500', '-p0.05', '-CMYK', '-png','-transparent')
export_fig(['Topo_hiGamma.pdf'], '-r500', '-p0.05', '-CMYK', '-pdf')

%% DELETE
% Below here not used...

% Single subs PD 1
subs = 1:length(PD_beta1bs);
cfg = [];
% cfg.baseline = [-inf 0];
cfg.zlim = [-0.5 0.5];
cfg.xlim = [-.5 2.5];
cfg.ylim = [4 30];
% cfg.baselinetype = 'relative';
cfg.colorbar = 'no';

figure;
for ss = 1:length(subs)
    subplot(4,3,ss); ft_singleplotTFR(cfg,PD_beta1bs{ss});
    title(subs(ss))
end
savefig(['/home/mikkel/PD_motor/rebound/export/TFR_PD1_subs.fig']);
print(['/home/mikkel/PD_motor/rebound/export/TFR_PD1_subs.png'], '-dpng')
    
% Single subs PD 2
figure;
for ss = 1:length(subs)
    subplot(4,3,ss); ft_singleplotTFR(cfg,PD_beta1bs{ss});
    title(subs(ss))
end
savefig(['/home/mikkel/PD_motor/rebound/export/TFR_PD2_subs.fig']);
print(['/home/mikkel/PD_motor/rebound/export/TFR_PD2_subs.png'], '-dpng')

% Single subs Ctrl 1
subs = 1:length(ctrl_beta1bs);

figure;
for ss = 1:length(subs)
    subplot(6,3,ss); ft_singleplotTFR(cfg,ctrl_beta1bs{ss});
    title(subs(ss))
end
savefig(['/home/mikkel/PD_motor/rebound/export/TFR_ctrl1_subs.fig']);
print(['/home/mikkel/PD_motor/rebound/export/TFR_ctrl1_subs.png'], '-dpng')

% Single subs PD 2
figure;
for ss = 1:length(subs)
    subplot(6,3,ss); ft_singleplotTFR(cfg,slct_data_ctrl2.(subs{ss}));
    title(subs(ss))
end
savefig(['/home/mikkel/PD_motor/rebound/export/TFR_ctrl2_subs.fig']);
print(['/home/mikkel/PD_motor/rebound/export/TFR_ctrl2_subs.png'], '-dpng')


export_fig