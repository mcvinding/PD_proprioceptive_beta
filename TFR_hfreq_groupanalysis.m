%%%%% GET TFR GRAND AVERAGE FOR HIGH FREQ DATA
%MOVED TO OTHER SCRIPT - DELETE THIS ONE

addpath /home/mikkel/PD_motor/global_scripts
[dirs, sub_info, lh_subs] = PD_proj_setup('rebound');

cd(dirs.megDir);

subs = dir(dirs.megDir);                                %Find subjects in folder
subs = {subs([subs.isdir]).name};                       %Make list
subs = subs(~(strcmp('.',subs)|strcmp('..',subs)));     %Remove dots
PD_subs = intersect(sub_info.PD,subs);
ctrl_subs = intersect(sub_info.ctrl,subs);

%% 

% MOVED


%% Below here not copied anywhere
%% Single plots
cfg = [];
% cfg.baseline = [-inf inf];
% cfg.baselinetype = 'relative';
cfg.layout = 'neuromag306cmb.lay';
cfg.zlim = [-.3 .3];
figure; ft_multiplotTFR(cfg,avgTFR_hfreq.PD1); title('TFR session 1');
figure; ft_multiplotTFR(cfg,avgTFR_hfreq.ctrl1); title('TFR session 2');



%% Plots
load([dirs.output,'/avgTFR_hfreq.mat']);
disp('done');

cfg = [];
% cfg.baseline = [-inf 0];
% cfg.baselinetype = 'relchange';
cfg.layout = 'neuromag306cmb.lay';
% cfg.zlim = [-3e-24 3e-24];
% cfg.zlim = [.6 1.6];

% cfg.zlim = [.2 1.8];
figure; ft_multiplotTFR(cfg,avgTFRbase_hfreq.PD1); title('PD session 1');
figure; ft_multiplotTFR(cfg,avgTFR_hfreq.PD2); title('PD session 2');
figure; ft_multiplotTFR(cfg,avgTFR_hfreq.ctrl1); title('ctrl session 1');
figure; ft_multiplotTFR(cfg,avgTFR_hfreq.ctrl2); title('ctrl session 2');

%%
conds = fields(avgTFR_hfreq);
cfg = [];
cfg.zlim = [.6 1.4];
cfg.xlim = [-.5 2];
cfg.ylim = [4 35];
cfg.baseline = [-inf inf];
cfg.baselinetype = 'absolute';
cfg.colorbar = 'no';

% All in one plot
figure;
for ii = 1:length(conds)
    subplot(2,2,ii); ft_singleplotTFR(cfg,avgTFR_hfreq.(conds{ii}));
    title(conds{ii});
    ylabel('Freq. (Hz)')
    ylabel('Time (s)')
end

% Seperate plots
for ii = 1:length(conds)
    figure;
    ft_singleplotTFR(cfg,avgTFR_hfreq.(conds{ii}));
    title(conds{ii});
    ylabel('Freq. (Hz)'); xlabel('Time (s)')
    savefig(['/home/mikkel/PD_motor/rebound/export/TFR_',conds{ii},'.fig']);
    print(['/home/mikkel/PD_motor/rebound/export/TFR_',conds{ii},'.png'], '-dpng')
    close
end

% Single subs PD 1
subs = fieldnames(slct_data_PD1);
cfg = [];
cfg.baseline = [-inf inf];
cfg.zlim = [.4 1.6];
cfg.xlim = [-.5 2.5];
cfg.ylim = [4 35];
cfg.baselinetype = 'relative';
cfg.colorbar = 'no';

figure;
for ss = 1:length(subs)
    subplot(4,3,ss); ft_singleplotTFR(cfg,slct_data_PD1.(subs{ss}));
    title(subs(ss))
end
savefig(['/home/mikkel/PD_motor/rebound/export/TFR_PD1_subs.fig']);
print(['/home/mikkel/PD_motor/rebound/export/TFR_PD1_subs.png'], '-dpng')
    
% Single subs PD 2
figure;
for ss = 1:length(subs)
    subplot(4,3,ss); ft_singleplotTFR(cfg,slct_data_PD2.(subs{ss}));
    title(subs(ss))
end
savefig(['/home/mikkel/PD_motor/rebound/export/TFR_PD2_subs.fig']);
print(['/home/mikkel/PD_motor/rebound/export/TFR_PD2_subs.png'], '-dpng')

% Single subs Ctrl 1
subs = fieldnames(slct_data_ctrl1);

figure;
for ss = 1:length(subs)
    subplot(6,3,ss); ft_singleplotTFR(cfg,slct_data_ctrl1.(subs{ss}));
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

%% Statistics %%
load([dirs.output,'/avgTFR_hfreq.mat'])
disp('done');
%% !MOVED!
% % baseline correct
% cfg = [];
% cfg.baseline        = [-inf 0];
% cfg.baselinetype    = 'relative';
% avgTFRbase_hfreq.PD1    = ft_freqbaseline(cfg,avgTFR_hfreq.PD1);
% avgTFRbase_hfreq.PD2    = ft_freqbaseline(cfg,avgTFR_hfreq.PD2);
% avgTFRbase_hfreq.ctrl1  = ft_freqbaseline(cfg,avgTFR_hfreq.ctrl1);
% avgTFRbase_hfreq.ctrl2  = ft_freqbaseline(cfg,avgTFR_hfreq.ctrl2);

%Log transform
cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'log10';
TFR_PD1log       = ft_math(cfg, avgTFR_hfreq.PD1);
TFR_ctrl1log     = ft_math(cfg, avgTFR_hfreq.ctrl1);

% Difference based on 1st session !MOVED!
% cfg = [];
% cfg.latency             = [-.5 2.5];
% cfg.frequency           = [31 100];
% cfg.avgoverchan         = 'no';
% cfg.method              = 'montecarlo';
% cfg.design = [ones(1,size(avgTFR_hfreq.PD1.powspctrm,1)) 2*ones(1,size(avgTFR_hfreq.ctrl1.powspctrm,1))];
% cfg.desing = [cfg.design; 1:length(cfg.design)];
% cfg.statistic           = 'ft_statfun_indepsamplesT';
% cfg.numrandomization    = 1000;
% cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
% cfg.tail                = 0;
% cfg.computeprob         = 'yes';
% cfg.computecritval      = 'yes';
% cfg.alpha               = .05;
% 
% stats = ft_freqstatistics(cfg, avgTFR_hfreq.PD1, avgTFR_hfreq.ctrl1);

% Interaction with session
%Get difference
cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = '(x1-x2)/(x1+x2)';
% cfg.operation = 'x1-x2';
TFR_PD_diff    = ft_math(cfg, avgTFR_hfreq.PD1,avgTFR_hfreq.PD2);
TFR_ctrl_diff    = ft_math(cfg, avgTFR_hfreq.ctrl1,avgTFR_hfreq.ctrl2);

cfg = [];
cfg.latency             = [-.5 2.5];
cfg.frequency           = [4 40];
cfg.avgoverchan         = 'yes';
cfg.method              = 'montecarlo';
cfg.design = [ones(1,size(avgTFR_hfreq.PD1.powspctrm,1)) 2*ones(1,size(avgTFR_hfreq.ctrl1.powspctrm,1))];
cfg.desing = [cfg.design; 1:length(cfg.design)];
cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';

stats_interact = ft_freqstatistics(cfg, TFR_PD_diff, TFR_ctrl_diff);

% Correlation with UPDRS score (PD only)
cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'log10';
% TFR_PD1log    = ft_math(cfg, avgTFR.PD1); %Already done above!
TFR_PD2log    = ft_math(cfg, avgTFR_hfreq.PD2);

% Scores from behavious assesment (manually entered here)
scores1 = [38 25 31 10 29 22 27 29 28 21 41 61];
scores2 = [21 7 16 5 19 6 9 14 15 12 31 39];
subs = [1:12, 1:12];
reg_design = [scores1,scores2; subs];

% Difference based UPDRS
cfg = [];
cfg.latency             = [-.5 2.5];
cfg.frequency           = [4 40];
cfg.avgoverchan         = 'yes';
cfg.method              = 'montecarlo';
cfg.design              =  reg_design;
cfg.statistic           = 'ft_statfun_depsamplesregrT';
cfg.numrandomization    = 1000;
cfg.ivar                = 1;            % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2;
cfg.tail                = 0;
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .05;

stats_UPDRS = ft_freqstatistics(cfg, TFR_PD1log, TFR_PD2log);

% Save stat outputs
save([dirs.output,'/TFRstats_hfreq.mat'],'stats') %,'stats_interact','stats_UPDRS');
disp('done')

%% PLOT STATS %%
load([dirs.output,'/TFRstats.mat'])
load([dirs.output,'/avgTFR.mat'])
%% cont.
% % Plot main effect (session 1 PD vs. Ctrl.)
% cfg = [];
% cfg.avgoverchan = 'no';
% cfg.latency             = [-.5 2.5];
% cfg.frequency           = [31 100];
% PD1_plt = ft_selectdata(cfg,avgTFRbase_hfreq.PD1);
% CTRL1_plt = ft_selectdata(cfg,avgTFRbase_hfreq.ctrl1);
% 
% PD1_plt.mask = stats.mask;
% CTRL1_plt.mask = stats.mask;
% 
% %Plot stat
% cfg               = [];
% cfg.marker        = 'on';
% cfg.layout = 'neuromag306cmb.lay';
% cfg.parameter     = 'powspctrm';
% cfg.maskparameter = 'mask';  % use the thresholded probability to mask the data
% cfg.maskstyle     = 'opacity';
% cfg.maskalpha       = .1;
% % cfg.baseline = [-inf inf];
% % cfg.baselinetype = 'relative';
% % cfg.xlim = [-.5 2.5];
% % cfg.ylim = [4 35];
% % cfg.zlim = [.5 1.5];
% 
% figure; ft_multiplotTFR(cfg, PD1_plt); title('PD')
% figure; ft_multiplotTFR(cfg, CTRL1_plt); title('Ctrl')

% Interaction
cfg = [];
cfg.parameter      = 'stat';
cfg.colorbar        = 'yes';
cfg.maskparameter   = 'mask';  % use the thresholded probability to mask the data
cfg.maskstyle     = 'outline';
figure; ft_singleplotTFR(cfg, stats_interact);

cfg = [];
cfg.avgoverchan = 'yes';
cfg.latency             = [-.5 2.5];
cfg.frequency           = [4 40];
conds = fields(avgTFR_hfreq);
for ii = 1:length(conds)
    avgTFRplt.(conds{ii}) = ft_selectdata(cfg,avgTFR_hfreq.(conds{ii}));
    avgTFRplt.(conds{ii}).mask = stats_interact.mask;
end

cfg = [];
cfg.avgoverchan = 'yes';
cfg.latency             = [-.5 2.5];
cfg.frequency           = [4 40];
TFR_PD_diff = ft_selectdata(cfg,TFR_PD_diff);
TFR_ctrl_diff = ft_selectdata(cfg,TFR_ctrl_diff);
TFR_PD_diff.mask = stats_interact.mask;
TFR_ctrl_diff.mask = stats_interact.mask;

% All in one plot
cfg               = [];
cfg.marker        = 'on';
cfg.layout = 'neuromag306cmb.lay';
cfg.parameter     = 'powspctrm';
cfg.maskparameter = 'mask';  % use the thresholded probability to mask the data
cfg.maskstyle     = 'opacity';
cfg.maskalpha     = .6;
cfg.baseline = [-inf inf];
cfg.baselinetype    = 'relative';
cfg.colorbar        = 'no';
cfg.xlim            = [-.5 2.5];
cfg.ylim            = [4 35];
cfg.zlim            = [.5 1.5];

conds = fields(avgTFR_hfreq);
figure;
for ii = 1:length(conds)
    subplot(2,2,ii); ft_singleplotTFR(cfg,avgTFRplt.(conds{ii}));
    title(conds{ii});
    ylabel('Freq. (Hz)')
    xlabel('Time (s)')
end

% Plot difference maps
cfg = [];
cfg.baseline = [-inf inf];
cfg.zlim            = [-.08 .05];
cfg.xlim            = [-.5 2.5];
cfg.ylim            = [4 35];
cfg.baselinetype    = 'absolute';
cfg.baseline = [-inf 0];
cfg.colorbar        = 'no';
cfg.maskparameter   = 'mask';  % use the thresholded probability to mask the data
cfg.maskstyle     = 'outline';
figure;
subplot(1,2,1); ft_singleplotTFR(cfg, TFR_ctrl_diff); title('Ctrl');
ylabel('Freq. (Hz)');xlabel('Time (s)');
subplot(1,2,2); ft_singleplotTFR(cfg, TFR_PD_diff); title('PD');
ylabel('Freq. (Hz)');xlabel('Time (s)');

% Correlation w/UPDRS
cfg               = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.parameter           = 'stat';
cfg.maskparameter       = 'mask';  % use the thresholded probability to mask the data
cfg.maskstyle           = 'opacity';
cfg.maskalpha           = .75;
cfg.colorbar            = 'no';
cfg.latency             = [-.5 2.5];
cfg.frequency           = [4 35];
% cfg.baseline = [-inf inf];
% cfg.baselinetype = 'relative';
figure; ft_singleplotTFR(cfg, stats_UPDRS); title('Correlation w/UPDRS-III score')
xlabel('Time (s)'); ylabel('Freq. (Hz)')

avgTFRplt.PD1.mask = stats_UPDRS.mask;
cfg.parameter           = 'powspctrm';
cfg.baseline = [-inf inf];
cfg.baselinetype = 'relative';
cfg.maskalpha           = .50;
cfg.zlim            = [.6 1.4];

figure; ft_singleplotTFR(cfg, avgTFRplt.PD1); title(' ')
xlabel('Time (s)'); ylabel('Freq. (Hz)')

% Correlation cluster plot
% Get data for dotplot
cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverfreq = 'yes';
cfg.latency     = [0.8 1];
cfg.frequency   = [11.5 16];
avg_cluster1 = ft_selectdata(cfg,avgTFR_hfreq.PD1);
avg_cluster2 = ft_selectdata(cfg,avgTFR_hfreq.PD2);
avg_clusterDiff = avg_cluster2.powspctrm-avg_cluster1.powspctrm;
scores_diff = scores1-scores2;

avg_cluster = [avg_cluster1.powspctrm', avg_cluster2.powspctrm'];
avg_cluster = [avg_cluster;reg_design];

figure; scatter(avg_cluster(2,:),avg_cluster(1,:), [], avg_cluster(3,:));
figure; scatter(scores_diff,avg_clusterDiff, 'filled','k');
lsline()
title('UPDRS-III improvement vs. PWR-change');
xlabel('UPDRS-III score diff.'); ylabel('POW diff.');

%% plot tfr %BELOW HERE NOT USED ... REDO PLOTTING!

cfg = [];
cfg.baseline = [-inf inf];
cfg.baselinetype = 'relative';
cfg.layout = 'neuromag306cmb.lay';
cfg.zlim = [.2 1.8];

figure; ft_multiplotTFR(cfg,test); title('TFR session 1');
figure; ft_multiplotTFR(cfg,PD_avg2); title('TFR session 2');

figure; ft_multiplotTFR(cfg,ctrl_avg1); title('TFR session 1');
figure; ft_multiplotTFR(cfg,ctrl_avg2); title('TFR session 2');

% Select data
% NEED L-R dif
% channels = {'MEG0412+0413','MEG0422+0423','MEG0432+0423','MEG0442+0443','MEG1812+1813','MEG1822+1823'};
channels = {'MEG0232+0233','MEG0442+0443','MEG0432+0433','MEG1622+1623','MEG1812+1813','MEG1822+1823'};
channels_L = {'MEG1142+1143','MEG1132+1133','MEG1342+1343','MEG2212+2213','MEG2222+2223','MEG2412+2413'};

cfg = [];
cfg.baseline = [-inf inf];
cfg.baselinetype = 'relative';
cfg.layout = 'neuromag306cmb.lay';
cfg.zlim = [.4 1.6];

% figure; ft_multiplotTFR(cfg, PD_avg1); title('AVG-1')
% figure; ft_multiplotTFR(cfg, PD_avg2); title('AVG-2')
% figure; ft_multiplotTFR(cfg, ctrl_avg1); title('AVG-1-ctrl')
% figure; ft_multiplotTFR(cfg, ctrl_avg2); title('AVG-2-ctrl')

cfg = [];
cfg.baseline = [-inf inf];
cfg.channel = channels;
cfg.baselinetype = 'relative';
cfg.zlim = [.4 1.6];
cfg.xlim = [-.5 2.5];
cfg.ylim = [5 35];
cfg.colorbar = 'no';
figure('Position', [300, 300, 1000, 400])
subplot(1,2,1); ft_singleplotTFR(cfg,PD_avg1); title('PD patients','fontsize',15) 
xlabel('Time (s)'); ylabel('Frequency (Hz)');
% subplot(1,2,1); ft_plot_line([0 0],[-inf, inf],'color','k','linestyle','--');
subplot(1,2,2); ft_singleplotTFR(cfg,ctrl_avg1); title('Healthy controls','fontsize',15)
xlabel('Time (s)');
savefig('/home/mikkel/PD_motor/rebound/export/TFR.fig');
print('/home/mikkel/PD_motor/rebound/export/TFR.pdf', '-dpdf','-bestfit')
% close

% figure; ft_multiplotTFR(cfg, data_PD1.a0327); title('0327')
% figure; ft_multiplotTFR(cfg, data_PD1.a0339); title('0339')
% figure; ft_multiplotTFR(cfg, data_PD1.a0340); title('0340')
% figure; ft_multiplotTFR(cfg, data_PD1.a0352); title('0352')
% figure; ft_multiplotTFR(cfg, data_PD1.a0353); title('0353')
% figure; ft_multiplotTFR(cfg, data_PD1.a0355); title('0355')

%% Extract beta trace
timeDim = ctrl2_bsPow{1}.time;
% idx = PD_avg1.freq > 12 & PD_avg1.freq < 25;
bline = timeDim < 0;

PD1_avg = ft_freqgrandaverage([],PD1_bsPow{:});
ctrl1_avg = ft_freqgrandaverage([],ctrl1_bsPow{:});
PD2_avg = ft_freqgrandaverage([],PD2_bsPow{:});
ctrl2_avg = ft_freqgrandaverage([],ctrl2_bsPow{:});

cfg = [];
% cfg.channel = channels;
% cfg.avgoverchan = 'yes';
cfg.frequency = [14 25];
cfg.avgoverfreq = 'yes';
cfg.avgovertime = 'no';

BdataPD1 = ft_selectdata(cfg,PD1_avg);
BtracePD1 = (squeeze(BdataPD1.powspctrm));

% BtracePD1bsln = BtracePD1-nanmean(BtracePD1(bline));
BdataPD2 = ft_selectdata(cfg,PD2_avg);
BtracePD2 = squeeze(BdataPD2.powspctrm);
% BtracePD2bsln = BtracePD2-nanmean(BtracePD2(bline));

BdataCtrl1 = ft_selectdata(cfg,ctrl1_avg);
BtraceCtrl1 = (squeeze(BdataCtrl1.powspctrm));
% BtraceCtrl1bsln = BtraceCtrl1-nanmean(BtraceCtrl1(bline));
BdataCtrl2 = ft_selectdata(cfg,ctrl2_avg);
BtraceCtrl2 = squeeze(BdataCtrl2.powspctrm);
% BtraceCtrl2bsln = BtraceCtrl2-nanmean(BtraceCtrl2(bline));

% Beta w/baseline
% figure;
% plot(timeDim,squeeze(BtracePD1bsln),'b-'); hold on
% plot(timeDim,squeeze(BtracePD2bsln),'b--');
% plot(timeDim,squeeze(BtraceCtrl1bsln),'r-');
% plot(timeDim,squeeze(BtraceCtrl2bsln),'r--');
% axis([-.5, 2.5, -inf, inf])
% line([0 0],[-9e-22, 9e-22],'color','k','LineStyle','--')
% xlabel('Time (s)');
% ylabel('Beta power')
% title('Beta-band power');
% legend('PD (off)','PD (on)', 'Ctrl (1st)','Ctrl (2nd)')
% savefig('/home/mikkel/PD_motor/rebound/export/Btraces.fig');
% print('/home/mikkel/PD_motor/rebound/export/Btraces.png', '-dpng')
% close

% Beta wo/baseline
figure;
plot(timeDim,squeeze(BtracePD1),'b-'); hold on
plot(timeDim,squeeze(BtracePD2),'b--');
plot(timeDim,squeeze(BtraceCtrl1),'r-');
plot(timeDim,squeeze(BtraceCtrl2),'r--');
axis([-.5, 2.5, -inf, inf])
line([0 0],[-inf, inf],'color','k','LineStyle','--')
xlabel('Time (s)'); ylabel('Beta power');
title('Beta-band power');
legend('PD (off)','PD (on)', 'Ctrl (1st)','Ctrl (2nd)')
savefig('/home/mikkel/PD_motor/rebound/export/Btraces_noBsln.fig');
print('/home/mikkel/PD_motor/rebound/export/Btraces_noBsln.png', '-dpng')

%% Topoplot

figure;
cfg = [];
cfg.baseline = [-inf inf];
cfg.baselinetype = 'relative';
cfg.layout = 'neuromag306cmb.lay';
cfg.comment = 'no';
cfg.marker = 'off';

cfg.zlim = [.5 1.5];
cfg.xlim = [.15 .6];
cfg.ylim = [12 25];

subplot(2,2,1); ft_topoplotTFR(cfg,PD_avg1); title('Beta desync. (PD)')
subplot(2,2,2); ft_topoplotTFR(cfg,ctrl_avg1); title('Beta desync. (ctrl.)')

cfg.zlim = [0.5 1.5];
cfg.xlim = [.75 1.20];
cfg.ylim = [12 25];

subplot(2,2,3); ft_topoplotTFR(cfg,PD_avg2); title('Beta rebound (PD)')
subplot(2,2,4); ft_topoplotTFR(cfg,ctrl_avg2); title('Beta rebound (ctrl.)')

savefig('/home/mikkel/PD_motor/rebound/export/topoplots.fig');
print('/home/mikkel/PD_motor/rebound/export/topoplots.png', '-dpng')
% close
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 1,'\bf Topographies','HorizontalAlignment' ,'center','VerticalAlignment', 'top')



%% EXIT
% exit