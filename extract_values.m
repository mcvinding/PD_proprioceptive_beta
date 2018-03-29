%%%%% Extract values for correlation analysis %%%%%

addpath /home/mikkel/PD_motor/global_scripts
[dirs, sub_info, lh_subs] = PD_proj_setup('rebound');

cd(dirs.megDir);

subs = dir(dirs.megDir);                             %Find subjects in folder
subs = {subs([subs.isdir]).name};                    %Make list
subs = subs(~(strcmp('.',subs)|strcmp('..',subs)));  %Remove dots

badsubs = {'0393'}; % Too few trials/bad data
subs(strcmp(badsubs,subs)) = [];

PD_subs = intersect(sub_info.PD,subs);
ctrl_subs = intersect(sub_info.ctrl,subs);

%% Get baseline power
% load([dirs.output,'/betaChan1log.mat']);
load([dirs.output,'/all_TFR.mat']);     %'data_ctrl1','data_ctrl2','data_PD1','data_PD2');
load([dirs.output,'/peak_chans.mat']);  %'channam_ctrl','channam_PD');
disp('done loading')

% Geat mean baseline absolute power
PD1_bsPow = zeros(1,length(data_PD1));
PD2_bsPow = zeros(1,length(data_PD2));
ctrl1_bsPow = zeros(1,length(data_ctrl1));
ctrl2_bsPow = zeros(1,length(data_ctrl2));

cfg = [];
cfg.frequency       = [14 25];
cfg.avgoverfreq     = 'yes';
cfg.latency         = [-1.25 -.2];
cfg.avgovertime     = 'yes';

for ii = 1:length(data_PD1)
    cfg.channel = channam_PD(ii);
    tempDat = ft_selectdata(cfg, data_PD1{ii});
    PD1_bsPow(ii) = tempDat.powspctrm;
%     PD1_bsPow{ii}.label = {'peak_channel'};
    
    tempDat = ft_selectdata(cfg, data_PD2{ii});
    PD2_bsPow(ii) = tempDat.powspctrm;
%     PD2_bsPow{ii}.label = {'peak_channel'};
end

for ii = 1:length(data_ctrl1)
    cfg.channel = channam_ctrl(ii);

    tempDat = ft_selectdata(cfg, data_ctrl1{ii});
    ctrl1_bsPow(ii) = tempDat.powspctrm;
%     ctrl1_bsPow{ii}.label = {'peak_channel'};

    tempDat = ft_selectdata(cfg, data_ctrl2{ii});
    ctrl2_bsPow(ii) = tempDat.powspctrm;
%     ctrl2_bsPow{ii}.label = {'peak_channel'};
    
end
disp('done')

%% movement related beta ERD

PD1_avgERD = zeros(1,length(PD_beta1bs));
PD2_avgERD = zeros(1,length(PD_beta2bs));
ctrl1_avgERD = zeros(1,length(PD_beta1bs));
ctrl2_avgERD = zeros(1,length(PD_beta2bs));

cfg = [];
cfg.frequency       = [14 25];
cfg.avgoverfreq     = 'yes';
cfg.latency         = [.1 .5];
cfg.avgovertime     = 'yes';

for ii = 1:length(PD_beta1bs)
%     cfg.channel = channam_PD(ii);
    tempDat = ft_selectdata(cfg, PD_beta1bs{ii});
    PD1_avgERD(ii) = tempDat.powspctrm;
    
    tempDat = ft_selectdata(cfg, PD_beta2bs{ii});
    PD2_avgERD(ii) = tempDat.powspctrm;
end

for ii = 1:length(data_ctrl1)
%     cfg.channel = channam_ctrl(ii);

    tempDat = ft_selectdata(cfg, ctrl_beta1bs{ii});
    ctrl1_avgERD(ii) = tempDat.powspctrm;

    tempDat = ft_selectdata(cfg, ctrl_beta2bs{ii});
    ctrl2_avgERD(ii) = tempDat.powspctrm;
end

disp('done')

%% Combine into one and save
PD_id = str2num(cell2mat(PD_subs'));
ctrl_id = str2num(cell2mat(ctrl_subs'));
all_id = [PD_id; PD_id; ctrl_id; ctrl_id];
all_bsPow = [PD1_bsPow PD2_bsPow ctrl1_bsPow ctrl2_bsPow]*10^20; %Change scale to avoid rounding error when saving
all_avgERD = [PD1_avgERD PD2_avgERD ctrl1_avgERD ctrl2_avgERD];
all_avgERS = [PD1_avgRebnd PD2_avgRebnd ctrl1_avgRebnd ctrl2_avgRebnd];
all_group = [ones(1,2*length(PD_id)), 2*ones(1,2*length(ctrl_id))];
all_session = [ones(1,length(PD_id)), 2*ones(1,length(PD_id)), ones(1,length(ctrl_id)), 2*ones(1,length(ctrl_id))];

M = [all_id'; all_group; all_session; all_bsPow; all_avgERD; all_avgERS]';

dlmwrite('summary_pow.txt',M,'delimiter',';')