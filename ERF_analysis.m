%%%%% GROUP SUMMARY ERF FOR PD-PROJ: BETA PART 
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

%% Load timelockeds -  combine into one structure, flip if necessary and baseline correct
ERF_PD1 = cell(1,length(PD_subs));
ERF_PD2 = cell(1,length(PD_subs));
ERF_ctrl1 = cell(1,length(ctrl_subs));
ERF_ctrl2 = cell(1,length(ctrl_subs));

for ii = 1:length(PD_subs)
    subID = PD_subs{ii};
    disp(subID);
    sub_dir = [dirs.megDir,'/',subID];
    
    %Session 1
    load([sub_dir,'/',subID,'_1-evoked.mat']);
    timelocked = ft_combineplanar([],timelocked);
    if any(strcmp(lh_subs,subID))
        disp('Flip')
        timelocked = flip_sens_neuromag(timelocked);
    end
    ERF_PD1{ii} = timelocked;
    
    %Session 2
    load([sub_dir,'/',subID,'_2-evoked.mat']);
    timelocked = ft_combineplanar([],timelocked);
    if any(strcmp(lh_subs,subID))
        disp('Flip')
        timelocked = flip_sens_neuromag(timelocked);
    end
    ERF_PD2{ii} = timelocked;
        
end

for ii = 1:length(ctrl_subs)
    subID = ctrl_subs{ii};
    disp(subID);
    sub_dir = [dirs.megDir,'/',subID];
    % Session 1
    load([sub_dir,'/',subID,'_1-evoked.mat']);
    timelocked = ft_combineplanar([],timelocked);
    ERF_ctrl1{ii} = timelocked;
    % Session 2
    load([sub_dir,'/',subID,'_2-evoked.mat']);
    timelocked = ft_combineplanar([],timelocked);
    ERF_ctrl2{ii} = timelocked;
end

save([dirs.output,'/all_ERF.mat'],'ERF_ctrl1','ERF_ctrl2','ERF_PD1','ERF_PD2');
disp('done')

%% Grand Average

cfg = [];
cfg.keepindividual = 'no';

avgERF.PD1 = ft_timelockgrandaverage(cfg, ERF_PD1{:});
avgERF.PD2 = ft_timelockgrandaverage(cfg, ERF_PD2{:});
avgERF.ctrl1 = ft_timelockgrandaverage(cfg, ERF_ctrl1{:});
avgERF.ctrl2 = ft_timelockgrandaverage(cfg, ERF_ctrl2{:});

% Save data
save([dirs.output,'/grandERF.mat'],'avgERF','-v7.3');
disp('done')

%% Plot grand avg.
cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.xlim = [-.1 1];
ft_multiplotER(cfg,avgERF.PD1,avgERF.ctrl1,avgERF.PD2,avgERF.ctrl2);

%% Plot individual
cfg = [];
cfg.layout = 'neuromag306cmb.lay';
ft_multiplotER(cfg,ERF_PD1{:});

%% Combined avarge across sessions (to find peak)
ERFpool = cell(1,length(subs));

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 90;

for ii = 1:length(subs)
    subID = subs{ii};
    sub_dir = [dirs.megDir,'/',subID];
    loadname1 = [dirs.megDir,'/',subID,'/',subID,'_1-epochs.mat'];
    load(loadname1);
    data1 = ft_preprocessing(cfg, data);
    loadname2 = [dirs.megDir,'/',subID,'/',subID,'_2-epochs.mat'];
    load(loadname1);
    data2 = ft_preprocessing(cfg, data);
    
    pooldata = ft_appenddata([],data1,data2);
    disp('Data appended...');
    
    poolERF = ft_timelockanalysis([],pooldata);
    ERFpool{ii} = ft_combineplanar([],poolERF);
    disp(['done with ',subID])    
end

save([dirs.output,'/ERFpooled.mat'],'ERFpool','-v7.3');
disp('done');



