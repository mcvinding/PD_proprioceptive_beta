%%%%% Get as summary of useful trials after epoching. 

addpath /home/mikkel/PD_motor/global_scripts
[dirs, sub_info, lh_subs] = PD_proj_setup('rebound');

cd(dirs.megDir);

subs = dir(dirs.megDir);                                %Find subjects in folder
subs = {subs([subs.isdir]).name};                       %Make list
subs = subs(~(strcmp('.',subs)|strcmp('..',subs)));     %Remove dots

%% Get n trials
trial_summary1 = zeros(length(subs),1);
trial_summary2 = zeros(length(subs),1);

for ii = 1:length(subs)
    subID = subs{ii};
    sub_dir = [dirs.megDir,'/',subID];
    loadname1 = [dirs.megDir,'/',subID,'/',subID,'_1-epochs.mat'];
    loadname2 = [dirs.megDir,'/',subID,'/',subID,'_2-epochs.mat'];
    load(loadname1);
    trial_summary1(ii) = length(data.trial);
    clear data
    load(loadname2);
    trial_summary2(ii) = length(data.trial);
    
end

save([dirs.output,'/trial_summary.mat'],'trial_summary1','trial_summary2');

csvwrite([dirs.output,'/ntrial1.csv'],trial_summary1)
csvwrite([dirs.output,'/ntrial2.csv'],trial_summary2)

disp('done');
