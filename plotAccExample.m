%% Plot example of acceprometer
addpath /home/mikkel/matlab/export_fig/
dirs.export = '/home/mikkel/PD_motor/rebound/export/publication';

data_file = '/home/mikkel/PD_motor/rebound/meg_data/0327/0327_1-ica_raw.fif';

hdr = ft_read_header(data_file);
misc_chan = find(~cellfun(@isempty, strfind(hdr.label, 'MISC')));

cfg = [];
cfg.dataset         = data_file;
cfg.continuous      = 'yes';
cfg.channel         = misc_chan(1:3);
cfg.bpfilter        = 'yes';
cfg.bpfreq          = [1 195];
alldataACC = ft_preprocessing(cfg);

temp_acc = alldataACC.trial{:};
euc_right = sqrt(sum(temp_acc.^2,1));

plot(euc_right(14000:32000), 'k', 'LineWidth', 1);
ylim([-0.1, 0.3])

set(gcf, 'Position', [100, 100, 1000, 150])
y2 = ylabel({'Acceleration'},'fontsize',6);
set(gca,'xtick',[],'xticklabel','','XColor','none','LineWidth', 1 )
set(gca,'ytick',[],'yticklabel','','YColor','none','LineWidth', 1 )
set(y2, 'HorizontalAlignment','center','VerticalAlignment','bottom')
title(' ');

cd(dirs.export);
export_fig('accExample.png', '-r500', '-p0.05', '-CMYK')
