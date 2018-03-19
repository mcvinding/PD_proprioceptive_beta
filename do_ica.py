# -*- coding: utf-8 -*-
"""
Updated Apr 26 2017

"""
import matplotlib
matplotlib.use('Agg')
from mne.preprocessing import ICA, create_ecg_epochs, create_eog_epochs, read_ica
from mne.io import Raw
from mne import pick_types
from os import listdir, chdir, mkdir, path
import numpy as np
import sys
sys.path.append('/home/mikkel/PD_motor/global_scripts')
from PD_motor_pyfun import run_ICA_wMNE

#import gc
#gc.enable()
#gc.DEBUG_COLLECTABLE()

# %% run options
project_part = 'rebound'
overwrite_old_files = True      # Wheter files should be overwritten if already exist

# Paths
raw_path = '/home/mikkel/PD_motor/'+project_part+'/raw' # Files sorted in raw folder by means of symbolic links
meg_path = '/home/mikkel/PD_motor/'+project_part+'/meg_data' # Work directory for MEG

sessions = '1', '2'

chdir(raw_path)
file_list = listdir(raw_path)

# %% Run through files
allReadyRunFiles = []
for ff in file_list:
    sub = ff[:4]

    #Find all files for same subject to run ICA on.  
    print('Now reading file: '+ff)
    prefx = ff[:8]       
 
    if ff in allReadyRunFiles:
        print('...Nope! I already have done that! Moving to next file.')
        continue
    else:
        allReadyRunFiles += [f for f in file_list if f.startswith(sub)]
        inFiles = [raw_path+'/'+f for f in file_list if f.startswith(sub)]
        print('Files to read = '+str(len(inFiles)))
        
    print 'Now running filter+ICA for sub = ' + sub   
    
    # % Make dirs for output
    out_path = meg_path+'/'+sub
#    out_icaFname = project_part+'-ica_comp.fif' # FOr now we run with the old name structure[!]
    out_icaFname = 'comp-ica.fif'
    
    if not path.exists(out_path):
        mkdir(out_path)
    if path.isfile(out_path+'/'+out_icaFname):
        print 'File "'+out_path+'/'+out_icaFname+'" already exists.'
        if not overwrite_old_files:
            print 'Do not overwrite. Pass.'
            continue
        elif overwrite_old_files:
            print 'Do overwrite! Running everything again!'

#% RUN ICA - save file in out folder
    chdir(out_path)
    
    if '0377' in sub:
        for kk, nn in enumerate(inFiles):
            run_ICA_wMNE(nn,out_path,ica_fname=str(kk)+'_comp-ica.fif')
    else:
        run_ICA_wMNE(inFiles,out_path)  #,ica_fname=out_fname)
    
    
#%% Find ECG and EOG artefacts                  
# Intit. variables
n_max_ecg = 3
n_max_eog = 2
title = 'Sources related to %s artifacts (red)'

for ff in file_list:
    print('Now reading file: '+ff)
    
    sub = ff[0:4]    
    session = ff[18]

    subDir = meg_path+'/'+sub
    in_fname = raw_path+'/'+ff
    out_fname = path.join(subDir,sub+'_'+session+'-ica_raw.fif')
    chdir(subDir)
    
    if path.isfile(out_fname) and not overwrite_old_files:
        print 'File '+out_fname+' aready exists. Skipping....'
        continue  
    
    #read ica from file
    ica_fname = subDir+'/comp-ica.fif'

    if not path.isfile(ica_fname):
        print('ICA file '+ica_fname+' does not exist')
        continue
    else:
        print('found ICA file. Proceeding...')
    
    ica = read_ica(ica_fname)
    ica.labels_ = dict()
    
    # %% Load data and remove components
    raw = Raw(in_fname, preload=True)
    
    picks_meg = pick_types(raw.info, meg=True, eeg=False, eog=False, emg=False, misc=False, 
                           stim=False, exclude='bads')
    picks_eXg = pick_types(raw.info, meg=False, eeg=False, eog=True, ecg = True, emg=False, misc=False, 
                           stim=False, exclude='bads')
                           
#DO not filter in this part. TFR later.                           
    raw.filter(1, 40, n_jobs=3, picks=picks_eXg)
    raw.notch_filter(50, n_jobs=3, picks=picks_eXg)
    
    ecg_epochs = create_ecg_epochs(raw, ch_name='ECG003', tmin=-.5, tmax=.5)    #, picks=picks)
    ecg_inds, ecg_scores = ica.find_bads_ecg(ecg_epochs, method='ctps', verbose=False)
    
    ecg_scores_fig = ica.plot_scores(ecg_scores, exclude=ecg_inds, title=title % 'ecg')
    ecg_scores_fig.savefig(sub+'_'+session+'_ICA_ecg_comp_score.png')

    if ecg_inds:
        show_picks = np.abs(ecg_scores).argsort()[::-1][:5]
        
        ecg_source_fig = ica.plot_sources(raw, show_picks, exclude=ecg_inds,
                         title=title % 'ECG')
        ecg_source_fig.savefig(sub+'_'+session+'_ICA_ecg_source.png')

        ecg_comp_fig = ica.plot_components(ecg_inds, title=title % 'ecg', colorbar=True)
        ecg_comp_fig.savefig(sub+'_'+session+'_ICA_ecg_comp.png')
        
    # estimate average artifact
    ecg_evoked = ecg_epochs.average()
    # plot ECG sources + selection
    ecg_evo_fig1 = ica.plot_sources(ecg_evoked, exclude=ecg_inds)
    # plot ECG cleaning
    ecg_evo_fig2 = ica.plot_overlay(ecg_evoked, exclude=ecg_inds)
    
    ecg_evo_fig1.savefig(sub+'_'+session+'_ICA_ecg_latSource.png')
    ecg_evo_fig2.savefig(sub+'_'+session+'_ICA_ecg_cleaning.png')
    
    ecg_inds = ecg_inds[:n_max_ecg]
    ica.exclude += ecg_inds
    
    #%% Find EOG artifacts
    
    eog_inds, eog_scores = ica.find_bads_eog(raw)
    
    eog_scores_fig = ica.plot_scores(eog_scores, exclude=eog_inds, title=title % 'eog')
    eog_scores_fig.savefig(sub+'_'+session+'ICA_eog_comp_score.png')

    if eog_inds:
        show_picks = np.abs(eog_scores[0]).argsort()[::-1][:5]
#        eog_source_fig = ica.plot_sources(raw, show_picks, exclude=eog_inds,
#                         title=title % 'EOG',show=False)
#                         
#        show_picks = np.abs(eog_scores[1]).argsort()[::-1][:5]
#        ica.plot_sources(raw, show_picks, exclude=eog_inds,
#                         title=title % 'EOG')
                         
        eog_comp_fig = ica.plot_components(eog_inds, title="Sources related to EOG artifacts",
                                colorbar=True,show=False)
        eog_comp_fig.savefig(sub+'_'+session+'ICA_eog_comp.png')

    # estimate average artifact
    eog_evoked = create_eog_epochs(raw, tmin=-.75, tmax=.75, picks=picks_meg).average()
                                   
#    ica.plot_sources(eog_evoked, exclude=eog_inds)

    # plot EOG sources + selection
    eog_evo_fig = ica.plot_overlay(eog_evoked, exclude=eog_inds)  # plot EOG cleaning
    eog_evo_fig.savefig(sub+'_'+session+'ICA_eog_cleaning.png')
    # check the amplitudes do not change
    ica.plot_overlay(raw, show=False)  # EOG artifacts remain
        
    eog_inds = eog_inds[:n_max_eog]
    ica.exclude += eog_inds
    
    #%% Apply the solution to Raw, Epochs or Evoked like this:
    raw_ica = ica.apply(raw, copy=False)
    raw_ica.save(out_fname, overwrite=overwrite_old_files)
    
    print('----------- DONE (Sub: '+sub+' Session: '+session+') -----------------')
    
print('----------- ALL DONE -----------------')

