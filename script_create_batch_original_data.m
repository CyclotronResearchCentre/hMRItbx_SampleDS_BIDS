%% SCRIPT generating the batch for the hMRI data processing
%
% The batch needs the correct path to the data but these are stored on the
% user's computer. Therefore this script simply
% 1/ tries to figure out where the data are stored or asks the user to 
%    specify the main folder manually.
% 2/ then it sets up the hMRI map creation batch, with the selection of the
%    provided 'hmri_local_defaults.m' file with default parameters.
% 3/ finally saves it a .mat file for further use, e.g. loading in the
%    batch-GUI
%__________________________________________________________________________
%
% 1/ Finding or defining the 'rooDir' of the data folder.
% The script successively tries a few options:
% a) the user explicitly sets up the folder path
% b) check if the data are in the current working directory
% c) check if the data can be found relative to this script's current
%    folder (one folder up actually)
% d) if all alse fails, ask the user to select the folder manually
%
% 2/ Create the hMRI batch
% The file sub-folder names are all explicitly stated, as they're those
% from the original data set.

%% 1/ Finding or defining the 'rooDir' of the data folder.
% The script successively tries a few options:
% a) the user explicitly sets up the folder path
% b) check if the data are in the current working directory
% c) check if the data can be found relative to this script's current
%    folder (one folder up actually)
% d) if all alse fails, ask the user to select the folder manually

rD_flag = false; % Meaning rootDir is not sorted yet

%%% ======================================================================
% SET UP THE ROOT FOLDER HERE
% a) manual setting of the rootDir by editing this file.
rootDir = '';
% -> MUST BE UPDATED by each user. For example:
% rootDir = 'C:\Dox\2_Data\hmri_sample_dataset_with_maps';
%%% ======================================================================

% Do the check, using the 1st echo MT image
fn_MT1 = fullfile(rootDir,'mtw_mfc_3dflash_v1i_R4_0012', ...
    'anon_s2018-02-28_18-26-190132-00001-00224-1.nii');
if exist(fn_MT1,'file')
    % sorted
    rD_flag = true;
else
    if isempty(rootDir), rootDir = '<EMPTY>'; end
    fprintf('\n''rootDir'' (%s) not correct.\n',rootDir)
end

% b) check current working directory
if ~rD_flag
    rootDir = pwd;
    % Do the check, using the 1st echo MT image
    fn_MT1 = fullfile(rootDir,'mtw_mfc_3dflash_v1i_R4_0012', ...
        'anon_s2018-02-28_18-26-190132-00001-00224-1.nii');
    if exist(fn_MT1,'file')
        % sorted
        rD_flag = true;
    else
        fprintf('''rootDir'' (%s) not correct.\n',rootDir)
    end
end

% c) check directory above the one containing this file
if ~rD_flag
    rootDir = spm_file(spm_file(mfilename('fullpath'),'path'),'fpath');
    % Do the check, using the 1st echo MT image
    fn_MT1 = fullfile(rootDir,'mtw_mfc_3dflash_v1i_R4_0012', ...
        'anon_s2018-02-28_18-26-190132-00001-00224-1.nii');
    if exist(fn_MT1,'file')
        % sorted
        rD_flag = true;
    else
        fprintf('''rootDir'' (%s) not correct.\n',rootDir)
    end   
end

% d) ask the the user to select the correct folder
while ~rD_flag
    rootDir = spm_select(1,'dir','Select main directory with hMRI data');
    % Do the check, using the 1st echo MT image
    fn_MT1 = fullfile(rootDir,'mtw_mfc_3dflash_v1i_R4_0012', ...
        'anon_s2018-02-28_18-26-190132-00001-00224-1.nii');
    if exist(fn_MT1,'file')
        % sorted
        rD_flag = true;
    else
        fprintf('''rootDir'' (%s) not correct.\n',rootDir)
    end   
end

% Display folder that will be used
fprintf('\n Using ''rootDir'' (%s).\n',rootDir)

%% 2/ Create the hMRI batch
% The file sub-folder names are all explicitly stated, as they're those
% from the original data set.

clear matlabbatch
% First with the loading of the hmri_local_defaults
matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {
    fullfile(rootDir,'Batch_Script_hMRI_Config','hmri_local_defaults.m')};
% Then the map creation with all the field maps
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.output.indir = 'yes';
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_per_contrast.raw_sens_MT = {
    fullfile(rootDir,'mfc_smaps_v1a_Array_0010','anon_s2018-02-28_18-26-190114-00001-00022-1.nii,1')
    fullfile(rootDir,'mfc_smaps_v1a_QBC_0011','anon_s2018-02-28_18-26-190121-00001-00022-1.nii,1')
    };
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_per_contrast.raw_sens_PD = {
    fullfile(rootDir,'mfc_smaps_v1a_Array_0007','anon_s2018-02-28_18-26-185318-00001-00022-1.nii,1')
    fullfile(rootDir,'mfc_smaps_v1a_QBC_0008','anon_s2018-02-28_18-26-185329-00001-00022-1.nii,1')
    };
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_per_contrast.raw_sens_T1 = {
    fullfile(rootDir,'mfc_smaps_v1a_Array_0013','anon_s2018-02-28_18-26-190901-00001-00022-1.nii,1')
    fullfile(rootDir,'mfc_smaps_v1a_QBC_0014','anon_s2018-02-28_18-26-190907-00001-00022-1.nii,1')
    };
%
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.b1input = {
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184837-00001-00001-1.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184837-00001-00049-2.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184849-00002-00097-1.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184849-00002-00145-2.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184901-00003-00193-1.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184901-00003-00241-2.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184913-00004-00289-1.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184913-00004-00337-2.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184925-00005-00385-1.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184925-00005-00433-2.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184937-00006-00481-1.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184937-00006-00529-2.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184949-00007-00577-1.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-184949-00007-00625-2.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-185001-00008-00673-1.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-185001-00008-00721-2.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-185013-00009-00769-1.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-185013-00009-00817-2.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-185025-00010-00865-1.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-185025-00010-00913-2.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-185037-00011-00961-1.nii,1')
    fullfile(rootDir,'mfc_seste_b1map_v1e_0004','anon_s2018-02-28_18-26-185037-00011-01009-2.nii,1')
    };
%
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.b0input = {
    fullfile(rootDir,'gre_field_mapping_1acq_rl_0005','anon_s2018-02-28_18-26-185100-00001-00001-1.nii,1')
    fullfile(rootDir,'gre_field_mapping_1acq_rl_0005','anon_s2018-02-28_18-26-185101-00001-00001-2.nii,1')
    fullfile(rootDir,'gre_field_mapping_1acq_rl_0006','anon_s2018-02-28_18-26-185101-00001-00001-2.nii,1')
    };
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.b1parameters.b1metadata = 'yes';
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT = {
    fullfile(rootDir,'mtw_mfc_3dflash_v1i_R4_0012','anon_s2018-02-28_18-26-190132-00001-00224-1.nii,1')
    fullfile(rootDir,'mtw_mfc_3dflash_v1i_R4_0012','anon_s2018-02-28_18-26-190132-00001-00448-2.nii,1')
    fullfile(rootDir,'mtw_mfc_3dflash_v1i_R4_0012','anon_s2018-02-28_18-26-190132-00001-00672-3.nii,1')
    fullfile(rootDir,'mtw_mfc_3dflash_v1i_R4_0012','anon_s2018-02-28_18-26-190132-00001-00896-4.nii,1')
    fullfile(rootDir,'mtw_mfc_3dflash_v1i_R4_0012','anon_s2018-02-28_18-26-190132-00001-01120-5.nii,1')
    fullfile(rootDir,'mtw_mfc_3dflash_v1i_R4_0012','anon_s2018-02-28_18-26-190132-00001-01344-6.nii,1')
    };
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD = {
    fullfile(rootDir,'pdw_mfc_3dflash_v1i_R4_0009','anon_s2018-02-28_18-26-185345-00001-00224-1.nii,1')
    fullfile(rootDir,'pdw_mfc_3dflash_v1i_R4_0009','anon_s2018-02-28_18-26-185345-00001-00448-2.nii,1')
    fullfile(rootDir,'pdw_mfc_3dflash_v1i_R4_0009','anon_s2018-02-28_18-26-185345-00001-00672-3.nii,1')
    fullfile(rootDir,'pdw_mfc_3dflash_v1i_R4_0009','anon_s2018-02-28_18-26-185345-00001-00896-4.nii,1')
    fullfile(rootDir,'pdw_mfc_3dflash_v1i_R4_0009','anon_s2018-02-28_18-26-185345-00001-01120-5.nii,1')
    fullfile(rootDir,'pdw_mfc_3dflash_v1i_R4_0009','anon_s2018-02-28_18-26-185345-00001-01344-6.nii,1')
    fullfile(rootDir,'pdw_mfc_3dflash_v1i_R4_0009','anon_s2018-02-28_18-26-185345-00001-01568-7.nii,1')
    fullfile(rootDir,'pdw_mfc_3dflash_v1i_R4_0009','anon_s2018-02-28_18-26-185345-00001-01792-8.nii,1')
    };
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1 = {
    fullfile(rootDir,'t1w_mfc_3dflash_v1i_R4_0015','anon_s2018-02-28_18-26-190921-00001-00224-1.nii,1')
    fullfile(rootDir,'t1w_mfc_3dflash_v1i_R4_0015','anon_s2018-02-28_18-26-190921-00001-00448-2.nii,1')
    fullfile(rootDir,'t1w_mfc_3dflash_v1i_R4_0015','anon_s2018-02-28_18-26-190921-00001-00672-3.nii,1')
    fullfile(rootDir,'t1w_mfc_3dflash_v1i_R4_0015','anon_s2018-02-28_18-26-190921-00001-00896-4.nii,1')
    fullfile(rootDir,'t1w_mfc_3dflash_v1i_R4_0015','anon_s2018-02-28_18-26-190921-00001-01120-5.nii,1')
    fullfile(rootDir,'t1w_mfc_3dflash_v1i_R4_0015','anon_s2018-02-28_18-26-190921-00001-01344-6.nii,1')
    fullfile(rootDir,'t1w_mfc_3dflash_v1i_R4_0015','anon_s2018-02-28_18-26-190921-00001-01568-7.nii,1')
    fullfile(rootDir,'t1w_mfc_3dflash_v1i_R4_0015','anon_s2018-02-28_18-26-190921-00001-01792-8.nii,1')
    };
matlabbatch{2}.spm.tools.hmri.create_mpm.subj.popup = true;

%% 3/ save the resulting matlabbatch in a .mat file
fn_mat = fullfile(rootDir,'Batch_Script_hMRI_Config','my_Create_Maps_Batch.mat');
save(fn_mat,'matlabbatch')
