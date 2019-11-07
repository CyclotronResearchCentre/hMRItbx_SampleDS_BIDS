%% Script to BIDSify the images from the hMRI example dataset
% Start with the anatomical images then the field and RF sensitivity maps!
%
% It is assumed that the data are organized as downloaded
% -> we know which specific folder contains what sort of data
%
% For sanity check parameters are available in the mpm_params.mat file
% previously created, so one can get the metadata directly from
% - the Dicom header (saved in the current .json file) or
% - via direct pulling from these mpm_params.mat as they're used in the
%   hMRI toolbox
%  This works similarly for the RFsens (B1-) and B1map (B0 and B1+) field
%  maps.
% 
% Still missing/to check:
% =======================
% - turn this into a function with a few option flags as to how things are
%   done, especially what should be saved (full DICOM header?) or not.
% - add the required top metadata files describing the dataset and 
%   participant
%   -> PARTLY DONE
% - for the 'fmap' images at the "IntendedFor" metadata
% - gzip all images ('gzip' function)?
% - check filenames convention                \_ according to BEP001 crowd
% - check metadata list (names and definition)/
%   -> PARTLY DONE
% - check the units of metadata extracted
%   -> LOOKS LIKE IT's SORTED
% - the '_mod-*' is used for the B1- maps. These are not defined in BEP001
%   so far. Moreover the filename '_mod-*' field is currently ONLY defined 
%   for the case of '_defacemask.nii' images!
%   -> COMBINE THINGS INTO THE '_acq-*' field
% - the B1+ maps is defined but no B1- map in the BEP001 so far
% - for B1+ (and maybe B1-) maps, there is a specific suffix, B1plusmap
%   (and maybe B1minusmap) but not for the B0 images.
%   -> ALREADY DISCUSSED, KEEP B0 STUFF AS IT IS FOR BACK COMPATIBILITY
% 
% Dependencies:
% =============
% - a function called "hmri_BIDSify_json.m" to refactor the metadata from 
%   the DICOM header. This itself relies on the get_metadata_val.m function
%   from the hMRI toolbox
% - a tab-separated-value file, JSONtabl_dcm2bids.tsv, with the list of 
%   metadata fields required.
% These 2 are already on the hMRI private github server in Leipzig
% -> need to keep track of these.
%__________________________________________________________________________
% Copyright (C) 2019 GIGA Institute

% Written by C. Phillips, 2019.
% Cyclotron Research Centre, University of Liege, Belgium

% NOTES:
% 1/ It could be useful to create a batch-GUI to select (semi-)manually the
% data in order to figure out what is really at hand
% -> one can reuse the "get_*_params.m" functions
% Sub-functions:
% - get_mpm_params.m in hmri_create_MTProt.m -> qMRI creation
% - get_rfsens_params.m in hmri_create_RFsens.m -> B1- (RF sensitivity)
% - get_b1map_params.m in hmri_create_b1map -> B0/B1+ field maps
%
% 2/ The RF sensitivity correction step DOES modify the whole data
% structure as available in the 'job' variable. This is a bad idea but this
% also makes the filenames available in the mpm_params.mat structure
% unusable since the "intermediate" files listed do not exist anymore at
% the end of the map creation process...

%% Some defaults

% Main folder with the data
rootPth = 'C:\Dox\2_Data\hmri_sample_dataset_with_maps';
% Subject label
subj_label = 'anon';
% Where the BIDS-ified data will be copied, bidsPth & subjPth
bidsPth = fullfile(rootPth,'BIDS_dataset_v3');
if ~exist(bidsPth,'dir'), mkdir(bidsPth), end;
subjPth = fullfile(bidsPth,sprintf('subj-%s',subj_label));
if ~exist(subjPth,'dir'), mkdir(subjPth), end;

%% Deal with the anatomical images

% 1/ Get the filenames for the 3 types (MTw, PDw, T1w) of anatomical images
% -------------------------------------------------------------------------
seq_label = { 'MTw' , 'PDw' , 'T1w' };
pth_MPM = { ...
    'mtw_mfc_3dflash_v1i_R4_0012' , ...
    'pdw_mfc_3dflash_v1i_R4_0009' , ...
    't1w_mfc_3dflash_v1i_R4_0015' };

fn_MPM = cell(1,3) ; % empty cell array for filenames
% get MTw
pth_MTw = fullfile(rootPth,pth_MPM{1});
fn_MPM{1} = spm_select('FPList',pth_MTw,'^.*\.nii$');
% get PDw
pth_PDw = fullfile(rootPth,pth_MPM{2});
fn_MPM{2} = spm_select('FPList',pth_PDw,'^.*\.nii$');
% get T1w
pth_T1w = fullfile(rootPth,pth_MPM{3});
fn_MPM{3} = spm_select('FPList',pth_T1w,'^.*\.nii$');

% 2/ Move and rename anatomical images + update corresponding json file
% ---------------------------------------------------------------------
% Filename structure
fn_bids_struct = 'sub-%s_echo-%d_acq-%s_MPM.%s';
% feed in "subject label", "echo index", "acqusition sequence", "extension"
% According to BEP001, multi-echi anatomical images acquired with the MPM
% protocal should be of type '_MPM'. The contrast (MT/PD/T1-weighted) is
% only an acquisition parameter -> in the '_acq-***_' part of the filename

% 'anat' folder: define & create as needed
anatPth = fullfile(subjPth,'anat');
if ~exist(anatPth,'dir'), mkdir(anatPth), end;

% MPM files, .nii and .json
for iseq = 1:3 % 3 types of sequences (MTw, PDw, T1w) *in that order*!
    nMPMw_ii = size(fn_MPM{iseq},1);
    if nMPMw_ii
        for ii=1:nMPMw_ii
            % a) deal with image files
            fn_bidsMPMw_ii = fullfile(anatPth, ...  % path
                sprintf( fn_bids_struct, ...            % filename
                subj_label, ii, seq_label{iseq}, 'nii' ) ); % feeds
            % copy file with name change
            copyfile(fn_MPM{iseq}(ii,:),fn_bidsMPMw_ii);
            % b) deal with JSON files
            % -> copy with right name then convert them all
            fn_json_in = spm_file(fn_MPM{iseq}(ii,:),'ext','json');
            fn_json_out = spm_file(fn_bidsMPMw_ii,'ext','json');
            copyfile(fn_json_in,fn_json_out)
        end
    else
        fprintf('\nNo %s filed.\n',seq_label{iseq});
    end
end
% Now update all JSON files
fn_jsonAnat = spm_select('FPList',anatPth,'^.*\.json$');
% In one go, using a single function to pull out all the potential fields
% necessary for the MPM data
fn_jsonAnat_BIDS = hmri_BIDSify_json(fn_jsonAnat);

% % or using a manually specified list and pre-extracted params
% fn_mpm_params = fullfile(rootPth,'MPM_params.mat');
% load(fn_mpm_params)
%
% for iseq = 1:3 % 3 types of sequences (MTw, PDw, T1w)
%     % find those .json files
%     filt_iseq = ['^.*_',seq_label{iseq},'\.json$'];
%     fn_jsonAnat_iseq = spm_select('FPList', anatPth, filt_iseq);
%     % check if there are such files
%     nJsonAnat_iseq = size(fn_jsonAnat_iseq,1);
%     if nJsonAnat_iseq
%         % get input parameters, account for sequence index!
%         params_iseq = mpm_params.input(mpm_params.([seq_label{iseq},'idx']));
%         for ii=1:nJsonAnat_iseq
%             fn_jsonAnat_ii = deblank(fn_jsonAnat(ii,:));
%             Jstr_in = spm_load(deblank(fn_jsonAnat(ii,:)));
%             nFields_in = numel(fieldnames(Jstr_in));
%             % add extrafields: EchoTime, FlipAngle,
%             % RepatitionTimeExcitiation, PulseSequenceType
%             Jstr_in.EchoTime = params_iseq.TE(ii);
%             Jstr_in.FlipAngle = params_iseq.fa;
%             Jstr_in.RepetitionTimeExcitiation = params_iseq.TR;
%             Jstr_in.PulseSequenceType = 'Flash';
%             % reorder fields -> put new fields first
%             nFields_out = numel(fieldnames(Jstr_in));
%             Jstr_in = orderfields(Jstr_in, ...
%                 [nFields_in+1:nFields_out 1:nFields_in]');
%             % save
%             spm_save(fn_jsonAnat_ii,Jstr_in, struct('indent','\t'))
%         end
%     end
% end

%% Deal with the field maps
% Work with the 3 different types of maps one at a time:
% - B1minus = RF sensitivity, 2 images per anatomical image type
%   -> modality (mod-MTw/PDw/T1w), acquisition (acq-head/body)
% - B1plus = emission with 22 images
%   -> echo index
% - B0 = classic field maps. 2 magnitude + 1 phase diff
%   -> part-magnitude1/magnitude2/phasediff

% 'fmap' folder: define & create as needed
fmapPth = fullfile(subjPth,'fmap');
if ~exist(fmapPth,'dir'), mkdir(fmapPth), end;

% 1/ Deal with 6 B1- maps(RF sensitivity), 3 mod-type x 2 acq-coils
% -----------------------------------------------------------------
% Filename structure
fn_bids_B1m_struct = 'sub-%s_acq-%s_mod-%s_B1minusmap.nii';
% feed in "subject label", "acquisition type" (head/body), "modality" (MTw/PDw/T1w)

% Filenames for the 3X2 types: (MTw, PDw, T1w) x (head/body)
mod_label = { 'MTw' , 'PDw' , 'T1w' };
acq_label = {'head', 'body'};
% and the corresponding folder in the data set structure
pth_B1m = { ...
    'mfc_smaps_v1a_Array_0010' , 'mfc_smaps_v1a_QBC_0011' ;...
    'mfc_smaps_v1a_Array_0007' , 'mfc_smaps_v1a_QBC_0008' ;...
    'mfc_smaps_v1a_Array_0013' , 'mfc_smaps_v1a_QBC_0014' };

for imod = 1:3
    for iacq = 1:2
        pth_iSA = fullfile(rootPth,pth_B1m{imod,iacq});
        % a) Deal with images
        % Select the image file requested.
        fn_iSA_nii = spm_select('FPList',pth_iSA,'^.*\.nii$');
        fn_bidsB1m_iAS_nii = fullfile(fmapPth, ...  % path
            sprintf( fn_bids_B1m_struct, ...            % filename
            subj_label, acq_label{iacq}, mod_label{imod} ) ); % feeds
        % copy file with name change
        copyfile(fn_iSA_nii,fn_bidsB1m_iAS_nii);
        % b) Deal with JSON
        fn_iSA_json = spm_file(fn_iSA_nii,'ext','json');
        fn_bidsB1m_iAS_json = spm_file(fn_bidsB1m_iAS_nii,'ext','json');
        % copy file with name change
        copyfile(fn_iSA_json,fn_bidsB1m_iAS_json);
    end
end

% 2/ Deal with B1+ (B1 bias correction), 11 FA x 2 acquisition types
% ------------------------------------------------------------------
% Filename structure
fn_bids_B1p_struct = 'sub-%s_acq-%s_fa-%02d_B1plusmap.nii';
% feed in "subject label", "acquisition type" (SE/STE), "flip angle" (index)
acq_label = {'SE', 'STE'};

% Select the images for the SE/STE acquisition
pth_B1p = fullfile(rootPth,'mfc_seste_b1map_v1e_0004');
fn_B1plus{1} = spm_select('FPList',pth_B1p,'^.*-1\.nii$');
fn_B1plus{2} = spm_select('FPList',pth_B1p,'^.*-2\.nii$');

for ifa = 1:11
    for iacq = 1:2
        % a) Deal with images
        fn_iFA_nii = deblank(fn_B1plus{iacq}(ifa,:));
        fn_bidsB1p_iFA_nii = fullfile(fmapPth, ...    % path
            sprintf( fn_bids_B1p_struct, ...      % filename
            subj_label, acq_label{iacq}, ifa ) ); % feeds
        % copy file with name change
        copyfile(fn_iFA_nii,fn_bidsB1p_iFA_nii);
        % b) Deal with JSON
        fn_iFA_json = spm_file(fn_iFA_nii,'ext','json');
        fn_bidsB1p_iFA_json = spm_file(fn_bidsB1p_iFA_nii,'ext','json');
        % copy file with name change
        copyfile(fn_iFA_json,fn_bidsB1p_iFA_json);
    end
end

% 3/ Deal with B0, 2 magnitude + 1 phase difference image
% -------------------------------------------------------
% First the 2 magnitude images
% Filename structure
fn_bids_B0mag_struct = 'sub-%s_magnitude%d.nii';
% feed in "subject label", "magnitude index"
B0mag_pth = fullfile(rootPth,'gre_field_mapping_1acq_rl_0005');
fn_B0mag_nii = spm_select('FPList',B0mag_pth,'^.*\.nii$');

for ii=1:2
    % a) Deal with images
    fn_bidsB0ii_nii = fullfile(fmapPth, ...    % path
        sprintf( fn_bids_B0mag_struct, ...      % filename
        subj_label, ii ) );         % feeds
    % copy file with name change
    copyfile(deblank(fn_B0mag_nii(ii,:)),fn_bidsB0ii_nii);
    % b) Deal with JSON
    fn_B0mag_json = spm_file(fn_B0mag_nii(ii,:),'ext','json');
    fn_bidsB0ii_json = spm_file(fn_bidsB0ii_nii,'ext','json');
    % copy file with name change
    copyfile(fn_B0mag_json,fn_bidsB0ii_json);
end

% Then the phase difference
% Filename structure
fn_bids_B0pdiff_struct = 'sub-%s_phasediff.nii';
B0pdiff_pth = fullfile(rootPth,'gre_field_mapping_1acq_rl_0006');
fn_B0pdiff_nii = spm_select('FPList',B0pdiff_pth,'^.*\.nii$');
% a) Deal with images
fn_bidsB0pdiff_nii = fullfile(fmapPth, ...    % path
    sprintf( fn_bids_B0pdiff_struct, subj_label ) );
% copy file with name change
copyfile(fn_B0pdiff_nii,fn_bidsB0pdiff_nii);
% b) Deal with JSON
fn_B0pdiff_json = spm_file(fn_B0pdiff_nii,'ext','json');
fn_bidsB0pdiff_json = spm_file(fn_bidsB0pdiff_nii,'ext','json');
% copy file with name change
copyfile(fn_B0pdiff_json,fn_bidsB0pdiff_json);

% Now update all JSON files
fn_jsonfmat = spm_select('FPList',fmapPth,'^.*\.json$');
% In one go, using a single function to pull out all the postential fields
% necessary for the MPM data
fn_jsonfmat_BIDS = hmri_BIDSify_json(fn_jsonfmat);
