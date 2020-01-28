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
% Operations
% ==========
% I. Deal with the anatomical images
% 1/ Get the filenames for the 3 types (MTw/PDw/T1w) of anatomical images
% 2/ Move and rename anatomical images
% 3/ update corresponding json file
%
% II. Deal with the field maps
% Work with the 3 different types of maps one at a time:
% 1/ Deal with 6 B1- maps(RF sensitivity), 3 mod-type x 2 acq-coils
%    B1minus = RF sensitivity, 2 images per anatomical image type
%    -> modality (mod-MTw/PDw/T1w), acquisition (acq-head/body)
% 2/ Deal with B1+ (B1 bias correction), 11 FA x 2 acquisition types
%    B1plus = emission with 22 images
%    -> echo index
% 3/ Deal with B0, 2 magnitude + 1 phase difference image
%    B0 = classic field maps. 2 magnitude + 1 phase diff
%    -> part-magnitude1/magnitude2/phasediff
% 
% The JSON files are BIDSified and updated for each type of fieldmap images
% because the 'IntendedFor' is fieldmap specific
%
% III. Create the general data set description
% Fill then create a 'dataset_description.json' file 
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

%__________________________________________________________________________
% Still missing/to check:
% =======================
% - turn this into a function with a few option flags as to how things are
%   done, especially what should be saved (full DICOM header?) or not.
% - add the required top metadata files describing the dataset and
%   participant
%   -> DONE
% - for the 'fmap' images at the "IntendedFor" metadata
%   -> DONE
% - gzip all images ('gzip' function)?
% - check filenames convention                \_ according to BEP001 crowd
% - check metadata list (names and definition)/
%   -> PARTLY DONE
% - check the units of metadata extracted
%   -> LOOKS LIKE IT's SORTED
% - the '_mod-*' is used for the B1- maps. These are not defined in BEP001
%   so far. Moreover the filename '_mod-*' field is currently ONLY defined
%   for the case of '_defacemask.nii' images!
%   -> COMBINE THINGS INTO THE '_acq-*' field -> DONE
% - the B1+ maps is defined but no B1- map in the BEP001 so far
% - for B1+ (and maybe B1-) maps, there is a specific suffix, B1plusmap
%   (and maybe B1minusmap). Or should it be 'B1plus' and B1minus'.
%__________________________________________________________________________
% NOTES on hMRI toolbox:
% ======================
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
% not usable since the "intermediate" files listed do not exist anymore at
% the end of the map creation process...
%__________________________________________________________________________

%% Some defaults

% Main folder with the data
rootPth = 'C:\Dox\2_Data\hmri_sample_dataset_with_maps';
% Subject label
subj_label = 'anon';
% Where the BIDS-ified data will be copied, bidsPth & subjPth
bidsPth = fullfile(rootPth,'BIDS_dataset_v3');
if ~exist(bidsPth,'dir'), mkdir(bidsPth), end;
subjPth = fullfile(bidsPth,sprintf('sub-%s',subj_label));
if ~exist(subjPth,'dir'), mkdir(subjPth), end;

%% I. Deal with the anatomical images

% 1/ Get the filenames for the 3 types (MTw, PDw, T1w) of anatomical images
% -------------------------------------------------------------------------
seq_label = { 'MTw' , 'PDw' , 'T1w' };
pth_MPM = { ...
    'mtw_mfc_3dflash_v1i_R4_0012' , ...
    'pdw_mfc_3dflash_v1i_R4_0009' , ...
    't1w_mfc_3dflash_v1i_R4_0015' };

% Get all the FP-filenames from the 3 folders for the 3 sequences
fn_MPM = cell(3,1) ; % empty cell array for filenames
for ii=1:3
    pth_ii = fullfile(rootPth,pth_MPM{ii});
    fn_MPM{ii} = spm_select('FPList',pth_ii,'^.*\.nii$');
end

% 2/ Move and rename anatomical images
% ------------------------------------
% Filename structure
fn_bids_struct = 'sub-%s_echo-%d_acq-%s_MPM.nii';
% feed in "subject label", "echo index", and "acqusition sequence".
% According to BEP001, multi-echo anatomical images acquired with the MPM
% protocal should be of type '_MPM'. The contrast (MT/PD/T1-weighted) is
% only an acquisition parameter -> in the '_acq-***_' part of the filename

% 'anat' folder: define & create as needed
anatPth = fullfile(subjPth,'anat');
if ~exist(anatPth,'dir'), mkdir(anatPth), end;

% created BIDS MPMw filenames
fn_bidsMPMw = cell(3,1);

% MPM files, .nii and .json
for iseq = 1:3 % 3 types of sequences (MTw, PDw, T1w) *in that order*!
    nMPMw_ii = size(fn_MPM{iseq},1);
    fn_bidsMPMw{iseq} = '';
    if nMPMw_ii
        for ii=1:nMPMw_ii
            % a) deal with image files
            % create fullpath filename for i^th anatomical image
            fn_bidsMPMw_ii = fullfile(anatPth, ...  % path
                sprintf( fn_bids_struct, ...            % filename
                subj_label, ii, seq_label{iseq} ) ); % feeds
            fn_bidsMPMw{iseq} = char(fn_bidsMPMw{iseq}, fn_bidsMPMw_ii);
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
    fn_bidsMPMw{iseq}(1,:) = [];
end

% 3/ update corresponding json file
% ---------------------------------
% Now update all JSON files
fn_jsonAnat = spm_select('FPList',anatPth,'^.*\.json$');
% In one go, using a single function to pull out all the potential fields
% necessary for the MPM data
fn_jsonAnat_BIDS = hmri_BIDSify_json(fn_jsonAnat);

% 4/ Need some hot-fix for the 'Manufacturer' field.
% --------------------------------------------------
% 3 names are returned -> keep 1st one only
for ii=1:size(fn_jsonAnat_BIDS)
    Json_ii = spm_load(fn_jsonAnat_BIDS(ii,:));
    if iscell(Json_ii.Manufacturer)
        Json_ii.Manufacturer = Json_ii.Manufacturer{1};
    end
    spm_save(deblank(fn_jsonAnat_BIDS(ii,:)),Json_ii,struct('indent','  '))
end

%% II. Deal with the field maps
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

% JSON files: keep history & acqpar + add 'IntendedFor'
BIDSjson_options = struct(...
    'keep_acqpar' , true, ...
    'keep_history', true, ...
    'addfield'    , 'IntendedFor') ;

% 1.a/ Deal with 6 B1- maps(RF sensitivity), 3 mod-type x 2 acq-coils
% -----------------------------------------------------------------
% Filename structure
fn_bids_B1m_struct = 'sub-%s_acq-%s_B1minus.nii';
% feed in "subject label", "acquisition type" (head/body), "modality" (MTw/PDw/T1w)

% Filenames for the 3X2 types: (MTw, PDw, T1w) x (head/body)
acq_label = { ...
    'MTwHead', 'MTwBody' ; ...
    'PDwHead', 'PDwBody' ; ...
    'T1wHead', 'T1wBody' };

% and the corresponding folder in the data set structure
pth_B1m = { ...
    'mfc_smaps_v1a_Array_0010' , 'mfc_smaps_v1a_QBC_0011' ;...
    'mfc_smaps_v1a_Array_0007' , 'mfc_smaps_v1a_QBC_0008' ;...
    'mfc_smaps_v1a_Array_0013' , 'mfc_smaps_v1a_QBC_0014' };
% new BIDS filenames
fn_bidsB1m_iAS_nii = cell(3,2);
fn_bidsB1m_iAS_json = cell(3,2);

for imod = 1:3 % (MTw, PDw, T1w)
    for iacq = 1:2 % (head/body)
        pth_iSA = fullfile(rootPth,pth_B1m{imod,iacq});
        % a) Deal with images
        % Select the image file requested.
        fn_iSA_nii = spm_select('FPList',pth_iSA,'^.*\.nii$');
        fn_bidsB1m_iAS_nii{imod,iacq} = ...
            fullfile( fmapPth, ...  % build fullpath filename
            sprintf( fn_bids_B1m_struct, ...            % filename
            subj_label, acq_label{imod,iacq} ) ); % feeds
        % copy file with name change
        copyfile(fn_iSA_nii,fn_bidsB1m_iAS_nii{imod,iacq});
        % b) Deal with JSON
        fn_iSA_json = spm_file(fn_iSA_nii,'ext','json');
        fn_bidsB1m_iAS_json{imod,iacq} = ...
            spm_file(fn_bidsB1m_iAS_nii{imod,iacq},'ext','json');
        % copy file with name change
        copyfile(fn_iSA_json,fn_bidsB1m_iAS_json{imod,iacq});
    end
end

% 1.b/ BIDSify JSON & fix it
% --------------------------
% 'Manufacturer' field has 3 names are returned -> keep 1st one only
% add "IntendedFor" fieldname, for each modality
for imod = 1:3 % (MTw, PDw, T1w)
    for iacq = 1:2 % (head/body)
        % BIDSify JSON
        fn_jsonfmat_BIDS = hmri_BIDSify_json( ...
            fn_bidsB1m_iAS_json{imod,iacq},BIDSjson_options);
        % Keep 1st name from Manufacturer
        Json_ii = spm_load(fn_jsonfmat_BIDS);
        if iscell(Json_ii.Manufacturer)
            Json_ii.Manufacturer = Json_ii.Manufacturer{1};
        end
        % Add IntendedFor information -> MPMw images per acuqisition
        Json_ii.IntendedFor = spm_file(fn_bidsMPMw{imod},'path','anat/');
        if strcmp(filesep,'\')
            % Make sure the relative path uses '/' and not '\'
            Json_ii.IntendedFor = ...
                char(regexprep(cellstr(Json_ii.IntendedFor),'\\','/'));
        end
        % Save the updated version
        spm_save(fn_jsonfmat_BIDS, Json_ii, struct('indent','  '))
    end
end

% 2/ Deal with B1+ (B1 bias correction), 11 FA x 2 acquisition types
% ------------------------------------------------------------------
% Filename structure
fn_bids_B1p_struct = 'sub-%s_acq-%s_fa-%02d_B1plus.nii';
% feed in "subject label", "acquisition type" (SE/STE), "flip angle" (index)
acq_label = {'SE', 'STE'};

% Select the images for the SE/STE acquisition
pth_B1p = fullfile(rootPth,'mfc_seste_b1map_v1e_0004');
fn_B1plus{1} = spm_select('FPList',pth_B1p,'^.*-1\.nii$');
fn_B1plus{2} = spm_select('FPList',pth_B1p,'^.*-2\.nii$');

% BIDSified filenames
fn_bidsB1p_iFA_nii = cell(1,2);
fn_bidsB1p_iFA_json = cell(1,2);

for iacq = 1:2
    fn_bidsB1p_iFA_nii{iacq} = '';
    fn_bidsB1p_iFA_json{iacq} = '';
    for ifa = 1:11
        % a) Deal with images
        fn_iFA_nii = deblank(fn_B1plus{iacq}(ifa,:));
        fn_bidsB1p_nii = fullfile( fmapPth, ... % path
            sprintf( fn_bids_B1p_struct, ...      % filename
            subj_label, acq_label{iacq}, ifa ) ); % feeds
        fn_bidsB1p_iFA_nii{iacq} = char(fn_bidsB1p_iFA_nii{iacq}, ...
            fn_bidsB1p_nii);
        % copy file with name change
        copyfile(fn_iFA_nii,fn_bidsB1p_nii);
        % b) Deal with JSON
        fn_iFA_json = spm_file(fn_iFA_nii,'ext','json');
        fn_bidsB1p_json = spm_file(fn_bidsB1p_nii,'ext','json');
        fn_bidsB1p_iFA_json{iacq} = char(fn_bidsB1p_iFA_json{iacq}, ...
            fn_bidsB1p_json);
        % copy file with name change
        copyfile(fn_iFA_json,fn_bidsB1p_json);
    end
    fn_bidsB1p_iFA_nii{iacq}(1,:) = [];
    fn_bidsB1p_iFA_json{iacq}(1,:) = [];
end

% 2.b/ BIDSify JSON & fix it
% --------------------------
% 'Manufacturer' field has 3 names are returned -> keep 1st one only
% add "IntendedFor" fieldname, for each modality

% BIDSify all at once
for iacq=1:size(fn_bidsB1p_iFA_json,2)
    % separate by acquisition type (STE vs. SE) to prevent character
    % padding
    fn_jsonfmat_BIDS_all = hmri_BIDSify_json( ...
        char(fn_bidsB1p_iFA_json(:,iacq)), BIDSjson_options);
    for ii=1:size(fn_jsonfmat_BIDS_all,1)
        % Keep 1st name from Manufacturer
        Json_ii = spm_load(fn_jsonfmat_BIDS_all(ii,:));
        if iscell(Json_ii.Manufacturer)
            Json_ii.Manufacturer = Json_ii.Manufacturer{1};
        end
        % Add IntendedFor information -> MPMw images
        Json_ii.IntendedFor = spm_file(char(fn_bidsMPMw(:)),'path','anat/');
        if strcmp(filesep,'\')
            % Make sure the relative path uses '/' and not '\'
            Json_ii.IntendedFor = ...
                char(regexprep(cellstr(Json_ii.IntendedFor),'\\','/'));
        end
        % Save the updated version
        spm_save(fn_jsonfmat_BIDS_all(ii,:), Json_ii, struct('indent','  '))
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

% Empty BIDSified filenames
fn_bidsB0_nii = '';
fn_bidsB0_json = '';

for ii=1:2
    % a) Deal with images
    fn_bidsB0ii_nii = fullfile( fmapPth, ...    % path
        sprintf( fn_bids_B0mag_struct, ...      % filename
                subj_label, ii ) );         % feeds
    % copy file with name change
    copyfile(deblank(fn_B0mag_nii(ii,:)),fn_bidsB0ii_nii);
    fn_bidsB0_nii = char(fn_bidsB0_nii, fn_bidsB0ii_nii);
    % b) Deal with JSON
    fn_B0mag_json = spm_file(fn_B0mag_nii(ii,:),'ext','json');
    fn_bidsB0ii_json = spm_file(fn_bidsB0ii_nii,'ext','json');
    % copy file with name change
    copyfile(fn_B0mag_json,fn_bidsB0ii_json);
    fn_bidsB0_json = char(fn_bidsB0_json, fn_bidsB0ii_json);
end

% Then the phase difference
% Filename structure
fn_bids_B0pdiff_struct = 'sub-%s_phasediff.nii';
B0pdiff_pth = fullfile(rootPth,'gre_field_mapping_1acq_rl_0006');
fn_B0pdiff_nii = spm_select('FPList',B0pdiff_pth,'^.*\.nii$');
% a) Deal with images
fn_bidsB0pdiff_nii = fullfile( fmapPth, ...    % path
    sprintf( fn_bids_B0pdiff_struct, subj_label ) ); % filename
% copy file with name change
copyfile(fn_B0pdiff_nii,fn_bidsB0pdiff_nii);
fn_bidsB0_nii = char(fn_bidsB0_nii, fn_bidsB0pdiff_nii);
fn_bidsB0_nii(1,:) = [];

% b) Deal with JSON
fn_B0pdiff_json = spm_file(fn_B0pdiff_nii,'ext','json');
fn_bidsB0pdiff_json = spm_file(fn_bidsB0pdiff_nii,'ext','json');
% copy file with name change
copyfile(fn_B0pdiff_json,fn_bidsB0pdiff_json);
fn_bidsB0_json = char(fn_bidsB0_json);
fn_bidsB0_json(1,:) = [];

fn_bidsB0pdiff_json = char(fn_bidsB0pdiff_json);

% 3.b/ BIDSify JSON & fix it
% --------------------------
% 'Manufacturer' field has 3 names are returned -> keep 1st one only
% add "IntendedFor" fieldname, for each modality
fn_bidsB0_json = hmri_BIDSify_json(fn_bidsB0_json, BIDSjson_options);
fn_bidsB0pdiff_json = hmri_BIDSify_json(fn_bidsB0pdiff_json, ...
    BIDSjson_options);

% convert to cell array of strings in order to prevent empty space padding
fn_bidsB0_json = [cellstr(fn_bidsB0_json); cellstr(fn_bidsB0pdiff_json)];

for ii=1:size(fn_bidsB0_json,1)
    % Keep 1st name from Manufacturer
    Json_ii = spm_load(char(fn_bidsB0_json(ii,:)));
    if iscell(Json_ii.Manufacturer)
        Json_ii.Manufacturer = Json_ii.Manufacturer{1};
    end
    % Add IntendedFor information -> B1minus images.
    Json_ii.IntendedFor = spm_file(char(fn_bidsB1m_iAS_json(:)),'path','fmap/');
    if strcmp(filesep,'\')
        % Make sure the relative path uses '/' and not '\'
        Json_ii.IntendedFor = ...
            char(regexprep(cellstr(Json_ii.IntendedFor),'\\','/'));
    end
    % Save the updated version
    spm_save(char(fn_bidsB0_json(ii,:)), Json_ii, struct('indent','  '))
end

% % 4/ Now update all JSON files from fmap folder
% % ---------------------------------------------
% fn_jsonfmat = spm_select('FPList',fmapPth,'^.*\.json$');
% % In one go, using a single function to pull out all the postential fields
% % necessary for the MPM data
% fn_jsonfmat_BIDS = hmri_BIDSify_json(fn_jsonfmat);

%% III. Create the general data set description
fn_dataset = fullfile(bidsPth,'dataset_description.json');

% Some references and links
RAL = { ...
    'http://www.hmri.info/' ; ...
    'https://doi.org/10.1016/j.dib.2019.104132' ; ...
    'https://doi.org/10.1016/j.neuroimage.2019.01.029' };
% Dataset info
ds_info = struct(...
    'Name', 'Example hMRI dataset', ...
    'BIDSVersion', 'based on BEP001 dev', ...
    'Authors', 'Martina Callaghan', ...
    'ReferencesAndLinks', {RAL} );
spm_save(fn_dataset, ds_info, struct('indent','  '))

% % Participant info
% fn_participant_tsv = fullfile(subjPth,'participants.tsv');
% fn_participant_json = fullfile(subjPth,'participants.json');
