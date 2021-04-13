%% Script to BIDSify the images from the hMRI example dataset
% 
% Work with the 2 subsets of data:
% - start with the anatomical images, 
% - then the field and RF sensitivity maps.
%
% It is assumed that the data are organized as downloaded
% -> we know which specific folder contains what sort of data
%
% For sanity check parameters are available in the mpm_params.mat file
% previously created, so one can get the metadata directly from
% - the Dicom header (saved in the current .json file) or
% - via direct pulling from these mpm_params.mat as they're used in the
%   hMRI toolbox
% This works similarly for the RFsens (B1-) and B1map (B0 and B1+) field
% maps. 
%
% Operations
% ==========
% I. Deal with the anatomical images
% ----------------------------------
% 1/ Get the filenames for the 3 types (MTw/PDw/T1w) of anatomical images
% 2/ Move and rename anatomical images
% 3/ update corresponding json file
%
% II. Deal with the field maps
% ----------------------------
% Work with the 3 different types of field maps one at a time:
% 1/ Deal with 6 B1- maps (RF sensitivity), 3 types x 2 coils
%    B1minus = RF sensitivity = suffix '_RB1COR'
%    -> spcecify 6 acquisitions : _acq-[MTw/PDw/T1w]/[head/body]
% 2/ Deal with B1+ (B1 emission bias correction), 11 FA x 2 echo types
%    B1plus = Transmit = suffix '_TB1EPI' 
%    -> specify _echo-[1/2] and _FA-[1...11] indexes, 
%       where echo-1/2 matches SE/STE echoes
% 3/ Deal with B0, 2 magnitude + 1 phase difference image
%    B0 = classic field maps with 2 magnitude + 1 phase diff images
%    -> specify _part-magnitude1/magnitude2/phasediff
% 
% The JSON files are BIDSified and updated for each type of fieldmap images
% because the 'IntendedFor' is fieldmap specific
%
% III. Create the general data set description
% --------------------------------------------
% Fill then create a 'dataset_description.json' file 
%
% Dependencies:
% =============
% - a function called "hmri_BIDSify_json.m" to refactor the metadata from
%   the DICOM header. This itself relies on the get_metadata_val.m function
%   from the hMRI toolbox
% - a tab-separated-value file, JSONtabl_dcm2bids.tsv, with the list of
%   metadata fields required.
% 
% Note:
% =====
% There is currently no clear guidelines on how to deal with the
% B1-/B1+/B0 "field maps" in BIDS 
% -> we will abide by the general BIDS principles, as well as those
% introduced for these multi-contrast/-echo structural images
% 
% MAJOR ASSUMPTIONS:
% ==================
% Since this is ad hoc code for the hMRI example Dataset, I make these
% assumptions about their acquistion order and properties:
% - for each contrast (MTw/PDw/T1w), the various echoes are acquired in
%   *ascending* order 
% - the flip angle will be low for the MTw/PDw images and high for the T1w
%   images
% => no need to check EchoTime/FlipAngle JSON values to define 'echo'/'fa' 
%    indexable metadata value when creating the BIDS-format JSON files and 
%    BIDS-renaming the image files
% Ideally these assumptions should be lifted: the actual echo time and flip
% angle would be read from the JSON files and the indexable filename
% metadata ordered accordingly.
%__________________________________________________________________________
% Copyright (C) 2019 GIGA Institute

% Written by C. Phillips, 2019.
% Cyclotron Research Centre, University of Liege, Belgium

%__________________________________________________________________________
% Still missing/to check:
% =======================
% - turn this into a function with a few option flags as to how things are
%   done, especially what should be saved (full DICOM header?) or not.
% - for the B1+/B1- images, the 'IntendedFor' field is repeated for all the 
%   files. It could be useful to put these in a separate JSON file, that 
%   would be applicable to all the B1plus.json files.
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

% Add subfolder with code to read metadata:
% addpath('C:\Dox\1_Code\hMRI_sampledataset_code\hMRItbx_SampleDS_BIDS\metadata')

% Main folder with the data
pth_Root = 'D:\ddd_Codes\Test_hMRI_BIDS\DATA\Source\ORIG';
% Subject label
subj_label = '01';
% Where the BIDS-ified data will be copied, pth_BIDS & pth_Subj
pth_BIDS = fullfile(pth_Root,'BIDS_dataset');
if ~exist(pth_BIDS,'dir'), mkdir(pth_BIDS), end;
pth_Subj = fullfile(pth_BIDS,sprintf('sub-%s',subj_label));
if ~exist(pth_Subj,'dir'), mkdir(pth_Subj), end;

% JSON files: keep 'history' but remove 'acqpar'
BIDSjson_options = struct(...
    'keep_acqpar' , false, ...
    'keep_history', true, ...
    'addfield'    , '') ;

%% I. Deal with the anatomical images

% 1/ Get the filenames for the 3 types (MTw, PDw, T1w) of anatomical images
% -------------------------------------------------------------------------
seq_label = { 'MTw' , 'PDw' , 'T1w' };
pth_MPM = { ...
    'mtw_mfc_3dflash_v1i_R4_0012' , ...
    'pdw_mfc_3dflash_v1i_R4_0009' , ...
    't1w_mfc_3dflash_v1i_R4_0015' };

% Get all the FP-filenames from the 3 folders for the 3 sequences
fn_MPM = cell(3,1) ; % empty cell array for filenames (MTw, PDw, T1w)
for ii=1:3
    pth_ii = fullfile(pth_Root,pth_MPM{ii});
    fn_MPM{ii} = spm_select('FPList',pth_ii,'^.*\.nii$');
end

% 2/ Move and rename anatomical images
% ------------------------------------
% Filename structure
fn_bids_struct = 'sub-%s_acq-%s_echo-%d_flip-%d_mt-%s_MPM.nii';
% Feed in 
% - 'sub' -> subject label, 
% - 'acq' -> free-form contrast info (MTw,PDw,T1w), kept for human readability
% - 'echo' -> echo index,
% - 'flip' -> flip angle index,
% - 'mt' -> magnetization transfer 'on' or 'off'
% 
% Since this is an MPM acquisition protocol, 
% - the filename suffix is '_MPM'. 
% - key "multiple parameters" info are in the indexable key/value pairs
%   'echo', 'fa', and 'mt'. 
% - image contrast (MTw/PDw/T1w) is only optional here and set in the 'acq'
%   key/value pair.
% 
% Assumptions:
% - for each contrast (MTw/PDw/T1w), the various echoes are acquired in
%   *ascending* order 
%   -> use their listing order for the 'echo' indexable metadata
% - the flip angle will be low for the MTw/PDw images and high for the T1w
%   images 
%   -> define a priori 'fa-1' and 'fa-2' for MTw/PDw and T1w resp.
% - magentization transfer only occurs for the MTw images
%   -> define a priori 'mt-on' and 'mt-off' for MTw and PDw/T1w resp.
fa_mt_values = { ...    % for
    1 , 'on'  ; ...     %   MTw
    1 , 'off' ; ...     %   PDw
    2 , 'off'};         %   T1w

% 'anat' folder: define & create as needed
pth_anat = fullfile(pth_Subj,'anat');
if ~exist(pth_anat,'dir'), mkdir(pth_anat), end;

% created BIDS MPMw filenames
fn_bidsMPMw = cell(3,1);

% MPM files, .nii and .json
for iseq = 1:3 % 3 types of sequences (MTw, PDw, T1w) *in that order*!
    nMPMw_ii = size(fn_MPM{iseq},1);
    fn_bidsMPMw{iseq} = '';
    if nMPMw_ii
        for ii=1:nMPMw_ii
            % a) deal with image files
            % create filename for i^th anatomical image
            fn_bidsMPMw_ii = sprintf( fn_bids_struct, ...  
                subj_label, seq_label{iseq}, ii, fa_mt_values{iseq,1}, ...
                fa_mt_values{iseq,2} ) ; % feeds
            % add the path
            fn_bidsMPMw_ii = fullfile(pth_anat, fn_bidsMPMw_ii ); 
            % stack up
            fn_bidsMPMw{iseq} = char(fn_bidsMPMw{iseq}, fn_bidsMPMw_ii);
            % copy file with name change
            copyfile(fn_MPM{iseq}(ii,:),fn_bidsMPMw_ii);
            % b) deal with JSON files
            % -> copy with right name then convert them all later on!
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
fn_jsonAnat = spm_select('FPList',pth_anat,'^.*\.json$');
% In one go, using a single function to pull out all the potential fields
% necessary for the MPM data
fn_jsonAnat_BIDS = hmri_BIDSify_json(fn_jsonAnat, BIDSjson_options);

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
% 1/ Deal with 6 B1- maps (RF sensitivity), 3 types x 2 coils
%    B1minus = RF sensitivity = suffix '_RB1COR'
%    -> spcecify 6 acquisitions : _acq-[MTw/PDw/T1w]/[head/body]
% 2/ Deal with B1+ (B1 emission bias correction), 11 FA x 2 echo types
%    B1plus = Transmit = suffix '_TB1EPI' 
%    -> specify _echo-[1/2] and _FA-[1...11] indexes, 
%       where echo-1/2 matches SE/STE echoes
% 3/ Deal with B0, 2 magnitude + 1 phase difference image
%    B0 = classic field maps with 2 magnitude + 1 phase diff images
%    -> specify _part-magnitude1/magnitude2/phasediff

% 0/ Create appropriate 'fmap' folder
% 'fmap' folder: define & create as needed
pth_fmap = fullfile(pth_Subj,'fmap');
if ~exist(pth_fmap,'dir'), mkdir(pth_fmap), end;

% JSON files: add 'IntendedFor' field
BIDSjson_options.addfield = 'IntendedFor' ;

% 1/ Deal with B1- maps(RF sensitivity), 
% --------------------------------------
% Keeping in mind there are 6 images: 3 types x 2 coils

% 1.a/ Deal with all the file moving & renaming
% .............................................
% Filename structure
fn_bids_B1m_struct = 'sub-%s_acq-%s_RB1COR.nii';
% feed in 
% - "subject label", 
% - "acquisition type" (head/body) \_ in the '_acq-' field
% - "modality" (MTw/PDw/T1w)       /
% One option could be to use the '_mod-' field but it is reserved for the
% "defacemask" suffix apparently. See issue #3:
% https://github.com/CyclotronResearchCentre/hMRItbx_SampleDS_BIDS/issues/3

% Filename tags for the 3X2 types: (MTw, PDw, T1w) x (head/body)
acq_label = { ...
    'headMTw', 'bodyMTw' ; ...
    'headPDw', 'bodyPDw' ; ...
    'headT1w', 'bodyT1w' };

% and the corresponding folders in the data set structure
pth_B1m = { ...
    'mfc_smaps_v1a_Array_0010' , 'mfc_smaps_v1a_QBC_0011' ;...
    'mfc_smaps_v1a_Array_0007' , 'mfc_smaps_v1a_QBC_0008' ;...
    'mfc_smaps_v1a_Array_0013' , 'mfc_smaps_v1a_QBC_0014' };
% new BIDS filenames
fn_bidsB1m_iAS_nii = cell(3,2);
fn_bidsB1m_iAS_json = cell(3,2);

for imod = 1:3 % (MTw, PDw, T1w)
    for iacq = 1:2 % (head/body)
        pth_iSA = fullfile(pth_Root,pth_B1m{imod,iacq});
        % a) Deal with images
        % Select the image file requested.
        fn_iSA_nii = spm_select('FPList',pth_iSA,'^.*\.nii$');
        fn_bidsB1m_iAS_nii{imod,iacq} = ...
            fullfile( pth_fmap, ...  % build fullpath filename
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
% ..........................
% Fixes:
% - 'Manufacturer' field has 3 names are returned -> keep 1st one only
% - add "IntendedFor" fieldname, for each modality
for imod = 1:3 % (MTw, PDw, T1w)
    for iacq = 1:2 % (head/body)
        % BIDSify JSON
        fn_jsonfmat_BIDS = hmri_BIDSify_json( ...
            fn_bidsB1m_iAS_json{imod,iacq}, BIDSjson_options);
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

% 2/ Deal with B1+ (B1 bias correction)
% -------------------------------------
% Keeping in mind there are 22 images: 11 FA x 2 echo types
% where echo-1/2 matches SE/STE echoes

% 2.a/ Deal with all the file moving & renaming
% .............................................
% Filename structure
fn_bids_B1p_struct = 'sub-%s_echo-%1d_flip-%02d_TB1EPI.nii';
% feed in 
% - "subject label", 
% - "echo (index) for SE/STE echoes (1/2)
% - "flip angle" (index), in increasing order of values for the FA!

% Select the images for the SE/STE acquisition
pth_B1p = fullfile(pth_Root,'mfc_seste_b1map_v1e_0004');
fn_B1p{1} = spm_select('FPList',pth_B1p,'^.*-1\.nii$');
fn_B1p{2} = spm_select('FPList',pth_B1p,'^.*-2\.nii$');

% Get the "B1map Nominal FA Values" (nomFA) then 
% 1) order the files in ascending order value
% 2) push these values in the FA field in JSON afterwards
% Note that all the values are saved in each original JSON file

nomFA = get_metadata_val(fn_B1p{1}(1,:), 'B1mapNominalFAValues');
[s_nomFA, ind_nomFA] = sort(nomFA); 
% use i_nomFA to properly index FA in BIDS-name

% BIDSified filenames
fn_bidsB1p_iFA_nii = cell(1,2);
fn_bidsB1p_iFA_json = cell(1,2);

for iacq = 1:2
    fn_bidsB1p_iFA_nii{iacq} = '';
    fn_bidsB1p_iFA_json{iacq} = '';
    for ifa = 1:11
        % a) Deal with images
        fn_iFA_nii = deblank(fn_B1p{iacq}(ifa,:));
        fn_bidsB1p_nii = fullfile( pth_fmap, ... % path
            sprintf( fn_bids_B1p_struct, subj_label, iacq, ind_nomFA(ifa) ) ); % fname
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
% ..........................
% Fixes:
% - 'Manufacturer' field has 3 names are returned -> keep 1st one only
% - add "IntendedFor" fieldname, for each modality
% - remove 'NominalFAValues' and put the corresponding nomFA in 'FlipAngle'
% - add 'FlipAngleSeries' with sorted 'NominalFAValues' (this is not 100% 
%   BIDS but useful for the hMRI toolbox)

% BIDSify all at once
fn_jsonfmat_BIDS_all = hmri_BIDSify_json( ...
    char(fn_bidsB1p_iFA_json(:)), BIDSjson_options);

% Fix JSON files
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
    
    % Remove 'NominalFAValues' and put correspondin nomFA in 'FlipAngle'
    ind_fn = rem(ii-1,11)+1; % figure out index [1 11]
    if isfield(Json_ii,'NominalFAValues')
        Json_ii = rmfield(Json_ii,'NominalFAValues');
    end
    Json_ii.FlipAngle = s_nomFA(ind_nomFA(ind_fn)); % as fa index are sorted
    % equivalent to Json_ii.FlipAngle = nomFA(ind_fn)
    
    % Add 'FlipAngleSeries' with sorted 'NominalFAValues' 
    Json_ii.FlipAngleSeries = s_nomFA;
    
    % Save the updated version
    spm_save(fn_jsonfmat_BIDS_all(ii,:), Json_ii, struct('indent','  '))
end

% 3/ Deal with B0
% ---------------
% Keeping in mind, there are only 3 images: 2 mag + 1 phase diff images
% Proceed 
% - first the 2 magnitude images
% - then the phase difference

% 3.a/ Deal with the 2 "magnitude" images moving & renaming
% .........................................................
% Filename structure
fn_bids_B0mag_struct = 'sub-%s_magnitude%d.nii';
% feed in "subject label", "magnitude index"
pth_B0mag = fullfile(pth_Root,'gre_field_mapping_1acq_rl_0005');
fn_B0mag_nii = spm_select('FPList',pth_B0mag,'^.*\.nii$');

% Empty BIDSified filenames
fn_bidsB0_nii = '';
fn_bidsB0_json = '';

for ii=1:2
    % i) Deal with images
    fn_bidsB0ii_nii = fullfile( pth_fmap, ...    % path
        sprintf( fn_bids_B0mag_struct, ...      % filename
                subj_label, ii ) );         % feeds
    % copy file with name change
    copyfile(deblank(fn_B0mag_nii(ii,:)),fn_bidsB0ii_nii);
    fn_bidsB0_nii = char(fn_bidsB0_nii, fn_bidsB0ii_nii);
    % ii) Deal with JSON
    fn_B0mag_json = spm_file(fn_B0mag_nii(ii,:),'ext','json');
    fn_bidsB0ii_json = spm_file(fn_bidsB0ii_nii,'ext','json');
    % copy file with name change
    copyfile(fn_B0mag_json,fn_bidsB0ii_json);
    fn_bidsB0_json = char(fn_bidsB0_json, fn_bidsB0ii_json);
end

% 3.b/ Deal with the "phasediff" image moving & renaming
% ......................................................
% Filename structure
fn_bids_B0pdiff_struct = 'sub-%s_phasediff.nii';
pth_B0pdiff = fullfile(pth_Root,'gre_field_mapping_1acq_rl_0006');
fn_B0pdiff_nii = spm_select('FPList',pth_B0pdiff,'^.*\.nii$');
% i) Deal with images
fn_bidsB0pdiff_nii = fullfile( pth_fmap, ...    % path
    sprintf( fn_bids_B0pdiff_struct, subj_label ) ); % filename
% copy file with name change
copyfile(fn_B0pdiff_nii,fn_bidsB0pdiff_nii);
fn_bidsB0_nii = char(fn_bidsB0_nii, fn_bidsB0pdiff_nii);
fn_bidsB0_nii(1,:) = [];

% ii) Deal with JSON
fn_B0pdiff_json = spm_file(fn_B0pdiff_nii,'ext','json');
fn_bidsB0pdiff_json = spm_file(fn_bidsB0pdiff_nii,'ext','json');
% copy file with name change
copyfile(fn_B0pdiff_json,fn_bidsB0pdiff_json);
fn_bidsB0_json = char(fn_bidsB0_json, fn_bidsB0pdiff_json);
fn_bidsB0_json(1,:) = [];

% 3.c/ BIDSify JSON & fix it
% ..........................
% Fixes:
% - 'Manufacturer' field has 3 names are returned -> keep 1st one only
% - add "IntendedFor" fieldname, for each modality
fn_bidsB0_json = hmri_BIDSify_json(fn_bidsB0_json, BIDSjson_options);

for ii=1:size(fn_bidsB0_json,1)
    % Keep 1st name from Manufacturer
    Json_ii = spm_load(fn_bidsB0_json(ii,:));
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
    spm_save(fn_bidsB0_json(ii,:), Json_ii, struct('indent','  '))
end

%% III. Create the general data set description
fn_dataset = fullfile(pth_BIDS,'dataset_description.json');

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
% fn_participant_tsv = fullfile(pth_Subj,'participants.tsv');
% fn_participant_json = fullfile(pth_Subj,'participants.json');
