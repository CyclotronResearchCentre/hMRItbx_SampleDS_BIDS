function fn_JSONtabl = hmri_create_JSONtabl(fn_JSONtabl, opt)
% Function to generate the list of BIDS metadata fields to retain for the
% MPM protocol in the JSON header.
% Since the conversion itself will rely on the 'get_metadata_val.m', it is
% important to define the scaling according to the units returned from
% 'get_metadata_val.m'.
%
% It works like this:
% 1/ list the metadata needed,
%   * by their 'search' name in the original Dicom-hdr JSON file
%   * with corresponding BIDS JSON file
%   in a single cell array,
% 2/ save it as a .tsv file for further use
% The table includes 2 main columns:
% - 'FieldnamesOriginal', fieldname in the Dicom header
% - 'FieldnamesBIDS', fieldname for BIDS
%
% NOTE:
% - Previous version included some scaling, which is dropped now.
% - there remain some issues wrt. some BIDS metadata, see further down.
%
% FORMAT
%   fn_JSONtabl = hmri_create_JSONtabl(fn_JSONtabl, opt)
%
% INPUT
%   fn_JSONtabl : filename of .tsv file to create.
%                 By default, <path_to_hmri>\JSONtabl_dcm2bids.tsv
%   opt : option field to create a full list of metadata ['full'] or just
%         those for the MPM protocal ['MPM', default]
%
%__________________________________________________________________________
% Copyright (C) 2018 GIGA Institute

% Written by C. Phillips, 2018.
% Cyclotron Research Centre, University of Liege, Belgium

%% Deal with input

% No filename passed -> create a generic one
fn_JSONtabl_defaults = 'JSONtabl_dcm2bids.tsv';
if nargin<1 || ~ischar(fn_JSONtabl) || isempty(fn_JSONtabl)
    P = mfilename('fullpath');
    fn_JSONtabl = fullfile(spm_file(P,'path'),fn_JSONtabl_defaults);
end

% Deal with option: 'full' or 'MPM' [default]
if nargin <2
    opt = 'MPM';
end

%% Create table
% =============
metadata_origin = {'FieldnamesOriginal', 'FieldnamesBIDS'} ; 

% 1/ MPM parameters only, 
% -----------------------
% from the calls to 'get_metadata_val.m' in the 'hmri_create_*.m' functions
% and which are used to create the qMRI.
metadata_MPM = {...
    ... % Some names for scanning, sequence, etc.
    'ScanningSequence',             'ScanningSequence'          ; ... % e.g. 'EP' for EPI
    'SequenceName',                 'SequenceName'              ; ...
    'ProtocolName',                 'PulseSequenceDetails'      ; ...
    ... % Main parameters: timing and RF
    'RepetitionTime',               'RepetitionTimeExcitation'  ; ... % TR [s]
    'EchoTime',                     'EchoTime'                  ; ... % TE [s]
    'FlipAngle',                    'FlipAngle'                 ; ... % FlipAngle [deg]
    ... % Sequence specific
    'MT',                           'MTState'                   ; ... % MT pulse [0/1]
    ... % In-Plane Spatial Encoding
    'NumberOfMeasurements',         'NumberShots'               ; ... % #RF excitations to reconstruct slove/volume
    'epiReadoutDuration',           'TotalReadoutTime'          ; ... % [s]
    'PhaseEncodingDirectionSign',   'PhaseEncodingDirectionSign'; ... % A>>P & R>>L = 1; P>>A & L>>R = 0
    ... % add for 3D_EPI used for al_B1mapping
    'B1mapNominalFAValues',         'NominalFAValues'           ; ... % [deg]
    'B1mapMixingTime',              'MixingTime'                ; ... % [s]
    ... % Further stuff for precalculation of RF correction, defined in al_B1mapping and mtflash3d sequences
    'RFSpoilingPhaseIncrement',     'SpoilingRFPhaseIncrement'  ; ... % [deg] 
    'spoilingGradientMoment',       'SpoilingGradientMoment'    ; ... % [T*s/m]
    'spoilingGradientDuration',     'SpoilingGradientDuration'    ... % [s]
    };

% BEP001 NOTE:
% ~~~~~~~~~~~~ 
% B1map parameters are not all defined in BEP001 yet or definition is not 
% clear enough:
% - new metadata fieldnames
%   * B1mapNominalFAValues -> NominalFAValues \_ B1 mapping specific in MPM
%   * B1mapMixingTime -> MixingTime           /
%   * PhaseEncodingDirectionSign -> new one. 
%     There exist a 'PhaseEncodingDirection' fieldname that encodes both 
%     direction and sign, with value i,j,k,i-,j-,k-, but it seems rather
%     complicated to extract....
% - should be "better defined"(?) fieldnames
%   * epiReadoutDuration -> TotalReadoutTime, assuming this fits the BEP001
%     definition but it seems a bit vague as total RO time might be defined
%     differently for different sequence. Plus what about parallel 
%     acceleration techniques, accounted for or not?

% 2/ Add some common stuff:
% -------------------------
extra_md = { ...
    'Manufacturer',          'Manufacturer' ; ...
    'ManufacturerModelName', 'ManufacturersModelName' ; ... % note the 's' difference !
    'DeviceSerialNumber',    'DeviceSerialNumber' ; ...
    'StationName',           'StationName' ; ...
    'MagneticFieldStrength', 'MagneticFieldStrength' ...
    };
% since they are general parameters, put them on top
metadata_MPM = [extra_md ; metadata_MPM ];

% 3/ Add a few more fieldnames, if required
% ----------------------------
% Those are explicitly listed in the 'get_metadata_val.m' function but not
% directly used in the current version of hMRI
if strcmpi(opt,'full')
    full_md = { ...
        'RepetitionTimes',      'RepetitionTimesExcitation' ; ... % TRs [s], note the extra 's' !
        'BandwidthPerPixelRO',  'BandwidthPerPixelRO'       ; ... % [Hz/pxl]
        'PELinesPAT',           'PELinesPAT'                  ...
        ... % size of the k-space PE dimension, taking into account Parallel
        ... % acceleration but not partial Fourier. Used to calculate the total
        ... % EPI Readout duration for FieldMap undistortion.
        };
    metadata_MPM = [metadata_MPM ; full_md];
end
% BEP001 NOTE:
% ~~~~~~~~~~~~ 
% These fields are not necessarily defined within BEP001.

% 4/ Add a few more recommended fields from BEP001
% extra_md = {};
% metadata_MPM = [metadata_MPM ; extra_md];

% 5/ Save in .tsv file to be used later on
% --------------------
% Put the top line with origin
metadata_MPM = [metadata_origin ; metadata_MPM];
spm_save(fn_JSONtabl, metadata_MPM)

end