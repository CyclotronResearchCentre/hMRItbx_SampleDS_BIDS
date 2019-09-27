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
%   * plus some scaling factor, e.g. [ms] -> [s], if necessary
%   in a single cell array,
% 2/ save it as a .tsv file for further use
% The table includes 3 main columns:
% - 'FieldnamesOriginal', fieldname in the Dicom header
% - 'FieldnamesBIDS', fieldname for BIDS
% - 'Scaling', scaling can be useful to match BIDS units (i.e. SI units),
%   for example [ms]->[s] means a scaling by .001.
%   Use 'NaN' if no scaling is required, e.g. for a char/string
%
% NOTE:
% the latest version of the DCM-to-NIfTI conversion works with SI units, so
% scaling values are set to 1.
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

% NOTE:
% most scaling has been included in the get_metadata_val.m function, so all
% values are set to 1 here. Still I keep this feature in as it could prove
% useful at some point in the future...

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

% MPM parameters only, from the calls to 'get_metadata_val.m' in the
% 'hmri_create_*.m' functions.
metadata_MPM = {...
    'FieldnamesOriginal',           'FieldnamesBIDS',           'Scaling' ; ...
    ... % Main parameters
    'RepetitionTime',               'RepetitionTimeExcitation', 1 ; ... % TR [s] -> [s]
    'EchoTime',                     'EchoTime',                 1 ; ... % TE [s] -> [s]
    'FlipAngle',                    'FlipAngle',                1 ; ...    % [deg]
    ... % add for 3D_EPI
    'B1mapNominalFAValues',         'B1mapNominalFAValues',     1 ; ...    % for al_B1mapping [deg]
    'B1mapMixingTime',              'MixingTime',               1 ; ... % for al_B1mapping [s] -> [s]
    'epiReadoutDuration',           'epiReadoutDuration',       1 ; ... % [s] -> [s]
    'PhaseEncodingDirectionSign',   'PhaseEncodingDirectionSign',NaN ; ...  % A>>P & R>>L = 1; P>>A & L>>R = 0
    ... % Some names for scanning, sequence, etc.
    'ScanningSequence',             'ScanningSequence',         NaN ; ...  % e.g. 'EP' for EPI
    'SequenceName',                 'SequenceName',             NaN ; ...
    'ProtocolName',                 'ProtocolName',             NaN ; ...
    ... % Further stuff for precalculation of RF correction
    'RFSpoilingPhaseIncrement',     'RFSpoilingPhaseIncrement', 1 ; ...   % [°] defined in al_B1mapping and mtflash3d sequences
    'spoilingGradientMoment',       'spoilingGradientMoment',   1 ; ...  % [T*s/m] defined in al_B1mapping and mtflash3d sequences
    'spoilingGradientDuration',     'spoilingGradientDuration', 1 ...   % [s] defined in al_B1mapping and mtflash3d sequences
    };

% If required, add a few more fieldnames
% Those are explicitly listed in the 'get_metadata_val.m' function
if strcmpi(opt,'full')
    full_md = { ...
        'RepetitionTimes', 'RepetitionTimesExcitation', 1 ; ... % TRs [s] -> [s], note the extra 's' !
        'BandwidthPerPixelRO',          'BandwidthPerPixelRO',          1 ; ...    % [Hz/pxl]
        'PELinesPAT',                   'PELinesPAT',                   1 ; ...
        ... % size of the k-space PE dimension, taking into account Parallel
        ... % acceleration but not partial Fourier. Used to calculate the total
        ... % EPI Readout duration for FieldMap undistortion.
        'NumberOfMeasurements',         'NumberOfMeasurements',         NaN ...  % Siemens-specific
        };
    metadata_MPM = [metadata_MPM ; full_md];
end

% Add some common stuff:
extra_md = { ...
    'Manufacturer', 'Manufacturer', NaN ; ...
    'ManufacturerModelName', 'ManufacturersModelName', NaN ; ... % note the 's' difference !
    'DeviceSerialNumber', 'DeviceSerialNumber', NaN ; ...
    'MagneticFieldStrength', 'MagneticFieldStrength', NaN };
metadata_MPM = [metadata_MPM ; extra_md];

% Save in .tsv file to be used later on
spm_save(fn_JSONtabl, metadata_MPM)

end