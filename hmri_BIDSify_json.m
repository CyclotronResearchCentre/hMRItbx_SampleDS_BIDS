function fn_BIDS_json = hmri_BIDSify_json(fn_json,opt)
% Turn existing DICOM-imported JSON files into BIDS-compliant JSON files.
% 
% The point is to extract from the full DICOM-header structure, the few
% parameters that are actually needed by the hMRI-toolbox: 
% 1/ load the 'JSONtabl_dcm2bids.tsv' table for the 'conversion' of the 
%    fieldnames and value scaling
% 2/ loop through the JSON files: load, update and rewrite them
%
% Note that the 'JSONtabl_dcm2bids.tsv' file can be generated from the
% function hmri_create_JSONtabl.m or be manually edited.
% 
% FORMAT
%   fn_BIDS_json = hmri_BIDSify_json(fn_json,opt)
% 
% INPUT
%   fn_json : char array of JSON files with the full Dicom header
%   opt : option flag to keep ['keep', default] or remove ['remove'] the
%         full Dicom header in the JSON structure.
% 
% OUTPUT
%   fn_json : char array of JSON files with the BIDS compliant format
% 
%__________________________________________________________________________
% Copyright (C) 2018 GIGA Institute

% Written by C. Phillips, 2018.
% Cyclotron Research Centre, University of Liege, Belgium

%% Deal with input and parameters
% Keep existing JSON file content
if nargin<2, opt = 'keep'; end
% Select some JSON files if not provided
if nargin<1
    fn_json = spm_select(Inf,'^.*\.json$','Select JSON file(s)');
end

% Default JSON file conversion table filename
fn_JSONtabl = 'JSONtabl_dcm2bids.tsv';

% The .tsv file should be in the main folder of hMRI-toolbox
P = mfilename('fullpath');
fn_JSONtabl = fullfile(spm_file(P,'path'),fn_JSONtabl);
if ~exist(fn_JSONtabl,'file')
    fn_JSONtabl = hmri_create_JSONtabl(fn_JSONtabl);
end
list_metadata_MPM = spm_load(fn_JSONtabl);
nMetadata = numel(list_metadata_MPM.FieldnamesOriginal);

% Load original JSON files
nJson = size(fn_json,1);
All_mdStruc = get_metadata(fn_json);
All_mdStruc_BIDS = cell(1,nJson);


%% Do the job!

% Loop over all the JSON files
for ijson=1:nJson
    All_mdStruc_BIDS{ijson} = struct;
    for ii=1:nMetadata
        val = get_metadata_val( All_mdStruc{ijson}, ...
                list_metadata_MPM.FieldnamesOriginal{ii});
        if ~isempty(val)
            % turn scaling from char into number...
            sc = eval(list_metadata_MPM.Scaling{ii});
            if ~isnan( sc )
                val = val * sc;
            end
            All_mdStruc_BIDS{ijson}.(list_metadata_MPM.FieldnamesBIDS{ii}) = val;
        end
    end
    % Stick in at the end: the full header and history 
    fnm = fieldnames(All_mdStruc{ijson});
    for ii = 1:numel(fnm)
        % Only remove if requested and keep 'history' whatever.
        if ~strcmpi(opt,'remove') || strcmpi(fnm{ii},'history')
            All_mdStruc_BIDS{ijson}.(fnm{ii}) = All_mdStruc{ijson}.(fnm{ii});
        end
    end

    % Save structure in original 
    spm_save(deblank(fn_json(ijson,:)),All_mdStruc_BIDS{ijson},struct('indent','  '))

end

fn_BIDS_json = fn_json;

end


