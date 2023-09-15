# hMRI toolbox Sample Data Set BIDS
Bits of code to BIDS-ify the example MPM data set, [reference paper](https://doi.org/10.1016/j.dib.2019.104132) and [original data](https://owncloud.gwdg.de/index.php/s/iv2TOQwGy4FGDDZ), from the [hMRI toolbox](https://hmri.info/).

## Original data description

These are images acquired from a single subject, following an MPM acquisition protocol, see the [parameters here](https://owncloud.gwdg.de/index.php/s/iv2TOQwGy4FGDDZ/download?path=%2F&files=hmri_sample_dataset_protocol_800um_64ch.pdf) from the server. These are also detailed and further explained in [this paper](https://doi.org/10.1016/j.dib.2019.104132). The data are provided in **NIfTI+JSON format** but importantly:

- the structural images have been defaced and the data are anonymized

- the original DICOM files are not available

- the different file types are saved in specific sub-folders but their naming is not obvious, except for the MT/PD/T1 weighted structural images (first 3 chars of last 3 folder names)

  ```
  gre_field_mapping_1acq_rl_0005
  gre_field_mapping_1acq_rl_0006
  mfc_seste_b1map_v1e_0004      
  mfc_smaps_v1a_Array_0007      
  mfc_smaps_v1a_Array_0010      
  mfc_smaps_v1a_Array_0013      
  mfc_smaps_v1a_QBC_0008        
  mfc_smaps_v1a_QBC_0011        
  mfc_smaps_v1a_QBC_0014        
  mtw_mfc_3dflash_v1i_R4_0012   
  pdw_mfc_3dflash_v1i_R4_0009               
  t1w_mfc_3dflash_v1i_R4_0015   
  ```

- the filenames are those from the console/DICOM files, with their subject and acquisition index, i.e. not easy to interpret. For example, the 6 MT-weighted images with JSON files:

  ```
  mtw_mfc_3dflash_v1i_R4_0012\
  	anon_s2018-02-28_18-26-190132-00001-00224-1.json
  	anon_s2018-02-28_18-26-190132-00001-00224-1.nii 
  	anon_s2018-02-28_18-26-190132-00001-00448-2.json
  	anon_s2018-02-28_18-26-190132-00001-00448-2.nii 
  	anon_s2018-02-28_18-26-190132-00001-00672-3.json
  	anon_s2018-02-28_18-26-190132-00001-00672-3.nii 
  	anon_s2018-02-28_18-26-190132-00001-00896-4.json
  	anon_s2018-02-28_18-26-190132-00001-00896-4.nii 
  	anon_s2018-02-28_18-26-190132-00001-01120-5.json
  	anon_s2018-02-28_18-26-190132-00001-01120-5.nii 
  	anon_s2018-02-28_18-26-190132-00001-01344-6.json
  	anon_s2018-02-28_18-26-190132-00001-01344-6.nii 
  ```
  
- the JSON files contain contain 2 main fields `history` and `acqpar` respectively providing 

  - information on the DICOM-to-NIfTI conversion step with SPM12 `spm_dicom_convert.m` function
  - a copy of the whole DICOM header, i.e. with all the meta data, included system/sequence specific ones

The organisation of the data in 12 folders is overall complex to handle, and there is little information in the folder-/filenames to know what they contain...



## Objective

The goal is thus to make this "raw" set of data BIDS compliant, according to the BIDS extension for qunatitative MRI, see the [qMRI BIDS specs](https://bids-specification.readthedocs.io/en/stable/appendices/qmri.html) and the [Karakuzu et al. (2022) reference paper](https://doi.org/10.1038/s41597-022-01571-4). This involves:

- reorganizing the files into BIDS-defined folders, here `anat` and `fmap`
- renaming the image+JSON files using BIDS nomenclature for multi-echo, multi-FA, image type, etc.
- rewrite the JSON files such that they contain explicitly defined meta-data keywords

The list of meta-data field names used in the hMRI toolbox/MPM protocol versus what's available in BIDS specs are compared in this [`MPM_parameter_list.md`](https://github.com/CyclotronResearchCentre/hMRItbx_SampleDS_BIDS/blob/master/MPM_parameter_list.md) file. Some are missing from the BIDS specs, especially relating to B1 maps.

The resulting BIDS-ified dataset is available in the OSF [qMRI-BIDS example dataset](https://osf.io/k4bs5/), see the `ds-mpm` sub-folder.

## Scripts and bits of code

The repository contains 2 scripts and 1 function:
- [`script_BIDS_MPM_exampleDataset.m`](https://github.com/CyclotronResearchCentre/hMRItbx_SampleDS_BIDS/blob/master/script_BIDS_MPM_exampleDataset.m), to BIDS-ify the example dataset
- [`hmri_BIDSify_json`](https://github.com/CyclotronResearchCentre/hMRItbx_SampleDS_BIDS/blob/master/hmri_BIDSify_json.m), to BIDS-ify the JSON files (which originally contain the DICOM header and history)
- [`hmri_create_JSONtabl`](https://github.com/CyclotronResearchCentre/hMRItbx_SampleDS_BIDS/blob/master/hmri_create_JSONtabl.m), to build the list of meta-data

These are further described here

### 1/ `script_BIDS_MPM_exampleDataset.m` to BIDS-ify the data

This [script](https://github.com/CyclotronResearchCentre/hMRItbx_SampleDS_BIDS/blob/master/script_BIDS_MPM_exampleDataset.m) BIDS-ify the example dataset, as long as it knows where the main folder is located.

There nevertheless remain a [few issues](https://github.com/CyclotronResearchCentre/hMRItbx_SampleDS_BIDS/issues) to sort as well as room for improvements, for example:

- check the list of meta-data and their definition. Some need to be added in the BEP001 proposal

- completing some JSON files, especially the fieldmap ones to link these with their corresponding anatomical images, using the `IntendedFor` field.

- add an option to zip image files as they're BIDS-ified. 

Note too that the script relies on :
- the `hmri_BIDSify_json.m` function to refactor the meta-data from  the DICOM header, especially the acquisition parameters. This itself relies on the `get_metadata_val.m` function from the hMRI toolbox
- a tab-separated-value file, `JSONtabl_dcm2bids.tsv`, with the list of meta-data fields required.

These could certainly be improved too.

### 2/ `hmri_BIDSify_json.m` to handle acquisition parameter metadata

The point is to extract the from the original JSON files, containing a "dump" of the DICOM header and the history of the file, the meta-data necessary to reconstruct the quantitative maps, i.e. mostly acquisition parameters, and store those according to BIDS specs.
The list of parameters used/needed by the hMRI toolbox is actually defined in the [`JSONtabl_dcm2bids.tsv`](https://github.com/CyclotronResearchCentre/hMRItbx_SampleDS_BIDS/blob/master/JSONtabl_dcm2bids.tsv) table, which itself is generated by the [`hmri_create_JSONtabl`](https://github.com/CyclotronResearchCentre/hMRItbx_SampleDS_BIDS/blob/master/hmri_create_JSONtabl.m) function.

### 3/ `hmri_create_JSONtabl.m` to list the acquisition meta-data

This is a simple and modular way to create a table of (acquisition) meta-data, necessary for the quantitative map calculation. Since meta-data fetching was already implemented in the hMRI toolbox (through the `get_metadata_val.m` function), this links the BIDS meta-data field names with those from the (original) hMRI toolbox.

Note this list is not definitive, so the function could be updated as seemed fit.

---

## Creating maps

Scripts to efficiently create the batch for map creation from the downloaded example dataset, in original or BIDS format, is available in a separate Github repo named [`hMRItoolbox_ExampleDataset_script`](https://github.com/CyclotronResearchCentre/hMRItoolbox_ExampleDataset_script)


Being able to quickly processing the original data is still useful, e.g. to set up a ground truth of the results. Unfortunately in the batch file provided with the example dataset **the path to the various images are hard-coded**  and therefore will not run as is.

