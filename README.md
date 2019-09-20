# hMRI toolbox Sample Data Set BIDS
Bits of code to BIDS-ify the [example data set](https://owncloud.gwdg.de/index.php/s/iv2TOQwGy4FGDDZ) from the [hMRI toolbox](https://hmri-group.github.io/hMRI-toolbox/).

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

## Objective

The goal is thus to make this "raw" set of data BIDS compliant, according to the [BEP001 proposal](https://github.com/bids-standard/bep001) that is currently being developed.  This involves:

- reorganizing the files into BIDS-defined folders, here `anat` and `fmap`
- renaming the image+JSON files using BIDS nomenclature for multi-echo, multi-FA, image type, etc.
- rewrite the JSON files such that they contain explicitly defined meta-data keywords

