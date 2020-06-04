# `metadata` folder

The content of this folder is actually a subset of functions from the full `hMRI` toolbox code, from a `BIDSify` branch. These functions are effectively available (or will be in a next release) located in this folder:

````
C:\Mycode\hMRI_toobox\spm12\metadata
````

## List of functions are their role

- `find_field_name.m`
- `get_metadata.m`
- `get_metadata_val.m` 

## Things to further check

- list of BIDS parameters, differentiate anatomical images from the field maps?
- use latest version of `hmri_BIDSify_json.m` not the one from GitLab, then remove the full `acqpar` structure for the demo bit.
- same for `JSONtabl_dcm2bids.tsv`  and `hmri_create_JSONtabl.m`