# READII RUN WITH IDC DATA

## Introduction

This repository contains the code to run the READII pipeline on data available from the [Imaging Data Commons (IDC)](https://imaging.datacommons.cancer.gov/).

The pipeline is designed to extract radiomic features from DICOM images representing CT scans paired with segmentation masks as RTSTRUCTs.;

The python packages used are:

- [readii](https://github.com/bhklab/readii)
- [med-imagetools](https://github.com/bhklab/med-imagetools)

## Notebooks

All notebooks are located in the `notebooks/` directory.

### `notebooks/0_DownloadPreparation.ipynb`

This will download the data from the IDC dataset and prepare it for the next steps

<details>
  <summary>Directory Structure</summary>

  ```console
  notebooks/data/images/nsclc_radiomics//
  ├── dicoms/
  │   ├── Patient-LUNG1-101/
  │   │   └── StudyUID-27911/
  │   │       ├── CT_SeriesUID-55665/
  │   │       └── RTSTRUCT_SeriesUID-25865/
  │   ├── Patient-LUNG1-108/
  │   │   └── StudyUID-62453/
  │   │       ├── CT_SeriesUID-81484/
  │   │       └── RTSTRUCT_SeriesUID-99496/
  │   └── Patient-LUNG1-162/
  │       └── StudyUID-21249/
  │           ├── CT_SeriesUID-72433/
  │           └── RTSTRUCT_SeriesUID-38612/
  └── niftis/
      ├── SubjectID-0_LUNG1-162/
      │   └── StudyUID-21249/
      │       ├── CT_SeriesUID-72433/
      │       │   ├── original.nii.gz
      │       │   ├── randomized_sampled_full.nii.gz
      │       │   ├── randomized_sampled_non_roi.nii.gz
      │       │   ├── randomized_sampled_roi.nii.gz
      │       │   ├── shuffled_full.nii.gz
      │       │   ├── shuffled_non_roi.nii.gz
      │       │   └── shuffled_roi.nii.gz
      │       └── RTSTRUCT_SeriesUID-38612/
      │           └── GTV.nii.gz
      ├── SubjectID-1_LUNG1-101/
      │   └── StudyUID-27911/
      │       ├── CT_SeriesUID-55665/
      │       │   ├── original.nii.gz
      │       │   ├── randomized_sampled_full.nii.gz
      │       │   ├── randomized_sampled_non_roi.nii.gz
      │       │   ├── randomized_sampled_roi.nii.gz
      │       │   ├── shuffled_full.nii.gz
      │       │   ├── shuffled_non_roi.nii.gz
      │       │   └── shuffled_roi.nii.gz
      │       └── RTSTRUCT_SeriesUID-25865/
      │           └── GTV.nii.gz
      └── SubjectID-2_LUNG1-108/
          └── StudyUID-62453/
              ├── CT_SeriesUID-81484/
              │   ├── original.nii.gz
              │   ├── randomized_sampled_full.nii.gz
              │   ├── randomized_sampled_non_roi.nii.gz
              │   ├── randomized_sampled_roi.nii.gz
              │   ├── shuffled_full.nii.gz
              │   ├── shuffled_non_roi.nii.gz
              │   └── shuffled_roi.nii.gz
              └── RTSTRUCT_SeriesUID-99496/
                  └── GTV.nii.gz
  ```

</details>
