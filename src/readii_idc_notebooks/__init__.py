import re
from collections import namedtuple
from itertools import product
from pathlib import Path
from typing import Generator, NamedTuple

import pandas as pd
from imgtools.autopipeline import ImageAutoInput
from imgtools.io import read_image
from imgtools.io.writers.nifti_writer import NiftiWriter
from readii import loaders as rdloaders
from readii.feature_extraction import generateNegativeControl
from tqdm import tqdm

subject_tuple = namedtuple(
    "Subject",
    [
        "PatientID",
        "image",
        "mask",
        "image_id",
        "mask_id",
        "image_modality",
        "mask_modality",
    ],
)


def generate_image_mask_pairs(
    nifti_dir: Path, index_filename: str = "dataset_index.csv"
) -> Generator[NamedTuple, None, None]:
    assert (
        nifti_dir.exists() and nifti_dir.is_dir()
    ), f"{nifti_dir} does not exist or is not a directory"

    index_path = nifti_dir / index_filename
    assert index_path.exists(), f"{index_filename} not found in {nifti_dir}"

    index = pd.read_csv(index_path)

    pairs = []
    grouped = index.groupby("PatientID")
    for _, group in grouped:
        rtstructs = group[group["Modality"] == "RTSTRUCT"].to_dict("records")
        cts = group[group["Modality"] == "CT"].to_dict("records")
        pairs.extend(product(rtstructs, cts))

    total = len(pairs)

    with tqdm(
        total=total,
        bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
    ) as pbar:
        for rt, ct in pairs:
            pbar.update(1)
            pbar.set_description(f"Processing {ct['PatientID']} : {ct['IMAGE_ID']:<30}")
            assert rt["PatientID"] == ct["PatientID"], "PatientID mismatch"
            yield subject_tuple(
                PatientID=ct["PatientID"],
                image=read_image(nifti_dir / ct["filepath"]),
                mask=read_image(nifti_dir / rt["filepath"]),
                image_id=ct["IMAGE_ID"],
                mask_id=rt["IMAGE_ID"],
                image_modality=ct["Modality"],
                mask_modality=rt["Modality"],
            )


def generate_and_save_negative_controls(
    patient: pd.Series,
    roi_match_pattern: dict,
    negative_control_list: list,
    writer: NiftiWriter,
    random_seed: int,
):
    base_image = rdloaders.loadDicomSITK(patient.folder_CT)
    ROI_NAME = list(roi_match_pattern.keys())[0]
    mask_image = rdloaders.loadRTSTRUCTSITK(
        rtstructPath=patient.folder_RTSTRUCT_CT,
        baseImageDirPath=patient.folder_CT,
        roiNames=roi_match_pattern,
    ).get(ROI_NAME)
    writer.save(
        image=base_image,
        PatientID=patient.Index,
        StudyInstanceUID=patient.study[-5:],
        SeriesInstanceUID=patient.series_CT[-5:],
        Modality="CT",
        IMAGE_ID="original",
    )
    writer.save(
        image=mask_image,
        PatientID=patient.Index,
        StudyInstanceUID=patient.study[-5:],
        SeriesInstanceUID=patient.series_RTSTRUCT_CT[-5:],
        Modality="RTSTRUCT",
        IMAGE_ID=ROI_NAME,
    )

    for i, NEGATIVE_CONTROL in enumerate(negative_control_list):
        print(
            f"\tGenerating negative control {i+1}/{len(negative_control_list)} {NEGATIVE_CONTROL}"
        )
        neg_control_image = generateNegativeControl(
            ctImage=base_image,
            alignedROIImage=mask_image,
            randomSeed=random_seed,
            negativeControl=NEGATIVE_CONTROL,
        )
        # Save the negative control image
        writer.save(
            image=neg_control_image,
            PatientID=patient.Index,
            StudyInstanceUID=patient.study[-5:],
            SeriesInstanceUID=patient.series_CT[-5:],
            Modality="CT",
            IMAGE_ID=NEGATIVE_CONTROL,
        )


# "SubjectID-{PatientID}/{Modality}_{SeriesInstanceUID}_{IMAGE_ID}.nii.gz",


def index_and_submit_saves(
    input_dir,
    modalities,
    roi_match_pattern,
    update_imgtools_index,
    n_jobs,
    nifti_output_dir,
    filename_format,
    overwrite,
    random_seed,
    negative_control_list,
):
    neg_nifti_writer = NiftiWriter(
        root_directory=nifti_output_dir,
        filename_format=filename_format,
        create_dirs=True,
        existing_file_mode="overwrite" if overwrite else "fail",
        sanitize_filenames=True,
    )
    dataset = ImageAutoInput(
        dir_path=input_dir,
        modalities=",".join(modalities),
        update=update_imgtools_index,
        n_jobs=n_jobs,
    )

    for p, patient in enumerate(dataset.df_combined.itertuples()):
        print(
            f"Loading data for subject_ID {p+1}/{len(dataset.df_combined)}: {patient.Index} (PatientID : {patient.patient_ID})"
        )

        generate_and_save_negative_controls(
            patient=patient,
            roi_match_pattern=roi_match_pattern,
            negative_control_list=negative_control_list,
            writer=neg_nifti_writer,
            random_seed=random_seed,
        )
    filename_pattern = neg_nifti_writer.pattern_resolver.formatted_pattern.replace(
        "%(", "(?P<"
    ).replace(")s", ">.*?)")

    datafiles = []
    for file_path in nifti_output_dir.rglob("*.nii.gz"):
        if match := re.search(filename_pattern, str(file_path).replace("\\", "/")):
            relative_path = file_path.absolute().relative_to(
                nifti_output_dir.absolute()
            )
            datafiles.append({**match.groupdict(), "filepath": relative_path})
    datafiles_df = pd.DataFrame(datafiles)
    csv_path = nifti_output_dir / "dataset_index.csv"
    # sort on PatientID and IMAGE_ID
    datafiles_df = datafiles_df.sort_values(
        by=["PatientID", "IMAGE_ID"], ignore_index=True
    )

    datafiles_df.to_csv(csv_path, index=False)
    return csv_path
