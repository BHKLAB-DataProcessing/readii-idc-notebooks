```python
from idc_index import index
import ipywidgets as widgets
from pathlib import Path
from pydicom import dcmread
from tqdm.notebook import tqdm
import subprocess
from imgtools.dicom.sort import DICOMSorter
from imgtools.logging import logger as imgtools_logger

import tempfile
imgtools_logger.setLevel("WARNING")
```


```python
client = index.IDCClient()
print(client.get_idc_version())
```

    v17



```python
collections = sorted(client.get_collections())
print(f"Found {len(collections)} collections")  
client.collection_summary
```

    Found 142 collections





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Modality</th>
      <th>series_size_MB</th>
    </tr>
    <tr>
      <th>collection_id</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>4d_lung</th>
      <td>[RTSTRUCT, CT]</td>
      <td>183054.14</td>
    </tr>
    <tr>
      <th>acrin_6698</th>
      <td>[MR, SEG]</td>
      <td>841956.27</td>
    </tr>
    <tr>
      <th>acrin_contralateral_breast_mr</th>
      <td>[MR, CR]</td>
      <td>199592.57</td>
    </tr>
    <tr>
      <th>acrin_flt_breast</th>
      <td>[PT, CT, OT]</td>
      <td>74235.64</td>
    </tr>
    <tr>
      <th>acrin_nsclc_fdg_pet</th>
      <td>[PT, CT, DX, CR, NM, MR, SC, SEG]</td>
      <td>145677.88</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>tcga_uvm</th>
      <td>[SM]</td>
      <td>102250.70</td>
    </tr>
    <tr>
      <th>upenn_gbm</th>
      <td>[MR]</td>
      <td>139399.35</td>
    </tr>
    <tr>
      <th>vestibular_schwannoma_mc_rc</th>
      <td>[MR]</td>
      <td>10015.97</td>
    </tr>
    <tr>
      <th>vestibular_schwannoma_seg</th>
      <td>[RTDOSE, RTSTRUCT, RTPLAN, MR]</td>
      <td>28194.20</td>
    </tr>
    <tr>
      <th>victre</th>
      <td>[MG]</td>
      <td>1027187.60</td>
    </tr>
  </tbody>
</table>
<p>142 rows × 2 columns</p>
</div>



# Filter Collections have both RTSTRUCTS and CTs


```python
# Group by collection_id and aggregate the Modality column
grouped = client.index.groupby('collection_id')['Modality'].apply(lambda x: set(x))

# Filter collection_ids where both RTSTRUCT and CT exist in the Modality set
rt_ct_collections = grouped[grouped.apply(lambda x: 'RTSTRUCT' in x and 'CT' in x)]

print(f"Found {len(rt_ct_collections)} collections with both RTSTRUCT and CT")
rt_ct_collections
```

    Found 13 collections with both RTSTRUCT and CT





    collection_id
    4d_lung                                           {CT, RTSTRUCT}
    cc_tumor_heterogeneity               {MR, CT, PT, REG, RTSTRUCT}
    cptac_ccrcc                               {CT, RTSTRUCT, SM, MR}
    cptac_pda                         {MR, CT, PT, SM, US, RTSTRUCT}
    cptac_ucec                            {CT, MR, PT, SM, RTSTRUCT}
    lctsc                                             {CT, RTSTRUCT}
    nsclc_radiomics                          {CT, SEG, SR, RTSTRUCT}
    nsclc_radiomics_interobserver1               {CT, SEG, RTSTRUCT}
    pancreatic_ct_cbct_seg                    {CT, RTDOSE, RTSTRUCT}
    pediatric_ct_seg                                  {CT, RTSTRUCT}
    prostate_anatomical_edge_cases                    {CT, RTSTRUCT}
    rider_lung_ct                        {CT, PR, SR, SEG, RTSTRUCT}
    soft_tissue_sarcoma                       {CT, MR, RTSTRUCT, PT}
    Name: Modality, dtype: object



# Select a Collection to explore from the dropdown menu below.

The default is `nsclc_radiomics`.


```python
collection_list = rt_ct_collections.index

collection_widget = widgets.Dropdown(
    options=collection_list,
    description='Collection:',
    value='nsclc_radiomics',
    disabled=False,
)
display(collection_widget)
```


    Dropdown(description='Collection:', index=6, options=('4d_lung', 'cc_tumor_heterogeneity', 'cptac_ccrcc', 'cpt…



```python
collection_id = collection_widget.value 

matching_series = client.index.loc[client.index.collection_id == collection_id, ['SeriesInstanceUID', 'Modality', "series_size_MB"]]
print(f"Found {len(matching_series)} series in collection {collection_id}")
options=[
  (
    f'SeriesUID-{row['SeriesInstanceUID'][-10:]} [Modality: {row["Modality"]}; Size: {row["series_size_MB"]}MB]',
    row["SeriesInstanceUID"]
  )
  for _, row in matching_series.iterrows()
  if row['Modality'] == 'RTSTRUCT'
]
print(f"Found {len(options)} RTSTRUCT series in collection {collection_id}")
rt_widget = widgets.Dropdown(
    options=options,
    description='Series:',
    layout={'width': 'max-content'},
    disabled=False,
)
display(rt_widget)
```

    Found 4926 series in collection nsclc_radiomics
    Found 422 RTSTRUCT series in collection nsclc_radiomics



    Dropdown(description='Series:', layout=Layout(width='max-content'), options=(('SeriesUID-9499301975 [Modality:…



```python
print(f"Selected RTSTRUCT Series: {rt_widget.value}")
```

    Selected RTSTRUCT Series: 1.3.6.1.4.1.32722.99.99.330349836421659898040510720189499301975


## Setting up Directories To Download



```python
# Find user's s5cmd path
s5cmd = client.s5cmdPath


# Save data to local directory
DATA_DIR = Path('data') 

# Create a temporary directory to store the downloaded files
TMP_DIR = Path(tempfile.mkdtemp())
TMP_DIR.mkdir(parents=True, exist_ok=True)


```


```python
rt_url = client.get_series_file_URLs(rt_widget.value)[0]
print(rt_url)
```

    s3://idc-open-data-cr/b113da7c-c8ef-4fcb-a6f8-a10e76f422f7/c7eef1b8-e04d-4fbb-b10f-6b421a23afe6.dcm


## Download chosen RTSTRUCT and CT files

1. download the RTSTRUCT
2. Query the RTSTRUCT's metadata for the CT `SeriesInstanceUID` it references
3. Download the CT files corresponding to the `SeriesInstanceUID`


```python
print("Downloading RTSTRUCT file...")
! $s5cmd --no-sign-request --endpoint-url https://s3.amazonaws.com cp --show-progress  $rt_url $TMP_DIR
rt_path = TMP_DIR.iterdir().__next__()
print(f"RTSTRUCT file downloaded to {rt_path}")
```

    Downloading RTSTRUCT file...


    [32m43.23%[0m [32m ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━──────────────────────────────────────────────── [0m [32m730.64 kB / 1.69 MB[0m [31m? p/s[0m [34m?[0m [33m(0/1)[0m

    [32m100.00%[0m [32m ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ [0m [32m1.69 MB / 1.69 MB[0m [31m(44.42 MB/s)[0m [34m200ms[0m [33m(1/1)[0m


    RTSTRUCT file downloaded to /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/c7eef1b8-e04d-4fbb-b10f-6b421a23afe6.dcm



```python
ds = dcmread(rt_path, stop_before_pixels=True, specific_tags=['ReferencedFrameOfReferenceSequence', 'StructureSetROISequence'])

referenced_ct = ds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].SeriesInstanceUID
print(f"Referenced CT Series: {referenced_ct}")
```

    Referenced CT Series: 1.3.6.1.4.1.32722.99.99.238725422689739575345635473736072968163



```python
ct_urls = client.get_series_file_URLs(referenced_ct)
for ct_url in tqdm(ct_urls, desc="Downloading CT files"):
  subprocess.run([s5cmd, "--no-sign-request", "--endpoint-url", "https://s3.amazonaws.com", "sync", ct_url, TMP_DIR])

```


    Downloading CT files:   0%|          | 0/134 [00:00<?, ?it/s]


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/004ac528-0bfc-4b17-a343-60c43b797f29.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/004ac528-0bfc-4b17-a343-60c43b797f29.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/01c8faae-d4fe-4d41-91cd-9e5862befb12.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/01c8faae-d4fe-4d41-91cd-9e5862befb12.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/051a5456-712d-4933-94c8-bae5c3bda4b8.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/051a5456-712d-4933-94c8-bae5c3bda4b8.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/053ee01a-5d70-454a-aed8-0f697686977c.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/053ee01a-5d70-454a-aed8-0f697686977c.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/0a5966ee-6e2d-4679-94e4-31f031e9ae8a.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/0a5966ee-6e2d-4679-94e4-31f031e9ae8a.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/0b3deafb-0b6f-48d0-93a0-2da76bbf2ff8.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/0b3deafb-0b6f-48d0-93a0-2da76bbf2ff8.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/0b563394-7579-4d91-b433-c35ec3cdb133.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/0b563394-7579-4d91-b433-c35ec3cdb133.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/1202c98d-5ab5-48a2-94ab-9fad0a4b44f0.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/1202c98d-5ab5-48a2-94ab-9fad0a4b44f0.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/155b1518-315e-4c9a-bb8d-c64de1aab82d.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/155b1518-315e-4c9a-bb8d-c64de1aab82d.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/1ddddc05-3341-4103-83dc-08cf8b725dc7.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/1ddddc05-3341-4103-83dc-08cf8b725dc7.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/1ded0b07-a566-4987-9d86-35b4837b9b2d.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/1ded0b07-a566-4987-9d86-35b4837b9b2d.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/1e83a2ef-c8fb-4a51-b7d6-4de180ec671e.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/1e83a2ef-c8fb-4a51-b7d6-4de180ec671e.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/1f10fb48-96d2-486d-b040-32dd881b07c3.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/1f10fb48-96d2-486d-b040-32dd881b07c3.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/2038f806-7e34-4e68-9256-2d363f9891b5.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/2038f806-7e34-4e68-9256-2d363f9891b5.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/20d1b8b6-5a47-4352-a5ff-1d56009c2787.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/20d1b8b6-5a47-4352-a5ff-1d56009c2787.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/21175ee4-52a6-4b98-a2cc-fb745f00a6a1.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/21175ee4-52a6-4b98-a2cc-fb745f00a6a1.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/22ac3131-b5ef-4142-b09f-fa2439eef180.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/22ac3131-b5ef-4142-b09f-fa2439eef180.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/22dd6afc-9d33-4beb-aa28-a4020f220a3d.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/22dd6afc-9d33-4beb-aa28-a4020f220a3d.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/256a4083-26ff-46b2-a8f6-c31ecdd68cad.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/256a4083-26ff-46b2-a8f6-c31ecdd68cad.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/259194d1-3dfe-43fd-a8c0-d920c388b9d2.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/259194d1-3dfe-43fd-a8c0-d920c388b9d2.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/2817d0e4-a159-4e77-986d-37c731af0120.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/2817d0e4-a159-4e77-986d-37c731af0120.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/29d0a855-38dc-47e4-b7a4-d49cf28c6b7c.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/29d0a855-38dc-47e4-b7a4-d49cf28c6b7c.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/2f1bd7e1-78b6-4737-9b40-26f26c50ba2a.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/2f1bd7e1-78b6-4737-9b40-26f26c50ba2a.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/30cdfcbd-7a68-48c6-b468-89fafb336b3b.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/30cdfcbd-7a68-48c6-b468-89fafb336b3b.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/336019a3-d697-4168-ae6b-1e7fbcadf9a7.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/336019a3-d697-4168-ae6b-1e7fbcadf9a7.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/344dee4f-4e8a-4372-bef1-400cc55b8127.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/344dee4f-4e8a-4372-bef1-400cc55b8127.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/35b3f25c-58cc-499d-ab57-4cb9c3ebd17e.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/35b3f25c-58cc-499d-ab57-4cb9c3ebd17e.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/37d1fa69-0f9e-4315-a71d-9606ffbe2d07.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/37d1fa69-0f9e-4315-a71d-9606ffbe2d07.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/37d7006c-7c7d-4041-9b10-c8e2b3b40f89.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/37d7006c-7c7d-4041-9b10-c8e2b3b40f89.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/3e3cdfeb-4e4f-4219-8d73-fada7502c317.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/3e3cdfeb-4e4f-4219-8d73-fada7502c317.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/3e4d285f-ab13-430c-9e46-8014cdfc5586.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/3e4d285f-ab13-430c-9e46-8014cdfc5586.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/3ecacd56-7e8a-46dc-9c67-1f8714d99467.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/3ecacd56-7e8a-46dc-9c67-1f8714d99467.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/3ed07449-f252-48d3-a426-1e0d7241fcb3.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/3ed07449-f252-48d3-a426-1e0d7241fcb3.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/437b0503-f8c9-4352-a3d5-9ca5b4cce1cf.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/437b0503-f8c9-4352-a3d5-9ca5b4cce1cf.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/46280344-1db6-48ef-bc29-d2655d071a17.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/46280344-1db6-48ef-bc29-d2655d071a17.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/4773c345-7c13-4b2a-be3d-2cd58ec4e8e3.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/4773c345-7c13-4b2a-be3d-2cd58ec4e8e3.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/4949a818-429d-47a9-978b-f92a4b816813.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/4949a818-429d-47a9-978b-f92a4b816813.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/4b0fc4e4-ed4e-4566-ada8-7eb772ad717d.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/4b0fc4e4-ed4e-4566-ada8-7eb772ad717d.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/4b5e714c-8013-408f-b90d-7f8daaabe81d.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/4b5e714c-8013-408f-b90d-7f8daaabe81d.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/4c2ee72c-0aad-4c61-8ca6-a8c603e2ec19.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/4c2ee72c-0aad-4c61-8ca6-a8c603e2ec19.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/4cdcdee8-5f7a-4d9e-a734-9ffab1d8b04c.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/4cdcdee8-5f7a-4d9e-a734-9ffab1d8b04c.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/4da34c04-d46c-4f80-b336-4fd46764fa00.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/4da34c04-d46c-4f80-b336-4fd46764fa00.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/4e29bf3a-5463-4ddb-9af5-b0e41b9c6ef0.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/4e29bf3a-5463-4ddb-9af5-b0e41b9c6ef0.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/4f9a4490-3051-49a2-b321-2c76694d5729.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/4f9a4490-3051-49a2-b321-2c76694d5729.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/5407ee7d-1512-476d-b26d-24b6aa91f76f.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/5407ee7d-1512-476d-b26d-24b6aa91f76f.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/545c883d-3cca-4ff9-a21f-c2f51c4f082a.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/545c883d-3cca-4ff9-a21f-c2f51c4f082a.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/5f17d22c-0b4b-4889-821c-e53818eb4b9c.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/5f17d22c-0b4b-4889-821c-e53818eb4b9c.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/6002bfb0-7d59-4309-bd8c-38ffbd0054ea.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/6002bfb0-7d59-4309-bd8c-38ffbd0054ea.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/6268dc8e-0328-492b-ac28-fd35bd961aba.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/6268dc8e-0328-492b-ac28-fd35bd961aba.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/65c470a0-09d2-4537-a823-0224181fbe6c.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/65c470a0-09d2-4537-a823-0224181fbe6c.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/699262f1-cbf7-4f61-b70f-6802c4fd7d34.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/699262f1-cbf7-4f61-b70f-6802c4fd7d34.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/69df9597-1310-4900-a578-0f5231c889c9.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/69df9597-1310-4900-a578-0f5231c889c9.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/6a7e0be2-b2c8-4b8a-b6d7-11004e90e8d4.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/6a7e0be2-b2c8-4b8a-b6d7-11004e90e8d4.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/748483ff-721b-4ebb-8c31-87577aa56f41.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/748483ff-721b-4ebb-8c31-87577aa56f41.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/75727310-5357-4d16-9c0a-d01203a76ecd.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/75727310-5357-4d16-9c0a-d01203a76ecd.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/75ef1df0-894e-4412-89c6-a60501f2e625.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/75ef1df0-894e-4412-89c6-a60501f2e625.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/76d1de9d-3ed1-4d30-ba38-4d1d830d6875.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/76d1de9d-3ed1-4d30-ba38-4d1d830d6875.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/7df6d390-d0cf-4d3d-ba2b-c8c246678070.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/7df6d390-d0cf-4d3d-ba2b-c8c246678070.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/7e4586cc-1234-462d-b9b0-55b9af5c1205.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/7e4586cc-1234-462d-b9b0-55b9af5c1205.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/80122308-4c24-44a2-977a-c24439163e52.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/80122308-4c24-44a2-977a-c24439163e52.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/82049c09-b518-4f38-b84f-1eb971810179.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/82049c09-b518-4f38-b84f-1eb971810179.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/82145482-24cc-4868-8f6a-de8873c7f43d.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/82145482-24cc-4868-8f6a-de8873c7f43d.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/82405f13-98af-428a-9266-8054ed540de0.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/82405f13-98af-428a-9266-8054ed540de0.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/825420d5-4c74-48e4-b9de-0d31aabb4f09.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/825420d5-4c74-48e4-b9de-0d31aabb4f09.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/86dc2441-de30-4129-a4f7-7b942bcb4236.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/86dc2441-de30-4129-a4f7-7b942bcb4236.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/86ee86e6-4844-4e78-9f58-10b75f6aee83.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/86ee86e6-4844-4e78-9f58-10b75f6aee83.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/89c6c545-fa5e-4080-9023-262ee1fc5ee1.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/89c6c545-fa5e-4080-9023-262ee1fc5ee1.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/8a9f5241-5d6c-41ac-bb9d-e59bf5e89092.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/8a9f5241-5d6c-41ac-bb9d-e59bf5e89092.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/8adc4ada-3b54-444a-9236-486e57e8917b.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/8adc4ada-3b54-444a-9236-486e57e8917b.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/8bc9ef2d-d24a-4055-90cd-9d3d7684c9a1.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/8bc9ef2d-d24a-4055-90cd-9d3d7684c9a1.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/8bdf9df7-091e-43eb-b8e0-5b95cdc9f525.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/8bdf9df7-091e-43eb-b8e0-5b95cdc9f525.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/8d11659b-53bf-4e98-99a5-0a46936490b1.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/8d11659b-53bf-4e98-99a5-0a46936490b1.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/8f049afb-c608-48d6-ad08-549615e99b10.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/8f049afb-c608-48d6-ad08-549615e99b10.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/955c9f3f-231e-4f6b-9261-86570bd7695c.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/955c9f3f-231e-4f6b-9261-86570bd7695c.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/9881bf31-924b-424b-885f-afc7e0f6ca33.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/9881bf31-924b-424b-885f-afc7e0f6ca33.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/99576d4f-d990-41d2-8c72-aefb02d4b6b4.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/99576d4f-d990-41d2-8c72-aefb02d4b6b4.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/99efad56-527a-456e-aa75-58c5d9325b07.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/99efad56-527a-456e-aa75-58c5d9325b07.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/9a943cbf-cda0-456f-bee5-0a4210b281a4.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/9a943cbf-cda0-456f-bee5-0a4210b281a4.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/9b1bb1d5-17eb-43b8-ba60-8367bf59e7be.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/9b1bb1d5-17eb-43b8-ba60-8367bf59e7be.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/9b24f980-97f4-41b0-b8ab-a29c34715992.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/9b24f980-97f4-41b0-b8ab-a29c34715992.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/9c1b37a3-1903-4ec9-bda6-d5bb1a08297d.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/9c1b37a3-1903-4ec9-bda6-d5bb1a08297d.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/9c8349a8-273c-4b23-b71d-8c1a82894211.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/9c8349a8-273c-4b23-b71d-8c1a82894211.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/9caa0265-0188-4732-ac0c-49fc5bc40b1d.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/9caa0265-0188-4732-ac0c-49fc5bc40b1d.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/9ee1715f-cb77-407d-acb5-196013160630.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/9ee1715f-cb77-407d-acb5-196013160630.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/a1aa09c2-bb73-4644-af2e-22da710f35b0.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/a1aa09c2-bb73-4644-af2e-22da710f35b0.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/a4d74efb-c2d6-4313-8ba3-0ac590072db9.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/a4d74efb-c2d6-4313-8ba3-0ac590072db9.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/a6e19055-0027-4557-8c23-8edd5240f306.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/a6e19055-0027-4557-8c23-8edd5240f306.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/a72f94f7-3424-4df4-a70f-ea8b4470a8d8.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/a72f94f7-3424-4df4-a70f-ea8b4470a8d8.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/ad11f401-3526-4b32-838e-dd10f783876b.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/ad11f401-3526-4b32-838e-dd10f783876b.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/b0b72889-1648-4a0d-9e6a-e815218969f8.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/b0b72889-1648-4a0d-9e6a-e815218969f8.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/b3290761-bf13-40a9-b394-3f9d85c29763.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/b3290761-bf13-40a9-b394-3f9d85c29763.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/b54e219a-add9-45c2-b106-f9d77f35c7f7.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/b54e219a-add9-45c2-b106-f9d77f35c7f7.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/b570f564-d468-4eff-9038-c7036cad66c1.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/b570f564-d468-4eff-9038-c7036cad66c1.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/b5c5f2f0-fa65-4029-b67c-9b5d55cbb822.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/b5c5f2f0-fa65-4029-b67c-9b5d55cbb822.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/b860f950-bd8c-4802-b579-380043f81ba8.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/b860f950-bd8c-4802-b579-380043f81ba8.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/b8797b82-8cc7-4674-a7ea-57529e4a1d0e.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/b8797b82-8cc7-4674-a7ea-57529e4a1d0e.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/b8dc7029-7f41-46cc-825a-3eb91070ca32.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/b8dc7029-7f41-46cc-825a-3eb91070ca32.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/bc8e22a9-da2a-47f1-b157-a2ca5b882757.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/bc8e22a9-da2a-47f1-b157-a2ca5b882757.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/bd83919a-a094-403f-8735-586571a5555a.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/bd83919a-a094-403f-8735-586571a5555a.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/bdb21e8f-c364-4730-9886-a53c4220fb00.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/bdb21e8f-c364-4730-9886-a53c4220fb00.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/be28e443-1cbf-40bf-ad55-bce56e21ca56.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/be28e443-1cbf-40bf-ad55-bce56e21ca56.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/bf313887-c738-4dd3-9f37-863ef67bb1e5.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/bf313887-c738-4dd3-9f37-863ef67bb1e5.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/bfdf30b6-f685-487a-925b-86443c90383f.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/bfdf30b6-f685-487a-925b-86443c90383f.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/c31be54a-b044-4949-bdb9-2b378133e673.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/c31be54a-b044-4949-bdb9-2b378133e673.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/c390caf2-4189-41ee-83ef-98b5be764e30.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/c390caf2-4189-41ee-83ef-98b5be764e30.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/c4c5309c-8632-412b-831e-8652318c5c85.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/c4c5309c-8632-412b-831e-8652318c5c85.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/c5496741-c20f-456d-ba02-d75a1de18a78.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/c5496741-c20f-456d-ba02-d75a1de18a78.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/c8ad46b6-8651-4df7-aaba-04ba5afde3af.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/c8ad46b6-8651-4df7-aaba-04ba5afde3af.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/c9f50185-6f91-439e-8db8-5c5ec8b7dd17.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/c9f50185-6f91-439e-8db8-5c5ec8b7dd17.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/ca0db255-2a7a-40cd-a1ae-99d57c7bf9ae.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/ca0db255-2a7a-40cd-a1ae-99d57c7bf9ae.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/cb8e85a8-31f4-4575-b2e9-9b9ee6297b3e.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/cb8e85a8-31f4-4575-b2e9-9b9ee6297b3e.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/d07647ca-ea1e-4314-bd6d-58ac8f35f2a3.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/d07647ca-ea1e-4314-bd6d-58ac8f35f2a3.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/d2bda643-ead9-4ce6-8967-31d295ed5cdd.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/d2bda643-ead9-4ce6-8967-31d295ed5cdd.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/d5c8f0e2-164e-4408-b7f8-e4f1a5fca5c8.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/d5c8f0e2-164e-4408-b7f8-e4f1a5fca5c8.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/d67663ef-9582-4660-9b78-c83fe859b428.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/d67663ef-9582-4660-9b78-c83fe859b428.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/d950cff1-c682-4874-bd61-8fa50321e962.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/d950cff1-c682-4874-bd61-8fa50321e962.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/da0b1fe2-82dc-4583-9705-3dcf34984d12.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/da0b1fe2-82dc-4583-9705-3dcf34984d12.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/daefaebc-ed23-4d41-a2d2-a5c724b5cbf5.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/daefaebc-ed23-4d41-a2d2-a5c724b5cbf5.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/dd03f444-b458-4189-bf79-9943d81b2a3e.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/dd03f444-b458-4189-bf79-9943d81b2a3e.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/dd3ad780-fe63-4801-93ef-adcee56e3c8e.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/dd3ad780-fe63-4801-93ef-adcee56e3c8e.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/dde1e8c9-5468-4bef-9580-7eae9efae2f7.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/dde1e8c9-5468-4bef-9580-7eae9efae2f7.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/dea32744-1a2a-4576-b3ca-979abac21f8f.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/dea32744-1a2a-4576-b3ca-979abac21f8f.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/e065a60e-2660-4b76-94a7-a13de8a16db9.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/e065a60e-2660-4b76-94a7-a13de8a16db9.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/e333047d-9adc-4bd5-b9ee-82c61eddc64c.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/e333047d-9adc-4bd5-b9ee-82c61eddc64c.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/e611f248-578e-48c8-89d2-109060172e9b.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/e611f248-578e-48c8-89d2-109060172e9b.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/ea150852-911d-4e1f-ba14-0225981748a5.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/ea150852-911d-4e1f-ba14-0225981748a5.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/ea999eea-624a-42bc-bc56-265bd9d906de.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/ea999eea-624a-42bc-bc56-265bd9d906de.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/eff4180a-392a-4675-ac5d-b956fbce61c5.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/eff4180a-392a-4675-ac5d-b956fbce61c5.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/f25b0ef2-afa0-40e1-8e76-cd4a6310843f.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/f25b0ef2-afa0-40e1-8e76-cd4a6310843f.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/f473cd3c-11ec-4a4c-83c0-4ee7da4f1504.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/f473cd3c-11ec-4a4c-83c0-4ee7da4f1504.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/f52f5b43-c6fc-4f77-9497-3941b89ae945.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/f52f5b43-c6fc-4f77-9497-3941b89ae945.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/f79c8d77-3073-4bb2-ab0f-069092453ae5.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/f79c8d77-3073-4bb2-ab0f-069092453ae5.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/f821bf0d-f131-4eb9-a410-810a94b3cd9f.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/f821bf0d-f131-4eb9-a410-810a94b3cd9f.dcm


    cp s3://idc-open-data-cr/43263a91-f1ae-4041-b9b4-c05b8a221f8a/fdbf99c0-9013-4431-8ad3-9d8ffd92ebf4.dcm /var/folders/8t/rwh6rzg93jxfqkb63gt2n4940000gn/T/tmpgg1p00fn/fdbf99c0-9013-4431-8ad3-9d8ffd92ebf4.dcm


# Sort the dicom files into an appropriate structure

The dicom files are all named with a unique UUID. 
This makes it difficult to understand which files are related to each other.

We will sort the files into a directory structure that makes it easier to understand the relationships between the files.

This uses `Med-ImageTools`' `DICOMSorter` class to sort the files into a directory structure.

The structure we are aiming for is:

```console
./data/<collection_name>/dicoms/sorted/
└── Patient-<PatientID>
    └── StudyUID-<StudyInstanceUID>
        ├── <Modality>_SeriesUID-<SeriesInstanceUID>
        └── <Modality2>_SeriesUID-<SeriesInstanceUID>
            ├── DICOM-FILE
            └── DICOM-FILE
```

**Note:**
Earlier, we downloaded the data to a temporary directory, so we will perform a `move` operation on the sorter
If you do not want to move your input data, you can use the `symlink` option to create symbolic links to the files instead of moving them.


```python
sorted_path = DATA_DIR / "images" / collection_id / "dicoms"

if sorted_path.exists():
  raise FileExistsError(f"Sorted Path Exists, please remove `{sorted_path}`")

dicomsorter = DICOMSorter(
  source_directory=TMP_DIR.absolute(),
  target_pattern=Path(
    sorted_path,
    "Patient-%PatientID/StudyUID-%StudyInstanceUID/%Modality_SeriesUID-%SeriesInstanceUID/"
  ).as_posix(),
)
dicomsorter.execute(action="move")

```


    Output()



<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"></pre>




    Output()



<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"></pre>




```python
print("New Directory Structure: ")
subprocess.run(["tree", "-d", sorted_path.absolute()])
```

    New Directory Structure: 
    [01;34m/Users/bhklab/dev/radiomics/readii-idc-notebooks/notebooks/data/images/nsclc_radiomics/dicoms[0m
    └── [01;34mPatient-LUNG1-303[0m
        └── [01;34mStudyUID-60013[0m
            ├── [01;34mCT_SeriesUID-68163[0m
            └── [01;34mRTSTRUCT_SeriesUID-01975[0m
    
    5 directories





    CompletedProcess(args=['tree', '-d', PosixPath('/Users/bhklab/dev/radiomics/readii-idc-notebooks/notebooks/data/images/nsclc_radiomics/dicoms')], returncode=0)




```python

```
