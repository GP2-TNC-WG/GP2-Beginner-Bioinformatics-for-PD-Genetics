## Introduction. Using Terra to Access Data and Perform Analyses

### Information
- **Created by:** GP2 Training and Networking
- **Don't forget** to go through the official AMP-PD 'Getting Started' notebooks for further information. 
    - [Release 1](https://app.terra.bio/#workspaces/amp-pd-release-1-stage/AMP%20PD%20Release%201%20Getting%20Started)
    - [Tier 1](https://app.terra.bio/#workspaces/amp-pd-release-v1/Getting%20Started%20Tier%201%20-%20Clinical%20Access) 
    - [Tier 2](https://app.terra.bio/#workspaces/amp-pd-release-v1/Getting%20Started%20Tier%202%20-%20Clinical%20and%20Omics%20Access)

## Table of Contents

#### [0. Using a Notebook](#0)

#### [1. Getting Started](#1)

#### [2. Querying Clinical Data](#2)

#### [3. Querying Variant Data](#3)

#### [4. Pulling data from Google Storage Buckets](#4)

#### [5. Quick Association Test with Plink](#5)


<a id="0"></a>
## 0. Using a Notebook

- If this is your first time running a Terra Notebook, you should run either the ```Py3 - Start Here``` or the ```R - Start Here``` depending on what language you want to use.



- Cells come in 2 types, markdown cells like this one and code cells like the cell below

    - To change the cell type, click on ```Cell``` in the top left and scroll to ```Cell Type```.



- Hit the ```Run``` button in the top left to move through the cells running the code sequentially as you go.


- You can also turn cells into bash cells to use commands like you would use in your terminal by using "%%bash"


```python
print(3 + 3)
```

    6



```python
print(20 * 2)
```

    40



```bash
%%bash

pwd
```

    /home/jupyter-user/notebooks/AMP PD WGS QC Collaboration/edit


<a id="1"></a>
## 1. Getting Started

- Run the cells below to set up libraries and functions.


- We will also set up the billing project for making queries and set a few path variables to pull data.

### Set up libraries


```python
# Use the os package to interact with the environment
import os

# Bring in Pandas for Dataframe functionality
import pandas as pd

# numpy for basics
import numpy as np

# Use StringIO for working with file contents
from io import StringIO

# Enable IPython to display matplotlib graphs
import matplotlib.pyplot as plt
%matplotlib inline

# Enable interaction with the FireCloud API
from firecloud import api as fapi

# Import the iPython HTML rendering for displaying links to Google Cloud Console
from IPython.core.display import display, HTML

# Import urllib modules for building URLs to Google Cloud Console
import urllib.parse

# BigQuery for querying data
from google.cloud import bigquery
```

### Set up billing project and data path variables


```python
BILLING_PROJECT_ID = os.environ['GOOGLE_PROJECT']
WORKSPACE_NAMESPACE = os.environ['WORKSPACE_NAMESPACE']
WORKSPACE_NAME = os.environ['WORKSPACE_NAME']
WORKSPACE_BUCKET = os.environ['WORKSPACE_BUCKET']

WORKSPACE_ATTRIBUTES = fapi.get_workspace(WORKSPACE_NAMESPACE, WORKSPACE_NAME).json().get('workspace',{}).get('attributes',{})

GS_RELEASE_PATH = 'gs://amp-pd-data/releases/2019_v1release_1015'
GS_CLINICAL_RELEASE_PATH = f'{GS_RELEASE_PATH}/clinical'

GS_WGS_RELEASE_PATH = 'gs://amp-pd-genomics/releases/2019_v1release_1015/wgs'
GS_WGS_RELEASE_PLINK_PATH = os.path.join(GS_WGS_RELEASE_PATH, 'plink')
GS_WGS_RELEASE_GATK_PATH = os.path.join(GS_WGS_RELEASE_PATH, 'gatk')

BQ_RELEASE_DATASET = 'amp-pd-research.2019_v1release_1015'


print(BILLING_PROJECT_ID)
print(GS_CLINICAL_RELEASE_PATH)
print(GS_WGS_RELEASE_PLINK_PATH)
print(GS_WGS_RELEASE_GATK_PATH)
```

    fc-amp-pd-alpha
    gs://amp-pd-data/releases/2019_v1release_1015/clinical
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk


### Some useful functions


```python
# Utility routine for printing a shell command before executing it
def shell_do(command):
    print(f'Executing: {command}')
    !$command

# Utility routines for reading files from Google Cloud Storage
def gcs_read_file(path):
    """Return the contents of a file in GCS"""
    contents = !gsutil -u {BILLING_PROJECT_ID} cat {path}
    return '\n'.join(contents)
    
def gcs_read_csv(path, sep=None):
    """Return a DataFrame from the contents of a delimited file in GCS"""
    return pd.read_csv(StringIO(gcs_read_file(path)), sep=sep, engine='python')

# Utility routine for display a message and a link
def display_html_link(description, link_text, url):
    html = f'''
    <p>
    </p>
    <p>
    {description}
    <a target=_blank href="{url}">{link_text}</a>.
    </p>
    '''

    display(HTML(html))
    
# Utility routine for displaying a message and link to Cloud Console
def link_to_cloud_console_gcs(description, link_text, gcs_path):
    url = '{}?{}'.format(
        os.path.join('https://console.cloud.google.com/storage/browser',
                     gcs_path.replace("gs://","")),
        urllib.parse.urlencode({'userProject': BILLING_PROJECT_ID}))

    display_html_link(description, link_text, url)
    
# Get the data from a query
def bq_query(query):
    """Return the contents of a query against BigQuery"""
    return pd.read_gbq(
        query,
        project_id=BILLING_PROJECT_ID,
        dialect='standard')
```

<a id="2"></a>
## 2. Querying Clinical Data

Let's look at all the available clinical tables in BigQuery. The following query will tell you the available dataset ID's, names, row count, and size.


```python
clinical_tables = f"""

SELECT 
project_id, dataset_id, table_id, row_count, size_bytes 

FROM `{BQ_RELEASE_DATASET}.__TABLES__`

"""


bq_query(clinical_tables)
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>project_id</th>
      <th>dataset_id</th>
      <th>table_id</th>
      <th>row_count</th>
      <th>size_bytes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>Biospecimen_analyses_CSF_abeta_tau_ptau</td>
      <td>9385</td>
      <td>458495</td>
    </tr>
    <tr>
      <th>1</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>Biospecimen_analyses_CSF_beta_glucocerebrosidase</td>
      <td>278</td>
      <td>24342</td>
    </tr>
    <tr>
      <th>2</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>Biospecimen_analyses_SomaLogic_plasma</td>
      <td>26100</td>
      <td>2273980</td>
    </tr>
    <tr>
      <th>3</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>Biospecimen_analyses_other</td>
      <td>6105</td>
      <td>559618</td>
    </tr>
    <tr>
      <th>4</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>Caffeine_history</td>
      <td>1222</td>
      <td>55586</td>
    </tr>
    <tr>
      <th>5</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>DTI</td>
      <td>1180</td>
      <td>132014</td>
    </tr>
    <tr>
      <th>6</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>DaTSCAN_SBR</td>
      <td>1976</td>
      <td>106071</td>
    </tr>
    <tr>
      <th>7</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>DaTSCAN_visual_interpretation</td>
      <td>842</td>
      <td>27048</td>
    </tr>
    <tr>
      <th>8</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>Demographics</td>
      <td>4298</td>
      <td>413191</td>
    </tr>
    <tr>
      <th>9</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>Enrollment</td>
      <td>4280</td>
      <td>447363</td>
    </tr>
    <tr>
      <th>10</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>Epworth_Sleepiness_Scale</td>
      <td>9353</td>
      <td>2614392</td>
    </tr>
    <tr>
      <th>11</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>Family_History_PD</td>
      <td>4255</td>
      <td>191512</td>
    </tr>
    <tr>
      <th>12</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>MDS_UPDRS_Part_I</td>
      <td>15971</td>
      <td>4455979</td>
    </tr>
    <tr>
      <th>13</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>MDS_UPDRS_Part_II</td>
      <td>16887</td>
      <td>4063524</td>
    </tr>
    <tr>
      <th>14</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>MDS_UPDRS_Part_III</td>
      <td>18624</td>
      <td>10595517</td>
    </tr>
    <tr>
      <th>15</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>MDS_UPDRS_Part_IV</td>
      <td>10137</td>
      <td>1330342</td>
    </tr>
    <tr>
      <th>16</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>MMSE</td>
      <td>914</td>
      <td>243770</td>
    </tr>
    <tr>
      <th>17</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>MOCA</td>
      <td>9470</td>
      <td>2774678</td>
    </tr>
    <tr>
      <th>18</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>MRI</td>
      <td>1163</td>
      <td>47510</td>
    </tr>
    <tr>
      <th>19</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>Modified_Schwab___England_ADL</td>
      <td>13127</td>
      <td>447501</td>
    </tr>
    <tr>
      <th>20</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>PDQ_39</td>
      <td>2993</td>
      <td>1373064</td>
    </tr>
    <tr>
      <th>21</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>PD_Medical_History</td>
      <td>17900</td>
      <td>1188904</td>
    </tr>
    <tr>
      <th>22</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>REM_Sleep_Behavior_Disorder_Questionnaire_Mayo</td>
      <td>2928</td>
      <td>352365</td>
    </tr>
    <tr>
      <th>23</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>REM_Sleep_Behavior_Disorder_Questionnaire_Stia...</td>
      <td>7674</td>
      <td>2288703</td>
    </tr>
    <tr>
      <th>24</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>Smoking_and_alcohol_history</td>
      <td>2796</td>
      <td>303897</td>
    </tr>
    <tr>
      <th>25</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>UPDRS</td>
      <td>668</td>
      <td>345588</td>
    </tr>
    <tr>
      <th>26</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>UPSIT</td>
      <td>4741</td>
      <td>322539</td>
    </tr>
    <tr>
      <th>27</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>amp_pd_case_control</td>
      <td>4298</td>
      <td>343991</td>
    </tr>
    <tr>
      <th>28</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>amp_pd_participants</td>
      <td>4298</td>
      <td>164495</td>
    </tr>
    <tr>
      <th>29</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>duplicate_subjects</td>
      <td>24</td>
      <td>665</td>
    </tr>
    <tr>
      <th>30</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>rna_sample_inventory</td>
      <td>8356</td>
      <td>326023</td>
    </tr>
    <tr>
      <th>31</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>wgs_gatk_joint_genotyping_samples</td>
      <td>3074</td>
      <td>74500</td>
    </tr>
    <tr>
      <th>32</th>
      <td>amp-pd-research</td>
      <td>2019_v1release_1015</td>
      <td>wgs_sample_inventory</td>
      <td>3941</td>
      <td>107446</td>
    </tr>
  </tbody>
</table>
</div>



Now we'll practice by querying the `Demographics` table. Let's say we are only interested in participant ID, baseline age, sex, and race. We can select just those columns in our query, and specify the table we are selecting from by appending the dataset name to the end of our `BQ_RELEASE_DATASET` path variable.

**TIP: Remember, if you don't want to select all available columns when querying, you can check out the column names and preview available datasets on the Google Cloud Platform**


```python
demographics = f"""

SELECT 
participant_id, age_at_baseline, sex, race

FROM `{BQ_RELEASE_DATASET}.Demographics`

"""


demo = bq_query(demographics)
```



```python
demo.info()
```

    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 4298 entries, 0 to 4297
    Data columns (total 4 columns):
    participant_id     4298 non-null object
    age_at_baseline    4298 non-null int64
    sex                4298 non-null object
    race               4297 non-null object
    dtypes: int64(1), object(3)
    memory usage: 134.4+ KB



```python
demo.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>participant_id</th>
      <th>age_at_baseline</th>
      <th>sex</th>
      <th>race</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>PP-41564</td>
      <td>22</td>
      <td>Male</td>
      <td>White</td>
    </tr>
    <tr>
      <th>1</th>
      <td>PD-PDZV843ATF</td>
      <td>24</td>
      <td>Male</td>
      <td>White</td>
    </tr>
    <tr>
      <th>2</th>
      <td>PD-PDCK871NBR</td>
      <td>24</td>
      <td>Male</td>
      <td>White</td>
    </tr>
    <tr>
      <th>3</th>
      <td>PP-51718</td>
      <td>25</td>
      <td>Male</td>
      <td>White</td>
    </tr>
    <tr>
      <th>4</th>
      <td>PP-56267</td>
      <td>25</td>
      <td>Male</td>
      <td>White</td>
    </tr>
  </tbody>
</table>
</div>



- You can manipulate datasets in the queries themselves, or simply edit them in the notebook.


```python
demo2 = demo[demo.age_at_baseline >45]
demo2.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>participant_id</th>
      <th>age_at_baseline</th>
      <th>sex</th>
      <th>race</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>138</th>
      <td>PD-PDTU629HKM</td>
      <td>46</td>
      <td>Male</td>
      <td>Black or African American</td>
    </tr>
    <tr>
      <th>139</th>
      <td>PD-PDVW790ELA</td>
      <td>46</td>
      <td>Male</td>
      <td>Black or African American</td>
    </tr>
    <tr>
      <th>140</th>
      <td>PD-PDJE363HFQ</td>
      <td>46</td>
      <td>Male</td>
      <td>Unknown</td>
    </tr>
    <tr>
      <th>141</th>
      <td>PD-PDRT993HEF</td>
      <td>46</td>
      <td>Male</td>
      <td>Asian</td>
    </tr>
    <tr>
      <th>142</th>
      <td>PD-PDWR910BPT</td>
      <td>46</td>
      <td>Male</td>
      <td>White</td>
    </tr>
  </tbody>
</table>
</div>



### **\*Important: clinical duplicates\***

In the clinical data there remain samples that have been determined to be genetically identical (either twins or simply sample duplication across studies).

The samples with lower coverage were removed from the WGS data. The two columns of the duplicated sample table below refer to the IDs of samples deleted from the WGS data (but still in the clinical data) and the IDs of samples that remain in the WGS data, respectively. 

**It is very important before running analyses using just the clinical data to investigate duplicates or you may be including the same sample twice. In addition, some of the duplicated samples from different studies have differing diagnoses.**


```python
dups = f"""

SELECT 
*

FROM `{BQ_RELEASE_DATASET}.duplicate_subjects`

"""


duplicates = bq_query(dups)
```



```python
duplicates
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>deleted_participant_id</th>
      <th>wgs_participant_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>BF-1091</td>
      <td>PD-PDEB612LPE</td>
    </tr>
    <tr>
      <th>1</th>
      <td>BF-1098</td>
      <td>PD-PDUA781LH0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>BF-1125</td>
      <td>PD-PDWY520TTB</td>
    </tr>
    <tr>
      <th>3</th>
      <td>HB-PD_INVFA733YEF</td>
      <td>BF-1130</td>
    </tr>
    <tr>
      <th>4</th>
      <td>HB-PD_INVDZ260AW1</td>
      <td>PD-PDFK551GW2</td>
    </tr>
    <tr>
      <th>5</th>
      <td>HB-PD_INVEP649ZE5</td>
      <td>PD-PDGJ885TBH</td>
    </tr>
    <tr>
      <th>6</th>
      <td>HB-PD_INVME305HYH</td>
      <td>PD-PDPV262JLC</td>
    </tr>
    <tr>
      <th>7</th>
      <td>HB-PD_INVWE905KPJ</td>
      <td>PD-PDKF434UCX</td>
    </tr>
    <tr>
      <th>8</th>
      <td>HB-PD_INVNE231ZUD</td>
      <td>PD-PDBR042ZHH</td>
    </tr>
    <tr>
      <th>9</th>
      <td>PD-PDBW494GHE</td>
      <td>BF-1100</td>
    </tr>
    <tr>
      <th>10</th>
      <td>HB-PD_INVRU618MM4</td>
      <td>PD-PDCE881ME6</td>
    </tr>
    <tr>
      <th>11</th>
      <td>HB-PD_INVXC454FPX</td>
      <td>PD-PDDW009MK9</td>
    </tr>
    <tr>
      <th>12</th>
      <td>PD-PDGU974LXR</td>
      <td>PP-3451</td>
    </tr>
    <tr>
      <th>13</th>
      <td>HB-PD_INVXA119CHN</td>
      <td>PD-PDGX770XYG</td>
    </tr>
    <tr>
      <th>14</th>
      <td>PD-PDJV686AAB</td>
      <td>BF-1088</td>
    </tr>
    <tr>
      <th>15</th>
      <td>PD-PDKV712AWZ</td>
      <td>PP-3615</td>
    </tr>
    <tr>
      <th>16</th>
      <td>PD-PDTN576BY4</td>
      <td>BF-1132</td>
    </tr>
    <tr>
      <th>17</th>
      <td>PD-PDUH505FEZ</td>
      <td>PP-3612</td>
    </tr>
    <tr>
      <th>18</th>
      <td>PD-PDVF405WVV</td>
      <td>PP-3307</td>
    </tr>
    <tr>
      <th>19</th>
      <td>PP-50463</td>
      <td>PD-PDZT881XTK</td>
    </tr>
    <tr>
      <th>20</th>
      <td>PP-50621</td>
      <td>PD-PDBD199NY6</td>
    </tr>
    <tr>
      <th>21</th>
      <td>PP-51632</td>
      <td>PD-PDAC481KNB</td>
    </tr>
    <tr>
      <th>22</th>
      <td>PP-53639</td>
      <td>PD-PDYP676BLT</td>
    </tr>
    <tr>
      <th>23</th>
      <td>PP-54916</td>
      <td>PD-PDMU756PTF</td>
    </tr>
  </tbody>
</table>
</div>



<a id="3"></a>
## 3. Querying Variant Data

Now we will practice querying variant data. You should be more careful when querying this data, as it is much larger than the clinical datasets and might get expensive if you make multiple naive queries. In most cases, you should be querying from the `gatk_passing_variants` table, as this has been reformated to make queries less expensive. More info on the `gatk_passing_variants` is in the AMP-PD Getting Started notebooks:


    The passing_variants table includes output from:

    - The Broad Institute's joint genotyping workflow
    - Variants annotated with the Variant Effect Predictor (VEP)
    - The resulting VCFs from the above process were imported to BigQuery using the Variant Transforms pipeline. 
    
    Further, the table was compacted and reshaped with:

    - Variant calls matching the reference (genotype 0, 0) moved to a column (hom_ref_call)
    - This column contains only the sample identifiers, read statistics are omitted
    - Variant calls which were no-calls (genotype -1, -1) moved to a column (no_call)
    - This column contains only the sample identifiers read statistics are omitted
    - Variants not matching FILTER = PASS have been omitted
    
    

- Best practice is to target your variant queries. This means looking at the specific regions or variants you're interested in rather than querying the whole datatset.

- The query below looks specifically at LRRK2 variants classified as missence variants.


```python
VARIANT_QUERY = f"""
    SELECT
      reference_name, start_position, end_position, reference_bases, names,
      ARRAY (
        SELECT AS STRUCT alt,
          ARRAY (
            SELECT DISTINCT AS STRUCT
               allele, ALLELE_NUM, VARIANT_CLASS, Consequence, IMPACT, SYMBOL, Gene,
               Protein_position, Amino_acids, Codons, Existing_variation, CLIN_SIG
            FROM UNNEST(CSQ)
            WHERE REGEXP_CONTAINS(Consequence, r'\\bmissense_variant\\b')
          ) AS CSQ
        FROM UNNEST(v.alternate_bases)
       ) AS alternate_bases,
      ARRAY (
        SELECT AS STRUCT name, genotype, AD FROM UNNEST(v.call)
      ) AS call
    FROM
      `amp-pd-research.2019_v1release_1015_genomics.gatk_passing_variants` v
    WHERE
      reference_name = 'chr12' AND
      EXISTS (SELECT 1 FROM UNNEST(v.alternate_bases) a
              WHERE
                EXISTS (SELECT 1 from a.CSQ
                        WHERE SYMBOL = 'LRRK2' AND
                              REGEXP_CONTAINS(Consequence, r'\\bmissense_variant\\b')))
"""
variant1 = bq_query(VARIANT_QUERY)
```



```python
variant1.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>reference_name</th>
      <th>start_position</th>
      <th>end_position</th>
      <th>reference_bases</th>
      <th>names</th>
      <th>alternate_bases</th>
      <th>call</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr12</td>
      <td>40251272</td>
      <td>40251273</td>
      <td>G</td>
      <td>[rs78501232]</td>
      <td>[{'alt': 'A', 'CSQ': [{'allele': 'A', 'ALLELE_...</td>
      <td>[{'name': 'BF-1253', 'genotype': [0, 1], 'AD':...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr12</td>
      <td>40309108</td>
      <td>40309109</td>
      <td>G</td>
      <td>[rs7133914]</td>
      <td>[{'alt': 'A', 'CSQ': [{'allele': 'A', 'ALLELE_...</td>
      <td>[{'name': 'BF-1009', 'genotype': [0, 1], 'AD':...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr12</td>
      <td>40320098</td>
      <td>40320099</td>
      <td>T</td>
      <td>[rs11564148]</td>
      <td>[{'alt': 'A', 'CSQ': [{'allele': 'A', 'ALLELE_...</td>
      <td>[{'name': 'BF-1001', 'genotype': [0, 1], 'AD':...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr12</td>
      <td>40313975</td>
      <td>40313976</td>
      <td>G</td>
      <td>[rs35507033]</td>
      <td>[{'alt': 'A', 'CSQ': [{'allele': 'A', 'ALLELE_...</td>
      <td>[{'name': 'BF-1038', 'genotype': [0, 1], 'AD':...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>chr12</td>
      <td>40225279</td>
      <td>40225280</td>
      <td>G</td>
      <td>[rs2256408]</td>
      <td>[{'alt': 'A', 'CSQ': [{'allele': 'A', 'ALLELE_...</td>
      <td>[{'name': 'BF-1001', 'genotype': [1, 1], 'AD':...</td>
    </tr>
  </tbody>
</table>
</div>



Many of the columns in the BigQuery variant tables are nested to save space. You can get a better look at the information in these columns by using something like the function below that prints out the contents in a more readable format.


```python
def print_variant(row):
    print(f"Variant reference information:")
    print(f"  Position: {row.reference_name}:{row.start_position}-{row.end_position}")
    print(f"  Reference base: {row.reference_bases}")
    
    alternate_names = row.names
    if alternate_names:
        print(f"  Names: {sorted(set(alternate_names))}")

    print()
    print(f"  Number of alternates (variants): {len(row.alternate_bases)}")
    for alternate in row.alternate_bases:
        print(f"    Alternate base(s): {alternate['alt']}")

        print(f"    Number of consequences(annotations): {len(alternate['CSQ'])}")
        for csq in alternate['CSQ']:
            print(f"    Consequence: {csq['Consequence']}")
            print(f"      Variant class: {csq['VARIANT_CLASS']}")
            print(f"      Ensembl gene: {csq['Gene']}")
            print(f"      AA position: {csq['Protein_position']}")
            print(f"      AA change: {csq['Amino_acids']}")
            print(f"      Codon change: {csq['Codons']}")
            print(f"      Known identifier: {csq['Existing_variation']}")
            if csq['CLIN_SIG']:
                print(f"      ClinVar clinical significance: {csq['CLIN_SIG'].split('&')}")

    print()
    print(f"  Number of samples with this call: {len(row.call)}")
    for call in row.call:
      print(f"    {call['name']}: genotype: {call['genotype']}, depth: {call['AD']}")

# Dump out all of the variant rows
for row in variant1.itertuples():
    print_variant(row)
    print()
```

- You can use the code in the cell below to get a list of samples IDs with calls for the variants you are interested in for future analyses.


```python
sample_ids = []
for row in variant1.itertuples():
    calls = row.call
    for call in calls:
        sample_ids.append(call['name'])

sample_ids = pd.DataFrame(sorted(set(sample_ids)))
```


```python
sample_ids.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>BF-1001</td>
    </tr>
    <tr>
      <th>1</th>
      <td>BF-1002</td>
    </tr>
    <tr>
      <th>2</th>
      <td>BF-1003</td>
    </tr>
    <tr>
      <th>3</th>
      <td>BF-1004</td>
    </tr>
    <tr>
      <th>4</th>
      <td>BF-1005</td>
    </tr>
  </tbody>
</table>
</div>



### Writing results to your Workspace Bucket


```python
with open('test_sample_ids.csv', 'w') as f:
    f.write(sample_ids.to_csv(index=False))
```


```bash
%%bash

ls
```

    BigQuerySampleQC.ipynb
    BigQueryVariantQC.ipynb
    compareNeuroXvsSeq.ipynb
    Example - Running a GWAS-Copy1.ipynb
    Example - Running a GWAS - deflaux.ipynb
    flagByGATKsummaryMetrics.ipynb
    [0m[01;34mgs:[0m/
    [01;34mhampton_temp[0m/
    Hampton_testing.ipynb
    HI0_getStarted.ipynb
    HI1_flagByGATKsummaryMetrics.ipynb
    HI2_compareNeuroXvsSeq.ipynb
    HI3_plinkQC.ipynb
    HI3_qcAmpPDplink.ipynb
    HI4_reportQCresults.ipynb
    Hirotaka_ClinicalData.ipynb
    HIX_HeterozygosityEvaluation.ipynb
    How_to_use_Terra_Training_Course.ipynb
    [01;34mhtslib[0m/
    Matt's copy of Hampton's notebook.ipynb
    Pulling LRRK2 GBA SNCA Samples.ipynb
    Py3 - Start Here (Hirotaka).ipynb
    Py3 - WGS - Start Here(copy).ipynb
    qcAmpPDplink2.ipynb
    QC Notebook (Hirotaka).ipynb
    QC Notebook (Matt) - Python3.ipynb
    Rbigquery.ipynb
    R Environment Setup .ipynb
    reportQCresults.ipynb
    R - Save and load R objects from workspace bucket.ipynb
    Taking a closer look at the star alternate allele.ipynb
    Terra_for_AMP_PD_Video1.ipynb
    testGWAS.ipynb
    Test - Running a GWAS.ipynb
    test_sample_ids.csv
    VCF_Sample_Filter.ipynb
    WGS_sample_IDs.csv



```python
shell_do(f'gsutil -u {BILLING_PROJECT_ID} cp -r ./test_sample_ids.csv {WORKSPACE_BUCKET}')
```

    Executing: gsutil -u fc-amp-pd-alpha cp -r ./test_sample_ids.csv gs://fc-secure-e4c969f8-ebed-4ab2-bd4c-a53cd897cddb
    Copying file://./test_sample_ids.csv [Content-Type=text/csv]...
    / [1 files][ 33.4 KiB/ 33.4 KiB]                                                
    Operation completed over 1 objects/33.4 KiB.                                     


<a id="4"></a>
## 4. Pulling Data from Google Storage Buckets

- Querying is a quick and easy way to pull the data you're interested in, but you can also read certain datasets directly from the Google Cloud Storage Buckets.
- You can submit the command below to the shell using our path variables to look at the Google Bucket locations for the clinical data.


```python
shell_do(f'gsutil -u {BILLING_PROJECT_ID} ls {GS_CLINICAL_RELEASE_PATH}')
```

    Executing: gsutil -u fc-amp-pd-alpha ls gs://amp-pd-data/releases/2019_v1release_1015/clinical
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Biospecimen_analyses_CSF_abeta_tau_ptau.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Biospecimen_analyses_CSF_abeta_tau_ptau_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Biospecimen_analyses_CSF_beta_glucocerebrosidase.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Biospecimen_analyses_CSF_beta_glucocerebrosidase_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Biospecimen_analyses_SomaLogic_plasma.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Biospecimen_analyses_SomaLogic_plasma_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Biospecimen_analyses_other.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Biospecimen_analyses_other_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Caffeine_history.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Caffeine_history_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/DTI.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/DTI_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/DaTSCAN_SBR.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/DaTSCAN_SBR_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/DaTSCAN_visual_interpretation.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/DaTSCAN_visual_interpretation_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Demographics.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Demographics_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Enrollment.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Enrollment_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Epworth_Sleepiness_Scale.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Epworth_Sleepiness_Scale_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Family_History_PD.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Family_History_PD_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MDS_UPDRS_Part_I.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MDS_UPDRS_Part_II.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MDS_UPDRS_Part_III.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MDS_UPDRS_Part_III_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MDS_UPDRS_Part_II_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MDS_UPDRS_Part_IV.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MDS_UPDRS_Part_IV_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MDS_UPDRS_Part_I_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MMSE.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MMSE_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MOCA.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MOCA_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MRI.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/MRI_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Modified_Schwab___England_ADL.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Modified_Schwab___England_ADL_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/PDQ_39.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/PDQ_39_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/PD_Medical_History.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/PD_Medical_History_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/REM_Sleep_Behavior_Disorder_Questionnaire_Mayo.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/REM_Sleep_Behavior_Disorder_Questionnaire_Mayo_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/REM_Sleep_Behavior_Disorder_Questionnaire_Stiasny_Kolster.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/REM_Sleep_Behavior_Disorder_Questionnaire_Stiasny_Kolster_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Smoking_and_alcohol_history.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/Smoking_and_alcohol_history_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/UPDRS.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/UPDRS_dictionary.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/UPSIT.csv
    gs://amp-pd-data/releases/2019_v1release_1015/clinical/UPSIT_dictionary.csv


- Let's read in the `Demographics` CSV


```python
demographics_df = gcs_read_csv(os.path.join(GS_CLINICAL_RELEASE_PATH, 'Demographics.csv'))

demographics_df.info()
```

    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 4298 entries, 0 to 4297
    Data columns (total 9 columns):
    participant_id           4298 non-null object
    GUID                     2688 non-null object
    visit_name               4298 non-null object
    visit_month              4298 non-null int64
    age_at_baseline          4298 non-null int64
    sex                      4298 non-null object
    ethnicity                4297 non-null object
    race                     4297 non-null object
    education_level_years    4297 non-null object
    dtypes: int64(2), object(7)
    memory usage: 302.3+ KB



```python
demographics_df.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>participant_id</th>
      <th>GUID</th>
      <th>visit_name</th>
      <th>visit_month</th>
      <th>age_at_baseline</th>
      <th>sex</th>
      <th>ethnicity</th>
      <th>race</th>
      <th>education_level_years</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>PD-PDNE299YPT</td>
      <td>PDNE299YPT</td>
      <td>M0</td>
      <td>0</td>
      <td>74</td>
      <td>Male</td>
      <td>Not Hispanic or Latino</td>
      <td>White</td>
      <td>Greater than 16 years</td>
    </tr>
    <tr>
      <th>1</th>
      <td>PD-PDYW828VAV</td>
      <td>PDYW828VAV</td>
      <td>M0</td>
      <td>0</td>
      <td>57</td>
      <td>Female</td>
      <td>Not Hispanic or Latino</td>
      <td>White</td>
      <td>12-16 years</td>
    </tr>
    <tr>
      <th>2</th>
      <td>PD-PDWV958JPG</td>
      <td>PDWV958JPG</td>
      <td>M0</td>
      <td>0</td>
      <td>68</td>
      <td>Male</td>
      <td>Not Hispanic or Latino</td>
      <td>White</td>
      <td>Greater than 16 years</td>
    </tr>
    <tr>
      <th>3</th>
      <td>PD-PDRF387ZV3</td>
      <td>PDRF387ZV3</td>
      <td>M0</td>
      <td>0</td>
      <td>68</td>
      <td>Male</td>
      <td>Not Hispanic or Latino</td>
      <td>White</td>
      <td>12-16 years</td>
    </tr>
    <tr>
      <th>4</th>
      <td>PD-PDAC268KWV</td>
      <td>PDAC268KWV</td>
      <td>M0</td>
      <td>0</td>
      <td>69</td>
      <td>Female</td>
      <td>Not Hispanic or Latino</td>
      <td>White</td>
      <td>Greater than 16 years</td>
    </tr>
  </tbody>
</table>
</div>



- Some data formats aren't queryable, such as CRAMS, VCF files, and FASTQs
- To use these data formats, we just need to copy the files to our cluster

First let's look at the locations in the cloud for these files:


```python
shell_do(f'gsutil -u {BILLING_PROJECT_ID} ls {GS_WGS_RELEASE_PATH}')
```

    Executing: gsutil -u fc-amp-pd-alpha ls gs://amp-pd-genomics/releases/2019_v1release_1015/wgs
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/wgs_samples.csv
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/


Let's take a look in the `gatk` directory


```python
shell_do(f'gsutil -u {BILLING_PROJECT_ID} ls {GS_WGS_RELEASE_GATK_PATH}')
```

    Executing: gsutil -u fc-amp-pd-alpha ls gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/metrics/
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/


- In the `gatk` directory, we have metrics files and per-chromosome VCFs.
- We also have per sample VCFs, CRAMS, and metrics in the `samples` directory: `gs://amp-pd-genomics/samples/wgs/gatk`.


```python
shell_do(f'gsutil -u {BILLING_PROJECT_ID} ls {GS_WGS_RELEASE_GATK_PATH}/vcf')
```

    Executing: gsutil -u fc-amp-pd-alpha ls gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr1.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr1.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr1.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr10.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr10.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr10.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr11.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr11.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr11.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr12.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr12.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr12.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr13.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr13.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr13.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr14.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr14.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr14.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr15.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr15.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr15.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr16.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr16.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr16.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr17.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr17.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr17.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr18.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr18.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr18.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr19.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr19.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr19.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr2.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr2.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr2.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr20.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr20.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr20.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr21.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr21.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr21.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr22.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr22.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr22.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr3.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr3.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr3.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr4.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr4.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr4.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr5.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr5.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr5.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr6.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr6.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr6.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr7.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr7.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr7.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr8.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr8.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr8.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr9.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr9.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chr9.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chrX.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chrX.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chrX.vcf.idx
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chrY.vcf.gz
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chrY.vcf.gz.tbi
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/gatk/vcf/chrY.vcf.idx



```python
shell_do(f'gsutil -u {BILLING_PROJECT_ID} ls {GS_WGS_RELEASE_PLINK_PATH}/bfile/all_vcfs*')
```

    Executing: gsutil -u fc-amp-pd-alpha ls gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/bfile
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/bfile/README.txt
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/bfile/all_vcfs.bed
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/bfile/all_vcfs.bim
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/bfile/all_vcfs.fam
    gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/bfile/all_vcfs.log


Now we can practice copying data from the cloud storage to our cluster. We'll need the Plink binary files for our analysis in a bit, so let's copy those.

- First let's make a directory to copy it to.


```bash
%%bash

mkdir ~/bin/data_temp
cd ~/bin/
ls
```

    bcftools
    data
    data_temp
    htslib
    king
    Linux-king
    plink
    plink_linux_x86_64_20190304


    mkdir: cannot create directory ‘/home/jupyter-user/bin/data_temp’: File exists



```python
shell_do(f'gsutil -mu {BILLING_PROJECT_ID} cp {GS_WGS_RELEASE_PLINK_PATH}/bfile/all_vcfs* ~/bin/data_temp')
```

    Executing: gsutil -mu fc-amp-pd-alpha cp gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/bfile/all_vcfs* ~/bin/data_temp
    Copying gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/bfile/all_vcfs.bim...
    Copying gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/bfile/all_vcfs.log...
    Copying gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/bfile/all_vcfs.bed...
    ==> NOTE: You are downloading one or more large file(s), which would            
    run significantly faster if you enabled sliced object downloads. This
    feature is enabled by default but requires that compiled crcmod be
    installed (see "gsutil help crcmod").
    
    Copying gs://amp-pd-genomics/releases/2019_v1release_1015/wgs/plink/bfile/all_vcfs.fam...
    | [4/4 files][ 21.4 GiB/ 21.4 GiB] 100% Done  54.7 MiB/s ETA 00:00:00           
    Operation completed over 4 objects/21.4 GiB.                                     



```bash
%%bash

cd ~/bin/data_temp
ls
```

    all_vcfs.bed
    all_vcfs.bim
    all_vcfs.fam
    all_vcfs.log



```bash
%%bash
cd ~/bin/data_temp

head all_vcfs.bim

```

    1	rs201752861	0	10177	*	A
    1	rs201694901	0	10180	C	T
    1	rs199706086	0	10250	C	A
    1	rs140194106	0	10254	*	TA
    1	rs145427775	0	10291	T	C
    1	rs112750067	0	10327	C	T
    1	rs55998931	0	10492	T	C
    1	rs62636508	0	10519	C	G
    1	rs58108140	0	10583	A	G
    1	rs189107123	0	10611	G	C


<a id="5"></a>
## 5. Quick Association Test with Plink

Now that we have the basics down, we'll run through a brief example of how to do an analysis on a Terra notebook.

Using the Plink binary files that we copied to our cluster previously, we'll do a Plink case/control association test. For this, we will be using the `fisher` flag, which performs Fisher's exact test. Fisher's test is used to determine if there are nonrandom associations between two variables. Our variables here are the presence of a Parkinon's Disease diagnosis (1 if yes, 0 if no), and the major and minor allele distribution for the tested variant. So for every variant in our data, we are testing whether the distribution of alleles is significantly different between the PD group and the non PD group. 


To do this you'll need a few extra things:

- Plink installed on your cluster
- Phenotype file
- Sex info file 

### Installing your own packages

- You can get package files ready for installation in different ways.
    - Get the files from online, upload them to the workspace bucket, copy them to the cluster.
    - Download straight from a code cell using commands like `git clone`.
- You will then follow the regular instructions for installation.

We are going to use `wget` to grab Plink 1.9 for Linux from http://s3.amazonaws.com


```bash
%%bash

if test -e ~/bin/plink; then

echo "Plink is already installed in ~/bin"
else
echo "Plink is not installed"
cd ~/bin/

wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190304.zip 

fi
```

    Plink is already installed in ~/bin



```bash
%%bash

ls ~/bin/

```

    bcftools
    data
    data_temp
    htslib
    king
    LICENSE
    Linux-king
    plink
    plink_linux_x86_64_20190304.zip
    prettify
    toy.map
    toy.ped



```python
!~/bin/plink 
```

    PLINK v1.90b6.9 64-bit (4 Mar 2019)            www.cog-genomics.org/plink/1.9/
    (C) 2005-2019 Shaun Purcell, Christopher Chang   GNU General Public License v3
    
      plink <input flag(s)...> [command flag(s)...] [other flag(s)...]
      plink --help [flag name(s)...]
    
    Commands include --make-bed, --recode, --flip-scan, --merge-list,
    --write-snplist, --list-duplicate-vars, --freqx, --missing, --test-mishap,
    --hardy, --mendel, --ibc, --impute-sex, --indep-pairphase, --r2, --show-tags,
    --blocks, --distance, --genome, --homozyg, --make-rel, --make-grm-gz,
    --rel-cutoff, --cluster, --pca, --neighbour, --ibs-test, --regress-distance,
    --model, --bd, --gxe, --logistic, --dosage, --lasso, --test-missing,
    --make-perm-pheno, --tdt, --qfam, --annotate, --clump, --gene-report,
    --meta-analysis, --epistasis, --fast-epistasis, and --score.
    
    'plink --help | more' describes all functions (warning: long).


Now we need our phenotype and covariate files. I'm going to just query a few data fields like we did before. I'm just going to grab diagnosis, sex, and age for this quick example. In reality, you would want to restrict analyses per ancestry and also want to include more covariates like PCs generated by Plink.

**Tip: Get to know the cohorts before running analyses, cohorts like PPMI have study arms that are enriched for certain variants and you may need to exclude those samples** 


We need to grab sex information to update the Plink fam file, as it currently does not have sex coded into it. Chi squared tests like Plink's association tests don't use covariates, but Plink still needs sex info or it sets those samples to missing. We'll grab the sex info from the demographics file and use it to update the Plink fam file.


```python
covariates = f"""

SELECT 
participant_id, sex 

FROM `{BQ_RELEASE_DATASET}.Demographics`

"""


covs = bq_query(covariates)
```



```python
covs.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>participant_id</th>
      <th>sex</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>PP-41564</td>
      <td>Male</td>
    </tr>
    <tr>
      <th>1</th>
      <td>PD-PDZV843ATF</td>
      <td>Male</td>
    </tr>
    <tr>
      <th>2</th>
      <td>PD-PDCK871NBR</td>
      <td>Male</td>
    </tr>
    <tr>
      <th>3</th>
      <td>PP-51718</td>
      <td>Male</td>
    </tr>
    <tr>
      <th>4</th>
      <td>PP-56267</td>
      <td>Male</td>
    </tr>
  </tbody>
</table>
</div>



There is a case control table available that includes baseline and most recent diagnosis. Healthy samples are controls, any Parkinson's disease or idiopathic Parkinson's diagnsoses are classified as cases, and any other related diagnoses are classified as other.


```python
phenotype = f"""

SELECT 
* 

FROM `{BQ_RELEASE_DATASET}.amp_pd_case_control`

"""

pheno = bq_query(phenotype)
```




```python
pheno.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>participant_id</th>
      <th>diagnosis_at_baseline</th>
      <th>diagnosis_latest</th>
      <th>case_control_other_at_baseline</th>
      <th>case_control_other_latest</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>PD-PDAB597LKL</td>
      <td>Parkinsonism</td>
      <td>Parkinsonism</td>
      <td>Other</td>
      <td>Other</td>
    </tr>
    <tr>
      <th>1</th>
      <td>PD-PDCT406CN9</td>
      <td>Parkinsonism</td>
      <td>Parkinsonism</td>
      <td>Other</td>
      <td>Other</td>
    </tr>
    <tr>
      <th>2</th>
      <td>PP-54461</td>
      <td>Prodromal motor PD</td>
      <td>Idiopathic PD</td>
      <td>Other</td>
      <td>Case</td>
    </tr>
    <tr>
      <th>3</th>
      <td>PP-54215</td>
      <td>Prodromal motor PD</td>
      <td>Idiopathic PD</td>
      <td>Other</td>
      <td>Case</td>
    </tr>
    <tr>
      <th>4</th>
      <td>PP-18567</td>
      <td>Prodromal motor PD</td>
      <td>Idiopathic PD</td>
      <td>Other</td>
      <td>Case</td>
    </tr>
  </tbody>
</table>
</div>



To use Plink, we need to get our sex info file and phenotype file in Plink format, which means 2 columns of IDs followed by the other information.


```python
covsp = covs
covsp['IID'] = covsp['participant_id']
covsp = covsp[['participant_id', 'IID', 'sex']]
covsp['sex'] = covsp['sex'].astype(str)
covsp.sex[(covsp.sex == "Male")] = 1
covsp.sex[(covsp.sex == "Female")] = 2
covsp.columns = ['FID','IID', 'sex']
```


```python
covsp.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>FID</th>
      <th>IID</th>
      <th>sex</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>PP-41564</td>
      <td>PP-41564</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>PD-PDZV843ATF</td>
      <td>PD-PDZV843ATF</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>PD-PDCK871NBR</td>
      <td>PD-PDCK871NBR</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>PP-51718</td>
      <td>PP-51718</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>PP-56267</td>
      <td>PP-56267</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




```python
phenop = pheno
phenop['IID'] = phenop['participant_id']
phenop = phenop[['participant_id', 'IID', 'case_control_other_latest']]
phenop.columns = ['FID','IID', 'PHENO']
phenop = phenop[(phenop.PHENO != 'Other')]
phenop.PHENO[(phenop.PHENO == "Case")] = 2
phenop.PHENO[(phenop.PHENO == "Control")] = 1
```


```python
phenop.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>FID</th>
      <th>IID</th>
      <th>PHENO</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2</th>
      <td>PP-54461</td>
      <td>PP-54461</td>
      <td>2</td>
    </tr>
    <tr>
      <th>3</th>
      <td>PP-54215</td>
      <td>PP-54215</td>
      <td>2</td>
    </tr>
    <tr>
      <th>4</th>
      <td>PP-18567</td>
      <td>PP-18567</td>
      <td>2</td>
    </tr>
    <tr>
      <th>5</th>
      <td>PP-40615</td>
      <td>PP-40615</td>
      <td>2</td>
    </tr>
    <tr>
      <th>6</th>
      <td>PP-40612</td>
      <td>PP-40612</td>
      <td>2</td>
    </tr>
  </tbody>
</table>
</div>



Let's write these to tab files so we can use them with Plink


```python
with open('plink_test_covs.tab', 'w') as f:
    f.write(covsp.to_csv(index=False, sep='\t'))
```


```python
with open('plink_test_pheno.tab', 'w') as f:
    f.write(phenop.to_csv(index=False, sep='\t'))
```


```bash
%%bash 

ls
pwd
```

Now we'll use Plink's --update-sex command to update the fam file with the information we queried from the demographics file.


```python
shell_do(f'~/bin/plink --bfile ~/bin/data_temp/all_vcfs \
            --update-sex /home/jupyter-user/notebooks/AMP*PD*WGS*QC*Collaboration/edit/plink_test_covs.tab \
            --make-bed \
            --out ~/bin/data_temp/plink_updated_sex \
            --threads 2')
```

    Executing: ~/bin/plink --bfile ~/bin/data_temp/all_vcfs             --update-sex /home/jupyter-user/notebooks/AMP*PD*WGS*QC*Collaboration/edit/plink_test_covs.tab             --make-bed             --out ~/bin/data_temp/plink_updated_sex             --threads 2
    PLINK v1.90b6.9 64-bit (4 Mar 2019)            www.cog-genomics.org/plink/1.9/
    (C) 2005-2019 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to /home/jupyter-user/bin/data_temp/plink_updated_sex.log.
    Options in effect:
      --bfile /home/jupyter-user/bin/data_temp/all_vcfs
      --make-bed
      --out /home/jupyter-user/bin/data_temp/plink_updated_sex
      --threads 2
      --update-sex /home/jupyter-user/notebooks/AMP PD WGS QC Collaboration/edit/plink_test_covs.tab
    
    15043 MB RAM detected; reserving 7521 MB for main workspace.
    28791969 variants loaded from .bim file.
    3074 people (0 males, 0 females, 3074 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to
    /home/jupyter-user/bin/data_temp/plink_updated_sex.nosex .
    --update-sex: 3074 people updated, 1225 IDs not present.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 3074 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Warning: 14673747 het. haploid genotypes present (see
    /home/jupyter-user/bin/data_temp/plink_updated_sex.hh ); many commands treat
    these as missing.
    Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    treat these as missing.
    Total genotyping rate is 0.993972.
    28791969 variants and 3074 people pass filters and QC.
    Note: No phenotypes present.
    --make-bed to /home/jupyter-user/bin/data_temp/plink_updated_sex.bed +
    /home/jupyter-user/bin/data_temp/plink_updated_sex.bim +
    /home/jupyter-user/bin/data_temp/plink_updated_sex.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.


### Running a Plink command

Now that our files are done we can perform our analysis. We will use the Plink `--assoc` command to run an association test and the `fisher` flag to perform Fisher's exact test instead of the default chi square test. We do this because the chi squared test is an approximate test and Fisher's test is exact, which is more accurate in cases where expected frequencies are low as can happen with rarer alleles.

We also will use the `--maf` flag to restrict our analysis to variants with a minor allele frequency greater that 5%.

The `--adjust` flag outputs an additional file that includes P values adjusted for multiple testing with a variety of methods.

**Tip: After running the command below, you will see the genomic inflation factor is higher than we want. Samples are not stratified by population in this example. In a real analysis you would need to run separate analyses per population or employ some of Plink's population clustering techniques and specify the analysis to run within clusters.**

**In addition, best practice would be to perform variant QC beforehand, such as removing missing by case/control and filtering controls by Hardy-Weinberg equilibrium. No variant QC has been done to these plink files other than removing variants that failed the GATK QC pipeline, so that specific choices could be left up to individual researchers. More information about to perform proper QC is provided in future modules of this course.**


```python
shell_do(f'~/bin/plink --bfile ~/bin/data_temp/plink_updated_sex \
            --pheno /home/jupyter-user/notebooks/AMP*PD*WGS*QC*Collaboration/edit/plink_test_pheno.tab \
            --pheno-name PHENO \
            --assoc fisher \
            --adjust \
            --maf 0.05 \
            --geno 0.05 \
            --out ~/bin/data_temp/test_PD \
            --set-hh-missing \
            --threads 2')
```

    Executing: ~/bin/plink --bfile ~/bin/data_temp/plink_updated_sex             --pheno /home/jupyter-user/notebooks/AMP*PD*WGS*QC*Collaboration/edit/plink_test_pheno.tab             --pheno-name PHENO             --assoc fisher             --adjust             --maf 0.05             --geno 0.05             --out ~/bin/data_temp/test_PD             --set-hh-missing             --threads 2
    PLINK v1.90b6.9 64-bit (4 Mar 2019)            www.cog-genomics.org/plink/1.9/
    (C) 2005-2019 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to /home/jupyter-user/bin/data_temp/test_PD.log.
    Options in effect:
      --adjust
      --assoc fisher
      --bfile /home/jupyter-user/bin/data_temp/plink_updated_sex
      --geno 0.05
      --maf 0.05
      --out /home/jupyter-user/bin/data_temp/test_PD
      --pheno /home/jupyter-user/notebooks/AMP PD WGS QC Collaboration/edit/plink_test_pheno.tab
      --pheno-name PHENO
      --set-hh-missing
      --threads 2
    
    15043 MB RAM detected; reserving 7521 MB for main workspace.
    Allocated 5640 MB successfully, after larger attempt(s) failed.
    28791969 variants loaded from .bim file.
    3074 people (1721 males, 1353 females) loaded from .fam.
    2835 phenotype values present after --pheno.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 3074 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Warning: 14673747 het. haploid genotypes present (see
    /home/jupyter-user/bin/data_temp/test_PD.hh ); many commands treat these as
    missing.
    Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    treat these as missing.
    Total genotyping rate is 0.993972.
    571564 variants removed due to missing genotype data (--geno).
    21146776 variants removed due to minor allele threshold(s)
    (--maf/--max-maf/--mac/--max-mac).
    7073629 variants and 3074 people pass filters and QC.
    Among remaining phenotypes, 1735 are cases and 1100 are controls.  (239
    phenotypes are missing.)
    Writing C/C --assoc report to
    /home/jupyter-user/bin/data_temp/test_PD.assoc.fisher ...
    10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788898990919293949596979899done.
    --adjust: Genomic inflation est. lambda (based on median chisq) = 1.31553.
    101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899--adjust values (7073399 variants) written to
    /home/jupyter-user/bin/data_temp/test_PD.assoc.fisher.adjusted .


Now we can look at our output file "test_PD.assoc.fisher" The fields are:

- Chromosome
- SNP identifier
- Code for allele 1 (the minor, rare allele based on the entire sample frequencies)
- The frequency of this variant in cases
- The frequency of this variant in controls
- Code for the other allele
- The asymptotic significance value for this test
- The odds ratio for this test


```bash
%%bash

head ~/bin/data_temp/test_PD.assoc.fisher

```

     CHR           SNP         BP   A1      F_A      F_U   A2            P           OR 
       1   rs200683566      13838    T   0.3796   0.3533    C      0.04984         1.12 
       1   rs370886505      14397    C  0.06297  0.05377 CTGT       0.1646        1.183 
       1   rs375086259      14653    T   0.3796   0.3863    C       0.6122       0.9718 
       1   rs201327123      14677    A   0.1068   0.1015    G       0.5333        1.059 
       1    rs79585140      14907    G   0.4017   0.4376    A     0.007857        0.863 
       1    rs75454623      14930    G    0.382   0.4262    A     0.000999       0.8325 
       1   rs199856693      14933    A  0.06102  0.07279    G      0.08758       0.8277 
       1    rs71252250      15118    G   0.3172   0.3339    A       0.1896       0.9265 
       1    rs78601809      15211    G   0.1518   0.2264    T    6.876e-12       0.6115 


Let's look at the adjusted file as well to get our adjusted P values. The fields are:

- Chromosome
- SNP identifier
- Unadjusted, asymptotic significance value
- Genomic control adjusted significance value. This is based on a simple estimation of the inflation factor based on   median chi-square statistic. These values do not control for multiple testing therefore.
- Bonferroni adjusted significance value
- Holm step-down adjusted significance value
- Sidak single-step adjusted significance value
- Sidak step-down adjusted significance value
- Benjamini & Hochberg (1995) step-up FDR control
- Benjamini & Yekutieli (2001) step-up FDR control


```bash
%%bash

head ~/bin/data_temp/test_PD.assoc.fisher.adjusted

```

     CHR           SNP      UNADJ         GC       BONF       HOLM   SIDAK_SS   SIDAK_SD     FDR_BH     FDR_BY
      12     rs7138781 8.713e-110  6.464e-84 6.163e-103 6.163e-103        INF        INF 6.163e-103 1.008e-101 
       6   rs369380078  7.436e-94  8.489e-72   5.26e-87   5.26e-87        INF        INF   2.63e-87    4.3e-86 
       6   rs372929245  1.121e-89   1.28e-68  7.928e-83  7.928e-83        INF        INF  2.643e-83  4.321e-82 
       6   rs376842100  4.227e-86  6.731e-66   2.99e-79   2.99e-79        INF        INF  7.474e-80  1.222e-78 
      16   rs181475825  1.158e-83  4.816e-64  8.192e-77  8.192e-77        INF        INF  1.638e-77  2.679e-76 
       6   rs368785910  7.327e-69  8.791e-53  5.183e-62  5.183e-62        INF        INF  8.638e-63  1.412e-61 
       6    rs28573898  6.131e-65   8.49e-50  4.337e-58  4.337e-58        INF        INF  6.196e-59  1.013e-57 
       5    rs10942417  1.224e-57  3.051e-44  8.658e-51  8.658e-51        INF        INF  1.082e-51  1.769e-50 
      20   rs371979500  2.579e-55  1.791e-42  1.825e-48  1.825e-48        INF        INF  2.027e-49  3.314e-48 


### Provenance


```python
import datetime
print(datetime.datetime.now())
```

    2020-05-18 15:51:25.611655

