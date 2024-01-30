# MIA-Akoya single cell knitting across different data modalities
Main contributor: Aditya Pratapa (spatial), Cam Quek and Ghamdan Al-Eryani (single-cell sequencing).

In this study, we constructed a Multimodal Integration Toolkit (MIT), an improved data analysis pipeline, to align each cellular phenotypes with their individual transcriptomic features, cellular epitopes and spatial information at a single cell level for deeper insights into the tumour microenvironment of melanoma in response to immunotherapy.

In order to run MIT, use the following code:
```
docker run -d --name sigcm -p 8787:8787  -v "\<source folder>\:/home/rstudio/singlecellmel" -e USER=rstudio -e PASSWORD=singlecellmel123 adyprat/singlecellmel
```

To view integrated .h5ad object, use cellxgene (https://github.com/interactivereport/cellxgene_VIP) to launch:
```
--host 0.0.0.0 --port 54890 inputs/P1_Transferred_092022_1.h5ad
```

Please note that the size of the remaining raw imaging and read count data are too large to be stored in a public repository and will therefore be stored in a private cloud-based server. Access to these data sets will be provided upon request.

