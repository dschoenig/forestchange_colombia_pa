[![DOI](https://zenodo.org/badge/287354868.svg)](https://zenodo.org/badge/latestdoi/287354868)


# The effect of post-conflict transition on deforestation of protected areas in Colombia

This repository contains analysis scripts and results that arose from
replicating and extending upon a recent study on deforestation in Colombian
protected areas. Our results strongly put in question the conclusions presented
in the original study: Reporting an increase in deforestation after ending armed
conflict, the authors of this study propose several drivers behind this trend
that are supposed to specifically affect protected areas and render them
particularly vulnerable to deforestation during post-conflict transition. In our
reanalysis, it has become apparent that the original study merely picked up a
national trend of increased deforestation and that forests in national protected
areas are actually much less affected by the transition than other forests in
Colombia -- it may even be argued that protected areas have become more
effective at slowing down deforestation when compared to non-protected forests.
The drivers and conservation lessons proposed in original study are therefore
highly speculative. We are deeply concerned by the general increase in
deforestation in Colombia. And we believe it is important to analyze potential
drivers of deforestation more comprehensively, so that adequate measures to
reduce forest loss can be identified.

## Code overview

* R scripts for geospatial and statistical analyses are placed in the `src`
  folder.
* Statistical analyses are presented in `src/3_analyses.R` and can be run using
  the results (i.e. `.csv` files) that are generated by the geospatial analysis,
  and which are included in this repository.
* Data preparation and geospatial analysis (in this order) are performed using
  the scripts `src/1_data_preperation.R`, and `src/2_forest_loss.R`,
  respecitvely. To run the scripts, it is necessary to download the raw data
  first and place it in the `data/` directory as described below.

## Sources for raw data

### Forest loss
We used version `1.6` of the global forest change data sets provided by [Hansen
et al. (2013)](https://doi.org/10.1126/science.1244693), which can be downloaded
from
<https://earthenginepartners.appspot.com/science-2013-global-forest/download_v1.6.html>.
To cover all of Colombia, the granules needed are: 20N, 90W; 20N 80W; 10N 90W;
10N 80W; 10N 70W; 0N 90W; 0N 80W. The data sets needed are `lossyear`,
`datamask`, and `treecover2000`. The granules for each data set should be placed
into corresponding subfolders within the `data/gfc/` directory, e.g. the file
`Hansen_GFC-2018-v1.6_lossyear_20N_090W` is expected in `data/gfc/lossyear/`.

### Protected areas
Shapes for Colombian protected areas were obtained from
<http://www.parquesnacionales.gov.co/portal/es/servicio-al-ciudadano/datos-abiertos/>.
The corresponding shapefile is expected in `data/colombia/runap2`. The
pre-extension shape of the *Serranía de Chiribiquete* was obtained from a
[different data
set](http://datosabiertos.esri.co/datasets/d4d80793ff604f7aa153f3cecbe0757e_0/data?geometry=-158.757%2C-10.476%2C8.938%2C19.819).
This data set is expected to reside in
`data/colombia/Parques_Nacionales_Naturales_de_Colombia-shp`.

### Administrative boundaries
The borders of Colombia are taken from version `3.6` of the *Database of Global
Administrative Boundaries* (GADM), which can be downloaded at
<https://biogeo.ucdavis.edu/data/gadm3.6/gpkg/gadm36_COL_gpkg.zip>. The
GeoPackage is expected in `data/colombia/`.
