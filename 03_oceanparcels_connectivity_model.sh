## in the bash terminal 

python 03_oceanparcels_connectivity_model.py \
--stations data/stations.csv \
--currents data/ocean_currents_surface_2020_2023.nc \
--outdir output/oceanparcels \
--days 365 \
--dt-hours 6 \
--output-hours 6 \
--cloud-side 7 \
--jitter-deg 0.08 \
--hit-radius-km 50