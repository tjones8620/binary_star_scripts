#!/bin/bash

datapath=/mnt/mimir/jm/groupspace/pion-sims/Ostar/
python3 wind-mass.py ${datapath}/mhd-2d Ostar_mhd_hlld_d2n0128l3 ./ Ostar_mhd_hlld_d2n0128l3 2
python3 wind-mass.py ${datapath}/mhd-2d Ostar_mhd_hlld_d2n0256l3 ./ Ostar_mhd_hlld_d2n0256l3 2
python3 wind-mass.py ${datapath}/mhd-2d Ostar_mhd_hlld_d2n0512l3 ./ Ostar_mhd_hlld_d2n0512l3 2
python3 wind-mass.py ${datapath}/mhd-2d Ostar_mhd_hlld_d2n1024l3 ./ Ostar_mhd_hlld_d2n1024l3 2

python3 wind-mass.py ${datapath}/mhd-2d Ostar_mhd_hll_d2n0128l3 ./ Ostar_mhd_hll_d2n0128l3 2
python3 wind-mass.py ${datapath}/mhd-2d Ostar_mhd_hll_d2n0256l3 ./ Ostar_mhd_hll_d2n0256l3 2
python3 wind-mass.py ${datapath}/mhd-2d Ostar_mhd_hll_d2n0512l3 ./ Ostar_mhd_hll_d2n0512l3 2
python3 wind-mass.py ${datapath}/mhd-2d Ostar_mhd_hll_d2n1024l3 ./ Ostar_mhd_hll_d2n1024l3 2

python3 total_lum.py ${datapath}/mhd-2d Ostar_mhd_hlld_d2n0128l3 ./ Ostar_mhd_hlld_d2n0128l3 2
python3 total_lum.py ${datapath}/mhd-2d Ostar_mhd_hlld_d2n0256l3 ./ Ostar_mhd_hlld_d2n0256l3 2
python3 total_lum.py ${datapath}/mhd-2d Ostar_mhd_hlld_d2n0512l3 ./ Ostar_mhd_hlld_d2n0512l3 2
python3 total_lum.py ${datapath}/mhd-2d Ostar_mhd_hlld_d2n1024l3 ./ Ostar_mhd_hlld_d2n1024l3 2

python3 total_lum.py ${datapath}/mhd-2d Ostar_mhd_hll_d2n0128l3 ./ Ostar_mhd_hll_d2n0128l3 2
python3 total_lum.py ${datapath}/mhd-2d Ostar_mhd_hll_d2n0256l3 ./ Ostar_mhd_hll_d2n0256l3 2
python3 total_lum.py ${datapath}/mhd-2d Ostar_mhd_hll_d2n0512l3 ./ Ostar_mhd_hll_d2n0512l3 2
python3 total_lum.py ${datapath}/mhd-2d Ostar_mhd_hll_d2n1024l3 ./ Ostar_mhd_hll_d2n1024l3 2


python3 plot_2D_lum-mass.py Ostar_mhd_hlld_d2n0128l3
python3 plot_2D_lum-mass.py Ostar_mhd_hlld_d2n0256l3
python3 plot_2D_lum-mass.py Ostar_mhd_hlld_d2n0512l3
python3 plot_2D_lum-mass.py Ostar_mhd_hlld_d2n1024l3

python3 plot_2D_lum-mass.py Ostar_mhd_hll_d2n0128l3
python3 plot_2D_lum-mass.py Ostar_mhd_hll_d2n0256l3
python3 plot_2D_lum-mass.py Ostar_mhd_hll_d2n0512l3
python3 plot_2D_lum-mass.py Ostar_mhd_hll_d2n1024l3
