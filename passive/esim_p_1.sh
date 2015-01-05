#!/bin/bash

python spom/simulator.py data/forest_cover_2000_13.csv com_D.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comD_2000_s500_1_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_D.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comD_2000_s500_2_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_D.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comD_2000_s500_3_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_G.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comG_2000_s500_1_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_G.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comG_2000_s500_2_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_G.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comG_2000_s500_3_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_I.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comI_2000_s500_1_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_I.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comI_2000_s500_2_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_I.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comI_2000_s500_3_  --threads=8

# ---

python spom/simulator.py data/forest_cover_2000_13.csv com_D.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comD_2000_s500_st_1_1_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_D.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comD_2000_s500_st_1_2_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_D.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comD_2000_s500_st_1_3_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_G.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comG_2000_s500_st_1_1_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_G.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comG_2000_s500_st_1_2_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_G.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comG_2000_s500_st_1_3_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_I.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comI_2000_s500_st_1_1_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_I.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comI_2000_s500_st_1_2_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_I.csv --steps 500 -vv -s data/sten00.csv -od data/ -p p_e_comI_2000_s500_st_1_3_  --threads=8 --stochasticity 0.0 1.0 1
