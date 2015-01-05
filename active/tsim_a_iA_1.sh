#!/bin/bash

python spom/simulator.py data/forest_cover_2000_13.csv com_D_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comD_2000_s2000_iA_1_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_D_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comD_2000_s2000_iA_2_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_D_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comD_2000_s2000_iA_3_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_G_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comG_2000_s2000_iA_1_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_G_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comG_2000_s2000_iA_2_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_G_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comG_2000_s2000_iA_3_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_I_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comI_2000_s2000_iA_1_  --threads=8

python spom/simulator.py data/forest_cover_2000_13.csv com_I_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comI_2000_s2000_iA_2_  --threads=8


python spom/simulator.py data/forest_cover_2000_13.csv com_I_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comI_2000_s2000_iA_3_  --threads=8

# - - - 

python spom/simulator.py data/forest_cover_2000_13.csv com_D_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comD_2000_s2000_iA_st_1_1_1_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_D_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comD_2000_s2000_iA_st_1_1_2_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_D_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comD_2000_s2000_iA_st_1_1_3_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_G_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comG_2000_s2000_iA_st_1_1_1_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_G_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comG_2000_s2000_iA_st_1_1_2_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_G_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comG_2000_s2000_iA_st_1_1_3_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_I_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comI_2000_s2000_iA_st_1_1_1_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_I_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comI_2000_s2000_iA_st_1_1_2_  --threads=8 --stochasticity 0.0 1.0 1

python spom/simulator.py data/forest_cover_2000_13.csv com_I_a.csv --steps 2000 -vv -s data/sten00.csv -l data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  data/initial_anp_13.csv  -od data/ -p a_t_comI_2000_s2000_iA_st_1_1_3_  --threads=8 --stochasticity 0.0 1.0 1

