import os

import BD_wrapper.BD_wrapper as bdw


table_name = 'Roman-Duval+09-GRS'
suffix = ''  # '_KD+GL'  # '_p_PS_zero'

b = bdw.BayesianDistance()
#  specify the path to the Bayesian Distance Estimator executable files
b.path_to_bde = '/disk1/riener/Bayesian_distance'

#  specify number of CPUs used for multiprocessing
b.use_ncpus = 30

#  specify the column names of the Galactic Longitude, Galactic Latitude and VLSR values. If there are no column names, this can also be specified in terms of the column number, i.e. colnr_lon = 0, colnr_lat = 1, etc.
b.colname_lon = 'GLON'
b.colname_lat = 'GLAT'
b.colname_vel = 'Vlsr'
# b.colname_kda = 'KDA'

#  set weight for spiral arm prior (default: 0.5)
# b.prob_sa = 0.0
#  set weight for maser parallax prior (default: 0.25)
# b.prob_ps = 0.0
#  set weight for Galactic Latitude prior (default: 1)
# b.prob_gl = 0.0
#  set weight for kinematic distance prior
# b.prob_kd = 0.0

b.path_to_input_table = '{}.dat'.format(table_name)
b.path_to_table = '{}_distance_results{}.dat'.format(table_name, suffix)

b.calculate_distances()
b.get_table_distance_max_probability()
