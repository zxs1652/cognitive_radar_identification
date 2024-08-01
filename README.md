# cognitive_radar_identification
Project FREEDOM 2.2.1: Cognitive Radar Identification

radarConfig.m:					set and config the radar and waveform parameters.
radarScenario_range_based.m: 	The most basic script for ranged-based radar tracking. 
radarScenario_range_based_policies_performance_fix_path.m:
								implement the cognitive radar environment.
								The cognitive behavior is for optimal power allocation strategy.
								The trahectory of the target is fixed.
radarScenario_range_based_policies_performance_rand_path.m:
								implement the cognitive radar environment.
								The trajectory of the target is randomly generated.
radarScenario_range_based_data_collection_pos_error.m:
								implement the beamforming and collect the data (power estimation and dis).
								The data is generated based on different uncertainties of the radar position.\
radarScenario_range_based_data_analysis_one_error_pos.m:
								use the	1) mutual information and AD-test
										2) casuality inference
								to make the decision of cognitive/non-cognitive radar.
								This script is for analyzing one of the realization.
radarScenario_range_based_data_analysis_all_error_pos.m:
								use the	1) mutual information and AD-test
										2) casuality inference
								to make the decision of cognitive/non-cognitive radar.
								This script is for analyzing all of Monte Carlo experiments
								and give the probability of the Type-I and Type-II error.


Extra packages:
All Scrtipts can be run on the MATLAB (version > 2023a)
You may install extra packages to fully run the scripts:
	-- Statistics and Machine Learning Toolbox (for KDE method)
	-- Optimization Toolbox (for power optimization)
	-- Parallel Computing Toolbox (for quick process when collecting data and compute mutual information)
