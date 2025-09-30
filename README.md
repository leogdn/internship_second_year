AIS spoofing is the act of faking or manipulating signals in the Automatic Identification System (AIS), which ships use to broadcast their location, identity, and course. By spoofing, attackers can make a vessel appear somewhere it isn’t, hide its true movements, or create fake “ghost ships.”
This is often used for smuggling, illegal fishing, sanctions evasion, or cyberattacks on maritime tracking systems.

AIS spoofing can be detected by comparing the AIS data with some radar data. If we are able to associate each AIS track to a single radar track and to estimates the parameters of the radar,
then we can detect spoofing. 

This internship was focusing on an **AIS spoofing detection and track association algorithm**, my task was to create scenarios and to perform tests of the algorithm.

A scenario is a set 

The function `main_T2T.m` is evaluating the track association and the spoofing detection, the user can plot the evolution of the detection error and the association error over time. 
This function is calling the main algorithm, implemented in the function `T2TA.m`, this function was developped by a by a researcher who has co-supervised my internship. 

The folder `scenario_generation` contains functions used in order to create scenarios for testing.

The folder `test` contains functions used to evaluate the quality of the scenarios which have been generated. 
