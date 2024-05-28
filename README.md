# atlantis_tools_data_processing
Misc. scripts I have developed to process data from Atlantis

Directory extract_biomass_from_nc_out:

R script for processing biomass data from the nc output file. The self-made functions are based on functions from ReactiveAtlantis, and can return either: total biomass, SSB, mature biomass, or harvestable biomass. Additionally, they can return biomass over time, or biomass-at-age over time, biomass-by-polygon over time, or biomass-at-age-at-polygon over time. 

Directory visualize_predator_prey_overlap:

R_example_Atlantis_spatialmaps_predatordistribution_against_preycomposition.R
(needs the following source - R_tools_from_ReactiveAtlantis.R)
- This script is for making a map the shows the spatial distribution of a predator group at time step t, against the polygon-specific pie charts showing the prey biomass composition.
![image](https://github.com/hollyannperryman/atlantis_tools_data_processing/assets/45412684/70c4e114-94fa-468f-9e3a-c2305770b0f0)
- I wanted to eventually make the image display the diet interactions, but I just never got around to it yet
