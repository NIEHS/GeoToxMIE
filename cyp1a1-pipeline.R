# CYP1A1-UP Tox-Geo Monte Carlo Analysis Pipeline
#  update: 11/30/21

#  MC-toxgeo-prepare-countystack.R : Brings in the EPA/ICE in-vitro data, 
# joins with county-level data data, NATA exposure data, 
#   and prepares the cyp1a1 data by county in a stacked format

# MC-ToxGeo-Run-up-to-ivive. R : Runs the MC for the parameters up to but not 
# including the IVIVE data

# MC-ToxGeo-IVIVE-Part1.R : Creates Css (steady state plasma concentrations)
# for the 10 age and obesity groups 

# MC-ToxGeo-IVIVE-Part2.R: Creates the Css by county and chemical by matching it with
# Css values that fit the group - created from part 1


# MC-ToxGeo-Run-Risk-Measure.R : Run the dose-response and hazard quotient 
#  MC analysis