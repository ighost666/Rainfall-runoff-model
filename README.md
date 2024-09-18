# Rainfall-runoff-model

## Model description

The model is in the file 'Main.ipynb'. The equations and the explanations of parameters are given in the model function. The "model_st" function represents the two-reservoir setup and the "model_st3" function is the three_reservoir setup.

### Input
For the current model, the only input is precipitation.
### Parameters
1. Surface areas of each reservoir (soil, road, tank)
2. Soil parameters: s_h, s_w, s*, s_fc, rooting depth, porosity, infiltration rate
3. Plant parameters: maximum evapotranspiration, evapotranspiration at wilting point
4. Model parameters: water policy ratio kp, maximum height of the tank

### Output
By inputting the precipitation, the model generates:
1. water height in road reservoir and tank reservoir, soil moisture
2. The overflow from soil, tank and road reservoir, plus the total runoff


## Data
1. "df_swstch.pickle" contains the data needed for running the short-term test, from July 17th to August 20th.
2. “df_precp_10yr.pickle" contains the long-term precipitation data sourced from Lausanne (LSN) MeteoSwiss station, from 2014 to 2024.

Both datasets were processed and could be used as direct input for the model.
