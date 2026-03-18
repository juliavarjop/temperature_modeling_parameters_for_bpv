# temperature_modeling_parameters_for_bpv

To run the present COMSOL Multiphysics simulations with MATLAB:

1.	Import AM1.5 Global spectrum data e.g. from https://www.pveducation.org/pvcdrom/appendices/standard-solar-spectra and save it as AM15.mat
   
2.	Download filtered datasets including experimental measurement data (panel temperature, weather, and solar radiation data) and simulated plane-of-array irradiance data for the (a) monofacial (MPV) and (b) vertical bifacial (VBPV) systems:

    -	(a) MPV_FMI_2015to2021_measdata_and_simirradiancedata.mat
    -	(b)	VBPV_TUAS_2018to2020_measdata_and_simirradiancedata_v2.mat

  - NOTE: Unfiltered datasets can be found at data repositories hosted by Finnish Meteorological Institute (DOI: 10.57707/fmi-b2share.tn0wb-as670, MPV data) and Turku University of Applied Science (https://nerc.turkuamk.fi/data-portal/, VBPV data). 

3.	Create COMSOL model for the (a) MPV and (b) VBPV panels using

    - (a) MPV_CBPV_SipanelTcomsolModel.m
    - (b)	VBPV_SipanelTcomsolModel.m

5.	Find the optimized convective heat transfer coefficient for the (a) MPV and (b) VBPV panels using
   
    - (a)	fitScript_MPV.m and calculateRMSdifference_MPV_v1.m
    - (b)	fitScript_VBPV.m and calculateRMSdifference_VBPV_v1.m

5.	Simulate the operating temperature of the (a) MPV and (b) VBPV panels in given weather conditions using optimized (and/or unoptimized) convective heat transfer coefficients
   
    - (a)	plotScript_MPV.m and calculateTmodsAtConditions_MPV_v1.m
    - (b)	plotScript_VBPV.m and calculateTmodsAtConditions_VBPV_v1.m

To determine the simulated Sandia, Faiman and PVsyst temperature model parameters:
1.	Simulate the panel operating temperatures as described above
   
2.	Store the simulated temperature data and weather conditions
   
    - (e.g. MPV_TmaxBack_fitCoeffs.mat and MPV_InputData.mat, or 
VBPV_TminBack_fitCoeffs.mat and VBPV_InputData.mat)

3.	Determine the model parameters for the (a) MPV and (b) VBPV panels using
   
    - (a)	temp_model_for_simulated_params_mpv.m
    - (b)	temp_model_for_simulated_params_vbpv.m
