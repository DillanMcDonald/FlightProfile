#Dillan McDonald
#Engine Calc Ver 1.0

import Flight_Profile
import math as m
import mpmath as mp

#Declaration of the parameters
    #Starting with the Lists
thrust = [] #list for the thrustcurve storage

    #The constants
Ratio_of_Specific_heats = 1.16  #non-dimensional
atmospheric_pressure = 14.7 #PSI
Wt = 3.98 #Kg/s
R = 65 #ft-ls/lb-deg
Gc = 32.2 #ft/s^2
Cf = 1.35
density_of_N20 = 750 #kg/m^3

    #User Defined Parameters
thrust_target = 8000 #N
Combustion_chamber_pressure = 450 #PSI - Should this be User defined?
Pa = 0.462815421    #psi
simple_epsilon = 7
Pa_for_complex_epsilon = 9.4

    #Calculated Parameters
Pc_Mpa = 0.00689476*Combustion_chamber_pressure


def main():
    print('Main')


def conical_nozzle_calculations_old():
    global Me
    global Pt
    global At
    global Dt
    global De
    global Re
    global Rt
    global Ln
    global L1
    global Pt
    global expirimental_epsilon

    Me = m.sqrt((2/(Ratio_of_Specific_heats-1))*(((Combustion_chamber_pressure/Pa)**((Ratio_of_Specific_heats-1)/Ratio_of_Specific_heats))-1))
    Pt = Combustion_chamber_pressure*(1+(Ratio_of_Specific_heats-1)/2)**(-Ratio_of_Specific_heats/(Ratio_of_Specific_heats-1))
    At = ((thrust_target*0.224808942443)/(Cf*Combustion_chamber_pressure))
    Dt = m.sqrt((4*At)/m.pi)
    De = m.sqrt(simple_epsilon)*Dt
    Re = De/2
    Rt = Dt/2
    Ln = (Rt*(m.sqrt(simple_epsilon)-1))+(((1.5*Rt)*(m.sec(15)-1))/m.tan(15))
    L1 = 1.5*Rt*m.sin(15)
    Pt = Pc_Mpa*((1+(Ratio_of_Specific_heats-1))/2)**(-Ratio_of_Specific_heats/(Ratio_of_Specific_heats-1))
    expirimental_epsilon = 1/((((Ratio_of_Specific_heats+1)/2)**(1/(Ratio_of_Specific_heats-1)))*((Pa_for_complex_epsilon/Combustion_chamber_pressure)**(1/Ratio_of_Specific_heats))*m.sqrt(((Ratio_of_Specific_heats+1)/(Ratio_of_Specific_heats-1))*(1-(Pa_for_complex_epsilon/Combustion_chamber_pressure)**((Ratio_of_Specific_heats-1)/Ratio_of_Specific_heats))))

def Post_and_Pre_combustion_chamber_calculations(airframe_diameter):
    post_length = .75*airframe_diameter #not sure where this is from
    pre_length = .5*airframe_diameter #also not sure where this is from


def simple_fuelgrain_calc_old(fuelgrain_mass, fuelgrain_casing_ID, htpb_percentage):  # input in kg,m,and decimal     also from previous tests 88% HTPB to 12% Papi Cures best
    # defined Constants for the fuel density
    HTPB_density = 930  # kg/m^3
    Curative_density = 1234  # kg/m^3 This is Papi 94
    Mass_of_HTPB_needed = fuelgrain_mass * htpb_percentage  # kg
    Mass_of_Curative_needed = fuelgrain_mass * (1 - htpb_percentage)  # kg
    volume_of_HTPB = Mass_of_HTPB_needed / HTPB_density  # m^3
    volume_of_Curative = Mass_of_Curative_needed / Curative_density  # m^3
    total_grain_density = fuelgrain_mass/total_fuelgrain_volume #kg/m^3
    airframe_inner_diameter = fuelgrain_casing_ID*.0254 #m
    airframe_inner_radius = airframe_inner_diameter/2   #m
    total_fuelgrain_volume = volume_of_HTPB * volume_of_Curative  # m^3
    porthole_diameter = 2*De*.0254  #m
    porthole_diameter_in = porthole_diameter/.0254 #in
    porthole_area = m.pi*(porthole_diameter_in/2)**2
    length_of_grain = total_fuelgrain_volume/((m.pi/4)*((airframe_inner_diameter**2)-(porthole_diameter**2)))   #m

def tank_calculations(oxidiser_mass, airframe_ID):  #mass in kg, and ID in inches
    volume_of_n20 = oxidiser_mass/density_of_N20    #m^3
    airframe_ID_m = airframe_ID*.0254   #m
    airframe_radius_m = airframe_ID_m/2 #m
    spherical_endcap_volume = 4/3*(m.pi*airframe_radius_m**3)   #m^3
    airframe_crossectional_area = m.pi*airframe_radius_m**2 #m^2
    cylindrical_tank_lenth = (volume_of_n20-spherical_endcap_volume)/airframe_crossectional_area    #m
    exact_total_tank_length = cylindrical_tank_lenth + airframe_ID_m    #m
    cylindrical_tank_lenth_safety = volume_of_n20/airframe_crossectional_area #m
    total_tank_length_safety = cylindrical_tank_lenth_safety + airframe_ID_m


if __name__ == '__main__':
    main()