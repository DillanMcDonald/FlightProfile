#Dillan McDonald
#Engine Calc Ver 1.0

import Flight_Profile
import math
import mpmath

#Declaration of the parameters
    #Starting with the Lists
thrust = [] #list for the thrustcurve storage

    #The constants
Ratio_of_Specific_heats = 1.16

    #User Defined Parameters

    #Calculated Parameters

def main():




def simple_fuelgrain_calc_old(fuelgrain_mass, fuelgrain_casing_ID, htpb_percentage):  # input in kg,m,and decimal     also from previous tests 88% HTPB to 12% Papi Cures best
    # defined Constants for the fuel density
    HTPB_density = 930  # kg/m^3
    Curative_density = 1234  # kg/m^3 This is Papi 94
    Mass_of_HTPB_needed = fuelgrain_mass * htpb_percentage  # kg
    Mass_of_Curative_needed = fuelgrain_mass * (1 - htpb_percentage)  # kg
    volume_of_HTPB = Mass_of_HTPB_needed / HTPB_density  # m^3
    volume_of_Curative = Mass_of_Curative_needed / Curative_density  # m^3
    total_fuelgrain_volume = volume_of_HTPB * volume_of_Curative  # m^3

if __name__ == '__main__':
    main()