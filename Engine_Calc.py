#Dillan McDonald
#Engine Calc Ver 1.0

import Flight_Profile
import math as m
import mpmath as mp
from matplotlib import pyplot as plt

#Declaration of the parameters
    #Starting with the Lists
thrust = [] #list for the thrustcurve storage
m_dot = [] #total propellant mass flow
m_dot_fuel = [] #mass flow of fuel
m_dot_oxidiser = [] #mass flow of oxidiser
m_dot_final = []    # a cross check measure
r = [] #o/f or mixture ratio
mf = [] #total fuel mass burned
time = [] #s
G_sub_o = []    #Free stream propellant mass velocity
r_dot = []  #regression rate
fuel_mass = []  #the mass of fuelgrain burned
port_area = []
port_radius = []    #document the port radius through the burn
throat_area = []    #this will be for nozzle erosion

    #The constants
Ratio_of_Specific_heats = 1.22  #non-dimensional this is for the TREL Operation
atmospheric_pressure = 14.7 #PSI
Wt = 3.98 #Kg/s
R = 65 #ft-ls/lb-deg
Gc = 32.2 #ft/s^2   gravitational constants
Cf = 1.35   #Vacuum thrust coefficient
density_of_N20 = 750 #kg/m^3
density_of_N20_I = 48.2 #lbs/ft^3
coeff_of_discharge = .6
k = 1.7 #idk what this is
grav_const_I = 32.2

    #User Defined Parameters
thrust_target = 8896 #N
Combustion_chamber_pressure = 500 #PSI
Pa = 0.462815421    #psi
Pa_for_complex_epsilon = 9.4
tank_pressure = 800 #psi
time_step = 0.1 #seconds
of_ratio = 2.3 #initial O/F ratio
burn_time = 87.498   #seconds
#combustion_chamber_id = 7.5 #random number in inches
specific_impulse = 250  #input from propep

Po = Combustion_chamber_pressure/atmospheric_pressure #chamber pressure in atmospheres
Pe = 1#pressure at nozzle exit in atmospheres
noz_length = 0.8 #this is the optimum nozzle length for a Rao Bell

    #Calculated Parameters
Pc_Mpa = 0.00689476*Combustion_chamber_pressure
total_impulse = thrust_target*.22481 * burn_time    #the odd number is Newton to lbf conversion
optimum_epsilon = 1/(pow((Ratio_of_Specific_heats+1)/2,1/(Ratio_of_Specific_heats-1))*pow(Pe/Po,1/Ratio_of_Specific_heats)*pow(((Ratio_of_Specific_heats+1)/(Ratio_of_Specific_heats-1))*(1-pow(Pe/Po,(Ratio_of_Specific_heats-1)/Ratio_of_Specific_heats)),1/2))

print("Optimum Expansion Ratio: ",optimum_epsilon)

    #just declared global storage
total_grain_density = 0
total_fuelgrain_volume = 0
length_of_grain = 0
porthole_area = 0
At = 0


def main():
    conical_nozzle_calculations_old()
    #tank_calculations()

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
    At = ((thrust_target*0.224808942443)/(Cf*Combustion_chamber_pressure))  #Area of the throat
    Dt = m.sqrt((4*At)/m.pi)    #inches
    De = m.sqrt(optimum_epsilon)*Dt  #inches
    Re = De/2   #inches
    Rt = Dt/2   #inches
    #Ln = (Rt*(m.sqrt(optimum_epsilon)-1))+(((1.5*Rt)*(mp.sec(15)-1))/m.tan(15))
    Ln = noz_length*((m.sqrt(optimum_epsilon)-1)*Rt/m.tan(m.radians(15)))
    L1 = 1.5*Rt*m.sin(15)
    #Pt = Pc_Mpa*((1+(Ratio_of_Specific_heats-1))/2)**(-Ratio_of_Specific_heats/(Ratio_of_Specific_heats-1))
    print("Throat radius (In): ", Rt)
    print("Exit radius(In): ",Re)
    #print("Pt : ",Pt)
    #print("Me : ",Me)

def tank_calculations(oxidiser_mass, airframe_ID):  #mass in kg, and ID in inches
    global volume_of_n20
    global exact_total_tank_length
    global total_tank_length_safety

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