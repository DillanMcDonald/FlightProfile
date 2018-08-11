#Dillan McDonald
#Engine Calc Ver 1.0

import Flight_Profile
import math as m
import mpmath as mp

#Declaration of the parameters
    #Starting with the Lists
thrust = [] #list for the thrustcurve storage
m_dot = [] #total propellant mass flow
m_dot_fuel = [] #mass flow of fuel
m_dot_oxidiser = [] #mass flow of oxidiser
r = [] #o/f or mixture ratio
mf = [] #total fuel mass burned
time = [] #s
G_sub_o = []    #Free stream propellant mass velocity
r_dot = []  #regression rate

    #The constants
Ratio_of_Specific_heats = 1.16  #non-dimensional
atmospheric_pressure = 14.7 #PSI
Wt = 3.98 #Kg/s
R = 65 #ft-ls/lb-deg
Gc = 32.2 #ft/s^2   gravitational constants
Cf = 1.35   #Vacuum thrust coefficient
density_of_N20 = 750 #kg/m^3
density_of_N20_I = 48.2 #lbs/ft^3
coeff_of_discharge = .6
k = 1.7 #not sure what this value is

    #User Defined Parameters
thrust_target = 8000 #N
Combustion_chamber_pressure = 450 #PSI - Should this be User defined?
Pa = 0.462815421    #psi
simple_epsilon = 7
Pa_for_complex_epsilon = 9.4
tank_pressure = 800 #psi
dp_wanted = 450 #psi
number_of_plates = 2
needed_dp_injector = 15 #psi
weight_flow = 7.6 #Not sure where this value came from
oriface_diameter = 0.0625   #this is the previous smallest value we could manufacture (inches)
c_star = 5273.921
time_step = 0.1 #seconds
of_ratio = 7 #initial O/F ratio
initial_port_size = 3.8829  #inches this is diameter
number_of_ports = 1 #this is in total
fuelgrain_diameter = 7.5    #inches, this needs to be updated with the phenolic liner ID
burn_time = 20.52   #seconds

    #Calculated Parameters
Pc_Mpa = 0.00689476*Combustion_chamber_pressure

    #just declared global storage
total_grain_density = 0
total_fuelgrain_volume = 0
length_of_grain = 0
porthole_area = 0
At = 0


def main():
    conical_nozzle_calculations_old()
    Post_and_Pre_combustion_chamber_calculations()
    simple_fuelgrain_calc_old()
    tank_calculations()
    injection_design_based_on_dP()
    diffusion_plate_based_on_dp()

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
    Dt = m.sqrt((4*At)/m.pi)
    De = m.sqrt(simple_epsilon)*Dt
    Re = De/2
    Rt = Dt/2
    Ln = (Rt*(m.sqrt(simple_epsilon)-1))+(((1.5*Rt)*(mp.sec(15)-1))/m.tan(15))
    L1 = 1.5*Rt*m.sin(15)
    Pt = Pc_Mpa*((1+(Ratio_of_Specific_heats-1))/2)**(-Ratio_of_Specific_heats/(Ratio_of_Specific_heats-1))
    expirimental_epsilon = 1/((((Ratio_of_Specific_heats+1)/2)**(1/(Ratio_of_Specific_heats-1)))*((Pa_for_complex_epsilon/Combustion_chamber_pressure)**(1/Ratio_of_Specific_heats))*m.sqrt(((Ratio_of_Specific_heats+1)/(Ratio_of_Specific_heats-1))*(1-(Pa_for_complex_epsilon/Combustion_chamber_pressure)**((Ratio_of_Specific_heats-1)/Ratio_of_Specific_heats))))

def Post_and_Pre_combustion_chamber_calculations(airframe_diameter):    #airframe diameter is actually the ID of the combustion chambers - so it needs the liners
    post_length = .75*airframe_diameter #not sure where this is from
    pre_length = .5*airframe_diameter #also not sure where this is from


def simple_fuelgrain_calc_old(fuelgrain_mass, fuelgrain_casing_ID, htpb_percentage):  # input in kg,m,and decimal     also from previous tests 88% HTPB to 12% Papi Cures best
    global total_grain_density
    global total_fuelgrain_volume
    global length_of_grain
    global porthole_area
    # defined Constants for the fuel density
    HTPB_density = 930  # kg/m^3
    Curative_density = 1234  # kg/m^3 This is Papi 94
    Mass_of_HTPB_needed = fuelgrain_mass * htpb_percentage  # kg
    Mass_of_Curative_needed = fuelgrain_mass * (1 - htpb_percentage)  # kg
    volume_of_HTPB = Mass_of_HTPB_needed / HTPB_density  # m^3
    volume_of_Curative = Mass_of_Curative_needed / Curative_density  # m^3
    airframe_inner_diameter = fuelgrain_casing_ID*.0254 #m
    airframe_inner_radius = airframe_inner_diameter/2   #m
    total_fuelgrain_volume = volume_of_HTPB * volume_of_Curative  # m^3
    total_grain_density = fuelgrain_mass/total_fuelgrain_volume #kg/m^3
    porthole_diameter = 2*De*.0254  #m
    porthole_diameter_in = porthole_diameter/.0254 #in
    porthole_area = m.pi*(porthole_diameter_in/2)**2
    length_of_grain = total_fuelgrain_volume/((m.pi/4)*((airframe_inner_diameter**2)-(porthole_diameter**2)))   #m

def fuelgrain_calc_new():
    global time
    global  m_dot
    global r
    global m_dot_fuel
    global m_dot_oxidiser
    global G_sub_o
    global r_dot
    global length_of_grain
    c = 1

    #set the initial values
    time.append(0)
    m_dot.append(Gc * Combustion_chamber_pressure*At / (.95*c_star))   #assuming 95% efficiency at best
    r.append(of_ratio)
    m_dot_fuel.append(m_dot[0] / (r[0] + 1))
    m_dot_oxidiser.append(m_dot[0] - m_dot_fuel[0])
    db = fuelgrain_diameter - (initial_port_size/2)*(number_of_ports*2) #burn distance required in inches with an assumed starting port size(probably wrong, I don't know)
    G_sub_o.append(m_dot_oxidiser[0]/(number_of_ports*(m.pi()*(initial_port_size/2)**2)))
    r_dot.append(0.104*G_sub_o[0]**0.681)
    length_of_grain = (m_dot_fuel/number_of_ports)/(2*m.pi()*(initial_port_size/2)*(0.0346)*r_dot[0])
    time.append(time_step*c)

    while time[c] <= burn_time:
        #create an assumed constant m_dot_oxidiser
        m_dot_oxidiser.append(m_dot_oxidiser[c-1])
        

        #progress onto the next iteration
        c = c + 1
        time.append(time_step*c)

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

def injection_design_based_on_dP():  #weight flow in  lbs/s , dp in PSI
    #injector defined as the pipe diameter coming from the valve
    injector_D = m.sqrt((weight_flow ** 2) / (0.525 * coeff_of_discharge * (density_of_N20_I * needed_dp_injector))) # this value is in inches
    return injector_D

def injection_dp_from_D(input_diameter):    #input in inches
    global pressure_drop
    pressure_drop = ((weight_flow**2)/((input_diameter**2)*(0.525*coeff_of_discharge))/density_of_N20_I)    #pressure drop is in PSI
    return pressure_drop

def diffusion_plate_based_on_dp():  #needs to be called after the injection_design_based_on_dP
    global number_of_orifaces
    dp_wanted_local = (dp_wanted - needed_dp_injector) / number_of_plates   #figure out the local drop needed from the previous pressure drop
    number_of_orifaces = m.sqrt((3.627*k*(weight_flow**2))/(dp_wanted_local*density_of_N20_I*oriface_diameter**4))
    return number_of_orifaces

def diffusion_plate_based_on_injector_D():  #needs to be called after the injection_dp_from_D
    dp_wanted_local = (dp_wanted - pressure_drop) / number_of_plates    #figure out the local drop needed from the previous pressure drop
    number_of_orifaces = m.sqrt((3.627 * k * (weight_flow ** 2)) / (dp_wanted_local * density_of_N20_I * oriface_diameter ** 4))
    return number_of_orifaces

if __name__ == '__main__':
    main()