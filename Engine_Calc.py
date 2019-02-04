#Dillan McDonald
#Engine Calc Ver 1.0

from rocketcea.cea_obj import CEA_Obj
#import Flight_Profile
import math as m
import mpmath as mp
import numpy as np
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
simple_gravity = 9.81 #m/s^2

    #User Defined Parameters
thrust_target = 11120.554070627251 #8896 #N
Combustion_chamber_pressure = 500 #PSI
Pa = 0.462815421    #psi
Pa_for_complex_epsilon = 9.4
tank_pressure = 800 #psi
time_step = 0.1 #seconds
of_ratio = 2.3 #initial O/F ratio
#burn_time = 100   #seconds
#combustion_chamber_id = 7.5 #random number in inches


#specific_impulse = 258.7  #input from propep

Po = Combustion_chamber_pressure/atmospheric_pressure #chamber pressure in atmospheres
Pe = 1#pressure at nozzle exit in atmospheres
noz_length = 0.8 #this is the optimum nozzle length for a Rao Bell
noz_inlet_angle = 15 #degrees

    #Calculated Parameters
ispObj = CEA_Obj( oxName='LOX', fuelName='RP1')

#optimum_epsilon = 1/(pow((Ratio_of_Specific_heats+1)/2,1/(Ratio_of_Specific_heats-1))*pow(Pe/Po,1/Ratio_of_Specific_heats)*pow(((Ratio_of_Specific_heats+1)/(Ratio_of_Specific_heats-1))*(1-pow(Pe/Po,(Ratio_of_Specific_heats-1)/Ratio_of_Specific_heats)),1/2))
optimum_epsilon = ispObj.get_eps_at_PcOvPe(Pc=Combustion_chamber_pressure, MR=of_ratio, PcOvPe=(Combustion_chamber_pressure/atmospheric_pressure))

specific_impulse = ispObj.get_Isp( Pc=Combustion_chamber_pressure, MR=of_ratio, eps=optimum_epsilon)
cstarcea = ispObj.get_Cstar(Pc=Combustion_chamber_pressure, MR=of_ratio)


total_impulse_base11 = 889600
start_vehicle_propellent_mass = 317.4 #kg
total_impulsem = 1000000

while total_impulsem > total_impulse_base11 :
    burn_time = (start_vehicle_propellent_mass*specific_impulse*simple_gravity/thrust_target) #Seconds
    total_impulse = thrust_target*.22481 * burn_time    #the odd number is Newton to lbf conversion
    total_impulsem = thrust_target* burn_time    #Requirement for base 11: 889600
    start_vehicle_propellent_mass = start_vehicle_propellent_mass - .01

propellent_mass_flow = start_vehicle_propellent_mass/burn_time #kg/s
Pc_Mpa = 0.00689476*Combustion_chamber_pressure
oxidiser_mass = (start_vehicle_propellent_mass/(of_ratio+1))*of_ratio
fuel_mass = start_vehicle_propellent_mass-oxidiser_mass
print("Thrust Target(N): ", thrust_target)
print("Vehicle Propellant Mass(Kg): ",start_vehicle_propellent_mass)
print("Burn Time(s): ",burn_time)
print("Total Impulse(lbs*s): ",total_impulse)
print("Total Impulse(N*s): ",total_impulsem, ", The Competition Maximum is 889600")
print("Specific Impusle: ",specific_impulse)
print("Oxidiser Mass(kg): ", oxidiser_mass)
print("Fuel Mass(kg): ", fuel_mass)
print("Propellant Mass Flow(kg/s): ", propellent_mass_flow)
print("Propellant Mass Flow(lbs/s): ", propellent_mass_flow*2.20462)
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

    Tt = ispObj.get_Tcomb(Pc=Combustion_chamber_pressure, MR=of_ratio) #Output is Kelvin
    Me = m.sqrt((2/(Ratio_of_Specific_heats-1))*(((Combustion_chamber_pressure/Pa)**((Ratio_of_Specific_heats-1)/Ratio_of_Specific_heats))-1))
    Me_CEA = ispObj.get_Chamber_SonicVel( Pc=Combustion_chamber_pressure, MR=of_ratio, eps=optimum_epsilon)
    Pt = Combustion_chamber_pressure*(1+(Ratio_of_Specific_heats-1)/2)**(-Ratio_of_Specific_heats/(Ratio_of_Specific_heats-1))
    At = ((thrust_target*0.224808942443)/(Cf*Combustion_chamber_pressure))  #Area of the throat   need to confirm
    Dt = m.sqrt((4*At)/m.pi)    #inches
    De = m.sqrt(optimum_epsilon)*Dt  #inches
    Re = De/2   #inches
    Rt = Dt/2   #inches
    #Ln = (Rt*(m.sqrt(optimum_epsilon)-1))+(((1.5*Rt)*(mp.sec(15)-1))/m.tan(15))
    Ln = noz_length*((m.sqrt(optimum_epsilon)-1)*Rt/m.tan(m.radians(15)))
    L1 = 1.5*Rt*m.sin(15)
    #cstar = (Combustion_chamber_pressure * At)/(0.0685218*start_vehicle_propellent_mass/burn_time) #home calc cstar
    cstar = cstarcea #their number is obviously better
    cstarm = cstar*0.3048
    #Pt = Pc_Mpa*((1+(Ratio_of_Specific_heats-1))/2)**(-Ratio_of_Specific_heats/(Ratio_of_Specific_heats-1))

    #for mass flow calc to double check thrust this stuff is from https://www.grc.nasa.gov/www/k-12/rocket/rktthsum.html
    mdot = ((At*Combustion_chamber_pressure)/m.sqrt(Tt))*m.sqrt((Ratio_of_Specific_heats/R))*((Ratio_of_Specific_heats+1)/2)**(-(Ratio_of_Specific_heats+1)/2*(Ratio_of_Specific_heats-1))
    Me = ispObj.get_MachNumber(Pc=Combustion_chamber_pressure, MR=of_ratio, eps=optimum_epsilon)
    Te = ((1 + Me**2 * (Ratio_of_Specific_heats-1)/2)**-1 )* Tt
    Ve = Me*m.sqrt(Ratio_of_Specific_heats*R*Te)
    #exit_gamma = ispObj.get_IvacCstrTc_exitMwGam(Pc=Combustion_chamber_pressure, MR=of_ratio, eps=optimum_epsilon)[4]

    pe = (ispObj.get_PcOvPe(Pc=Combustion_chamber_pressure, MR=of_ratio, eps=optimum_epsilon)/Combustion_chamber_pressure)**(-1)
    #pe = ((1 + Me ** 2 * (exit_gamma - 1) / 2) ** -(exit_gamma / (exit_gamma - 1)) )* Pt
    Ae = m.pi*(De/2)**2
    thrust_out = mdot*2.20462*Ve*0.3048 + (pe-Pa)*Ae
    mdot = (thrust_target*0.2248 - (pe-Pa)*Ae)/Ve
    #Combustion Chamber sizing the Alt Way
    lstar = 50 #high end
    Vc = lstar*At
    Lc = np.exp(0.029*np.log(Dt*2.54)**2 + 0.47*np.log(Dt*2.54)+1.94)/ 2.54 #this is a rough approximation converted to inches
    Rc = m.sqrt((Vc/Lc)/m.pi) #this is my assumed cylinder
    err = 100

    while err > .001 :
        Dci = Rc*2 #the initial DC for the iteration
        Dc = m.sqrt((Dt**3 + (24/m.pi)*m.tan(noz_inlet_angle*m.pi/180)*Vc)/(Dci+6*m.tan(noz_inlet_angle*m.pi/180)*Lc))
        err = Dc-Dci
        Dci=Dc

    CTR = (m.pi*(Dc/2)**2)/At

    Te_f = (Te -273.15)*(9/5)+32
    print("Exit Temperature(F): ", Te_f)
    #print("Mach number: ", Me)
    #print("Mach number CEA: ", Me_CEA)
    print("Throat radius (In): ", Rt)
    print("Throat Pressure(PSI): ",Pt) #should be about .56 combustion chamber pressure
    print("Exit radius(In): ",Re)
    print("Nozzle Length(In): ", Ln)
    #print("cstar(ft/s): ",cstar)
    #print("cstar(m/2): ",cstarm)
    print("cstar CEA(ft/s): ",cstarcea)
    print("cstar CEA(m/s): ", cstarcea*0.3048)
    print("M dot(kg/s): ",mdot)
    print("Exit Velocity(ft/s): ", Ve)
    print("lstar: ",lstar)
    print("Chamber Volume(In^3): ",Vc)
    print("Chamber Length(In): ", Lc )
    #print("Chamber Radius(In): ",Rc)
    print("Chamber Radius(In): ",Dc/2) #this one includes the nozzle inlet
    print("Chamber to Throat Area Ratio: ",CTR)
    #print("Pt : ",Pt)
    #print("Me : ",Me)
    print("Pressure Ratio PC/PE: ",Combustion_chamber_pressure/pe)
    print("Thrust Out(N): ",thrust_out)

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