#Dillan McDonald
#Engine Calc Ver 1.0

from rocketcea.cea_obj import CEA_Obj
#import Flight_Profile
import math as m
import mpmath as mp
import numpy as np
from matplotlib import pyplot as plt

def main():
    Pc = 500  #Combustion Chamber Pressure in PSI
    expansion_ratio = 5.7
    # optimum_epsilon = 1/(pow((Ratio_of_Specific_heats+1)/2,1/(Ratio_of_Specific_heats-1))*pow(Pe/Po,1/Ratio_of_Specific_heats)*pow(((Ratio_of_Specific_heats+1)/(Ratio_of_Specific_heats-1))*(1-pow(Pe/Po,(Ratio_of_Specific_heats-1)/Ratio_of_Specific_heats)),1/2))
    thrust_target = 2500 #in lb force
    of_ratio = 2.35  # initial O/F ratio
    altitude = 0 #m this is MSL from Spaceport
    thrust_chamber_ideal(Pc,expansion_ratio,thrust_target,of_ratio,altitude)


def thrust_chamber_ideal(Combustion_chamber_pressure, Epsilon, thrust_targetlbs, of_ratio,altitude):
    # Declaration of the parameters

    # Constants
    R = 65  # ft-ls/lb-deg
    air_pressure_spaceport = 14.6  # PSI
    Ratio_of_Specific_heats = 1.22  # Non-dimensional

    # Updated "constants"
    atmospheric_pressure = air_pressure_spaceport*mp.exp((-0.0001)*altitude)   #output is PSI
    grav_const_M =0.00000000006673*(5.97219E+24/(6.378*10**6+altitude)**2)    #variable gravitational acceleration based on altitude
    grav_const_I = grav_const_M*3.28084 #ft/s^2
    Cf = 1.532  # Vacuum thrust coefficient

    # User Defined Parameters
    Pa = 0.462815421  # psi
    cstar_effeciency = .975
    cf_effeciency = .98
    lstar = 42
    noz_length = 0.8  # this is the optimum nozzle length for a Rao Bell
    noz_inlet_angle = 45  # degrees

    # Calculated Parameters
    thrust_target = thrust_targetlbs / 0.2248  # 8896 #N
    ispObj = CEA_Obj(oxName='LOX', fuelName='RP1')
    CEA_epsilon = ispObj.get_eps_at_PcOvPe(Pc=Combustion_chamber_pressure, MR=of_ratio, PcOvPe=(Combustion_chamber_pressure / atmospheric_pressure))
    #print("Expansion Ratio from CEA: ",CEA_epsilon)
    specific_impulse = ispObj.get_Isp(Pc=Combustion_chamber_pressure, MR=of_ratio, eps=Epsilon)
    #print("CEA ISP: ",specific_impulse)
    cstarcea = ispObj.get_Cstar(Pc=Combustion_chamber_pressure, MR=of_ratio) * cstar_effeciency

    #Mass Iterator to match Base 11 Max total impulse
    total_impulse_base11 = 889600
    start_vehicle_propellent_mass = 400  # kg
    total_impulsem = 1000000
    while total_impulsem > total_impulse_base11:
        burn_time = (start_vehicle_propellent_mass * specific_impulse * grav_const_M / thrust_target)  # Seconds
        total_impulsem = thrust_target * burn_time  # Requirement for base 11: 889600
        start_vehicle_propellent_mass = start_vehicle_propellent_mass - .01

    #Total Masses for Oxidiser and Fuel
    oxidiser_mass = (start_vehicle_propellent_mass / (of_ratio + 1)) * of_ratio #kg
    fuel_mass = (start_vehicle_propellent_mass / (of_ratio + 1)) #kg
    oxidiser_mass_lbs = oxidiser_mass*2.20462 #lbs
    fuel_mass_lbs = fuel_mass*2.20462 #lbs

    #Nozzle Sizing
    Tt = ispObj.get_Tcomb(Pc=Combustion_chamber_pressure, MR=of_ratio) #Temperature at throat Kelvin
    Pt = Combustion_chamber_pressure*(1+(Ratio_of_Specific_heats-1)/2)**(-Ratio_of_Specific_heats/(Ratio_of_Specific_heats-1)) #Pressure at Throat PSI
    At = ((thrust_target*0.224808942443)/(Cf*cf_effeciency*Combustion_chamber_pressure))  #Area of the throat
    Dt = m.sqrt((4*At)/m.pi)    #inches Diameter of throat
    De = m.sqrt(Epsilon)*Dt  #inches Diameter of Exit
    Re = De/2   #inches Radius of the exit
    Rt = Dt/2   #inches Radius of the Throat
    Ae = m.pi*(De/2)**2     #inches^2 Area of the Exit
    Ln = noz_length*((m.sqrt(Epsilon)-1)*Rt/m.tan(m.radians(15))) #inches Lenght of the nozzle
    L1 = 1.5*Rt*m.sin(m.radians(15))
    cstar = cstarcea*cstar_effeciency #their number is obviously better, but here is a hand calc
    Me = ispObj.get_MachNumber(Pc=Combustion_chamber_pressure, MR=of_ratio, eps=Epsilon)    # Mach at exit
    Me_hand = m.sqrt((2/(Ratio_of_Specific_heats-1))*(((Combustion_chamber_pressure/Pa)**((Ratio_of_Specific_heats-1)/Ratio_of_Specific_heats))-1)) #hand calc for Exit Mach
    Te = ((1 + Me**2 * (Ratio_of_Specific_heats-1)/2)**-1 )* Tt #Temperature at the exit
    Tc = ispObj.get_Tcomb(Pc=Combustion_chamber_pressure, MR=of_ratio)#Temperture of Combustion Chamber
    Ve = m.sqrt(((2*grav_const_I*Ratio_of_Specific_heats)/(Ratio_of_Specific_heats-1))*R*Tc*(1-(atmospheric_pressure/Combustion_chamber_pressure)**((Ratio_of_Specific_heats-1 )/Ratio_of_Specific_heats)))
    Pe = (ispObj.get_PcOvPe(Pc=Combustion_chamber_pressure, MR=of_ratio, eps=Epsilon)/Combustion_chamber_pressure)**(-1) #Pressure at exit PSI
    mdot = (thrust_targetlbs * grav_const_I)/Ve - ((Pe-atmospheric_pressure)*Ae)/Ve #Mass flow in lbs/s

    #Combustion Chamber sizing
    Vc = lstar*At
    Lc = np.exp(0.029*np.log(Dt*2.54)**2 + 0.47*np.log(Dt*2.54)+1.94)/ 2.54 #this is a rough approximation converted to inches
    Rc = m.sqrt((Vc/Lc)/m.pi) #this is my assumed cylinder
    err = 100

    #Iteration for the diameter of the combustion chanber
    while err > .001 :
        Dci = Rc*2 #the initial DC for the iteration
        Dc = m.sqrt((Dt**3 + (24/m.pi)*m.tan(noz_inlet_angle*m.pi/180)*Vc)/(Dci+6*m.tan(noz_inlet_angle*m.pi/180)*Lc))
        err = Dc-Dci
        Dci=Dc

    #Assorted calculated output parameters
    CTR = (m.pi*(Dc/2)**2)/At
    c_surface_area = 2 * Lc * m.sqrt(m.pi * CTR * At) + (1 / m.sin((noz_inlet_angle*m.pi)/180)) * (CTR - 1) * At  # H&H Eqn.4.6
    Te_f = (Te -273.15)*(9/5)+32

# Print Statment Hell
    #print("Engine(Pressure Fed)")
    #print("Thrust Target(N): ", thrust_target)
    #print("Burn Time(s): ", burn_time)
    #print("Specific Impusle: ", specific_impulse)
    #print("Thrust Chamber(Tubular Wall construction regeneratively cooled by fuel)")
    #print("Thrust Target(lbs): ",thrust_targetlbs)
    #print("Combustion Chamber Pressure: ", Combustion_chamber_pressure)
    #print("O/F ratio: ",of_ratio)
    #print("Cstar Efficiency: ",cstar_effeciency)
    #print("cstar CEA(ft/s): ",cstarcea)
    #print("cstar CEA(m/s): ", cstarcea*0.3048)
    #print("Cf Efficiency: ", cf_effeciency)
    #print("Cf: ",Cf)
    #print("Contraction Ratio: ",CTR)
    #print("Expansion Ratio: ",Epsilon)
    #print("lstar: ",lstar)
    print("Total Propellant Mass(lbm): ",start_vehicle_propellent_mass*2.20462)
    print("Total Mass Flow Rate(lbm/s): ",mdot)
    print("Oxidiser Flow Rate(lbm/s): ", (mdot / (of_ratio + 1)) * of_ratio)
    print("Fuel Flow Rate(lbm/s): ", (mdot / (of_ratio + 1)) * 1)
    #print("Chamber Volume: ", Vc)
    #print("Chamber to Throat Area Ratio: ",CTR)
    #print("Chamber Length(In): ", Lc )
    #print("Chamber Diameter(In): ",Dc) #this one includes the nozzle inlet
    #print("Chamber Surface Area(In^2): ", c_surface_area )



if __name__ == '__main__':
    main()