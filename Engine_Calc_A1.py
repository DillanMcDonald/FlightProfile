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
    Epsilon = 5.7
    # optimum_epsilon = 1/(pow((Ratio_of_Specific_heats+1)/2,1/(Ratio_of_Specific_heats-1))*pow(Pe/Po,1/Ratio_of_Specific_heats)*pow(((Ratio_of_Specific_heats+1)/(Ratio_of_Specific_heats-1))*(1-pow(Pe/Po,(Ratio_of_Specific_heats-1)/Ratio_of_Specific_heats)),1/2))
    thrust_chamber_ideal(Pc)


def thrust_chamber_ideal(Combustion_chamber_pressure,Epsilon):
    # Declaration of the parameters

    # The constants
    Ratio_of_Specific_heats = 1.22  # non-dimensional this is for the TREL Operation

    atmospheric_pressure = 14.7  # PSI
    Wt = 3.98  # Kg/s
    R = 65  # ft-ls/lb-deg
    Gc = 32.2  # ft/s^2   gravitational constants
    Cf = 1.532  # Vacuum thrust coefficient
    density_of_N20 = 750  # kg/m^3
    density_of_N20_I = 48.2  # lbs/ft^3
    coeff_of_discharge = .6
    k = 1.7  # idk what this is
    grav_const_I = 32.2
    simple_gravity = 9.81  # m/s^2

    # User Defined Parameters
    thrust_targetlbs = 2500
    thrust_target = thrust_targetlbs / 0.2248  # 8896 #N
    Pa = 0.462815421  # psi
    Pa_for_complex_epsilon = 9.4
    # tank_pressure = 800 #psi
    time_step = 0.1  # seconds
    of_ratio = 2.35  # initial O/F ratio
    # burn_time = 100   #seconds
    # combustion_chamber_id = 7.5 #random number in inches
    cstar_effeciency = .975
    cf_effeciency = .98
    lstar = 42  # high end
    Po = Combustion_chamber_pressure / atmospheric_pressure  #chamber pressure in atmospheres
    Pe = 1  # pressure at nozzle exit in atmospheres
    noz_length = 0.8  # this is the optimum nozzle length for a Rao Bell
    noz_inlet_angle = 45  # degrees

    # Calculated Parameters
    ispObj = CEA_Obj(oxName='LOX', fuelName='RP1')
    optimum_epsilon = ispObj.get_eps_at_PcOvPe(Pc=Combustion_chamber_pressure, MR=of_ratio, PcOvPe=(Combustion_chamber_pressure / atmospheric_pressure))
    specific_impulse = ispObj.get_Isp(Pc=Combustion_chamber_pressure, MR=of_ratio, eps=optimum_epsilon)
    # specific_impulse = 262.4
    cstarcea = ispObj.get_Cstar(Pc=Combustion_chamber_pressure, MR=of_ratio) * cstar_effeciency

    total_impulse_base11 = 889600
    start_vehicle_propellent_mass = 317.4  # kg
    total_impulsem = 1000000

    while total_impulsem > total_impulse_base11:
        burn_time = (start_vehicle_propellent_mass * specific_impulse * simple_gravity / thrust_target)  # Seconds
        # burn_time = 165
        total_impulse = thrust_target * .22481 * burn_time  # the odd number is Newton to lbf conversion
        total_impulsem = thrust_target * burn_time  # Requirement for base 11: 889600
        start_vehicle_propellent_mass = start_vehicle_propellent_mass - .01

    propellent_mass_flow = start_vehicle_propellent_mass / burn_time  # kg/s
    Pc_Mpa = 0.00689476 * Combustion_chamber_pressure
    oxidiser_mass = (start_vehicle_propellent_mass / (of_ratio + 1)) * of_ratio
    fuel_mass = start_vehicle_propellent_mass - oxidiser_mass
    print("Engine(Pressure Feed)")
    print("Thrust Target(N): ", thrust_target)
    print("Burn Time(s): ", burn_time)
    print("Specific Impusle: ", specific_impulse)
    print("\n\n\n\n\n\n")
    print("Thrust Chamber(Tubular Wall construction regeneratively cooled by fuel)")

    Tt = ispObj.get_Tcomb(Pc=Combustion_chamber_pressure, MR=of_ratio) #Output is Kelvin
    Me = m.sqrt((2/(Ratio_of_Specific_heats-1))*(((Combustion_chamber_pressure/Pa)**((Ratio_of_Specific_heats-1)/Ratio_of_Specific_heats))-1))
    Me_CEA = ispObj.get_Chamber_SonicVel( Pc=Combustion_chamber_pressure, MR=of_ratio, eps=optimum_epsilon)
    Pt = Combustion_chamber_pressure*(1+(Ratio_of_Specific_heats-1)/2)**(-Ratio_of_Specific_heats/(Ratio_of_Specific_heats-1)) #psi
    At = ((thrust_target*0.224808942443)/(Cf*cf_effeciency*Combustion_chamber_pressure))  #Area of the throat   need to confirm
    Dt = m.sqrt((4*At)/m.pi)    #inches
    De = m.sqrt(optimum_epsilon)*Dt  #inches
    Re = De/2   #inches
    Rt = Dt/2   #inches
    #Ln = (Rt*(m.sqrt(optimum_epsilon)-1))+(((1.5*Rt)*(mp.sec(15)-1))/m.tan(15))
    Ln = noz_length*((m.sqrt(optimum_epsilon)-1)*Rt/m.tan(m.radians(15)))
    L1 = 1.5*Rt*m.sin(15)
    #cstar = (Combustion_chamber_pressure * At)/(0.0685218*start_vehicle_propellent_mass/burn_time) #home calc cstar
    cstar = cstarcea*cstar_effeciency #their number is obviously better
    #cstarm = cstar*0.3048
    #Pt = Pc_Mpa*((1+(Ratio_of_Specific_heats-1))/2)**(-Ratio_of_Specific_heats/(Ratio_of_Specific_heats-1))

    #for mass flow calc to double check thrust this stuff is from https://www.grc.nasa.gov/www/k-12/rocket/rktthsum.html
    #mdot = ((At*Combustion_chamber_pressure)/m.sqrt(Tt))*m.sqrt((Ratio_of_Specific_heats/R))*((Ratio_of_Specific_heats+1)/2)**(-(Ratio_of_Specific_heats+1)/2*(Ratio_of_Specific_heats-1))
    Me = ispObj.get_MachNumber(Pc=Combustion_chamber_pressure, MR=of_ratio, eps=optimum_epsilon)
    Te = ((1 + Me**2 * (Ratio_of_Specific_heats-1)/2)**-1 )* Tt
    #Ve = Me*m.sqrt(Ratio_of_Specific_heats*R*Te)
    #exit_gamma = ispObj.get_IvacCstrTc_exitMwGam(Pc=Combustion_chamber_pressure, MR=of_ratio, eps=optimum_epsilon)[4]
    #print(Ve)
    pe = (ispObj.get_PcOvPe(Pc=Combustion_chamber_pressure, MR=of_ratio, eps=optimum_epsilon)/Combustion_chamber_pressure)**(-1)
    #pe = ((1 + Me ** 2 * (exit_gamma - 1) / 2) ** -(exit_gamma / (exit_gamma - 1)) )* Pt
    Ae = m.pi*(De/2)**2
    #thrust_out = mdot*2.20462*Ve*0.3048 + (pe-Pa)*Ae

    #New Exit Velocity Equation
    Rstar = 8314.4621 #J/kmol-K
    Pc = Combustion_chamber_pressure*0.068046
    Tc = ispObj.get_Tcomb(Pc=Combustion_chamber_pressure, MR=of_ratio)#Temperture of Combustion Chamber
    M = 21.40#average molecular weight of the exhaust gases
    Pe1=1
#    Ve = m.sqrt(((2*Ratio_of_Specific_heats)/(Ratio_of_Specific_heats-1))*((Rstar*Tc)/(M))*(1-(Pe1/Pc)**((Ratio_of_Specific_heats-1)/Ratio_of_Specific_heats)))
#    Ve = Ve*3.28084 #convert to ft/s
    #R_imp = ispObj.
#    print("Combustion Chamber Temperature: ",Tc )
    Ve = m.sqrt(((2*grav_const_I*Ratio_of_Specific_heats)/(Ratio_of_Specific_heats-1))*R*Tc*(1-(atmospheric_pressure/Combustion_chamber_pressure)**((Ratio_of_Specific_heats-1 )/Ratio_of_Specific_heats)))
#    print(Ve)

    mdot = (thrust_targetlbs*grav_const_I)/Ve #- (pe-Pa)*Ae)/Ve
    mdot2 = 1941 +827
    At2 = mdot2 /(Pt*m.sqrt((R*Tt)/(Ratio_of_Specific_heats*grav_const_I)))

    mdot3 = ((At*Pt)/m.sqrt(Tt))*m.sqrt(Ratio_of_Specific_heats/R)*((Ratio_of_Specific_heats+1)/2)**-((Ratio_of_Specific_heats+1)/(2*(Ratio_of_Specific_heats-1)))
    #       (in^2 psi /     k )                                 ?
    #print(mdot)
    #print(mdot2)
    #print(mdot3)


    #Combustion Chamber sizing the Alt Way
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
    c_surface_area = 2 * Lc * m.sqrt(m.pi * CTR * At) + (1 / m.sin((noz_inlet_angle*m.pi)/180)) * (CTR - 1) * At  # H&H Eqn.4.6
    Te_f = (Te -273.15)*(9/5)+32
    print("Thrust Target(lbs): ",thrust_targetlbs)
    print("\n")
    print("Combustion Chamber Pressure: ", Combustion_chamber_pressure)
    print("Oxidiser flow Rate(lb/s): ",(mdot/(of_ratio+1))*of_ratio)
    print("Fuel flow Rate(lb/s): ",(mdot/(of_ratio+1))*1)
    print("O/F ratio: ",of_ratio)
    print("Cstar Efficiency: ",cstar_effeciency)
    print("cstar CEA(ft/s): ",cstarcea)
    #print("cstar CEA(m/s): ", cstarcea*0.3048)
    print("Cf Efficiency: ", cf_effeciency)
    print("Cf: ",Cf)
    print("Contraction Ratio: ",CTR)
    print("Expansion Ratio: ",optimum_epsilon)
    print("Throat Area(in^2): ",At)
    print("Throat Area(in^2): ",At2)
    print("lstar: ",lstar)
    print("Mass Flow (lb/s): ",mdot)
    print("Chamber Volume: ", Vc)
    print("Chamber to Throat Area Ratio: ",CTR)
    print("Chamber Length(In): ", Lc )
    print("Chamber Diameter(In): ",Dc) #this one includes the nozzle inlet
    print("Chamber Surface Area(In^2): ", c_surface_area )

if __name__ == '__main__':
    main()