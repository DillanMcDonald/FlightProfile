#Dillan McDonald
#Flight Profile Ver 1.0

import Engine_Calc
import Parachute_Parameters

import math as m
import mpmath as mp
from matplotlib import pyplot as plt

#Declaration of the parameters
    #starting with the lists
time = []   #seconds
altitude = [] #meters
velocity = [] # m/s
thrust = [] #thrust curve data in N
coefficient_of_drag = [] #generated coefficient of Drag data
dynamic_pressure = [] #Pa
drag_force = [] #N
vehicle_mass = [] #kg
thrust_accel = []   #m/s^2
gravity_accel = []  #m/s^2
drag_accel = [] #m/s^2
net_accel = []  #m/s^2
temp = [] #ambient temperature in K
air_pressure = [] #in Pascals
air_density = [] #kg/m^3

    #The constants
specific_gas_constant = 287.05 #J/kJ*K
rho_sea_level = 1.177
simple_gravity = 9.81 #m/s^2
temp_sea_level = 299 #K
air_pressure_sea_level = 101325 #Pa

    #The User Defined Parameters
bodytube_od = 0.2094484 #meters
bodytube_id = 0.2046224 #meters
start_vehicle_dry_mass = 330.589099999 #kg
start_vehicle_propellent_mass = 317.4 #kg
of_ratio = 2.3 #should be investigated from trades
time_step = 0.1 #seconds
simple_cd = .45 #generalized Cd given with nosecone
thrust_goal = 8896.443 #N ,the number that we will design for
is_parachute = 0 #is there a parachute?
parachute_diameter = 0 #parachute diameter in meters
cd_for_parachute_design = 0.9   #this ain't no cone
parachute_deployment_altitude = 30000   #in meters

    #Propep Parameters based on first initialization
specific_impulse = 250 #m/s

    #The Calculated Parameters
vehicle_crossectional_area = m.pi*(bodytube_od/2)**2    #assuming constant
start_vehicle_wet_mass = start_vehicle_dry_mass + start_vehicle_propellent_mass #kg
burn_time = (start_vehicle_propellent_mass*specific_impulse*simple_gravity/thrust_goal) #Seconds
nitrous_oxide_mass = (start_vehicle_propellent_mass/(of_ratio+1))*of_ratio  #nitrous mass from o/f Ratio
fuelgrain_mass = (start_vehicle_propellent_mass/(of_ratio+1))   #Grain mass from o/f Ratio
propellent_mass_flow = start_vehicle_propellent_mass/burn_time #kg/s
nitrous_oxide_mass_flow = nitrous_oxide_mass/burn_time #kg/s
fuelgrain_mass_flow = fuelgrain_mass/burn_time #kg/s

print("Burn Time: ", burn_time)

    #Just Programming Parameters
past_apogee = 0 #boolean variable for past apogee
hit_ground = 0  #boolean variable for hitting the ground
apogee = 0  #storage of apogee

def main():
    global hit_ground
    count = 1
    first_pass()
    while hit_ground != 1:
        flight_loop(count)
        count = count + 1
        #print("Altitude : ", altitude[count-1])
        #print("ANet : ", net_accel[count - 1])
    print("Done")
    print("Max Altitude(m): ", max(altitude))
    print("Max Velocity(m/s): ", max(velocity))
    plt.subplot(2,2,1)
    plt.plot(time,net_accel)
    plt.title('Flight Profile')

    plt.ylabel('Net Acceleration (m/s^2)')
    #plt.xlabel('time (s)')
    plt.subplot(2,2,2)
    plt.plot(time,velocity)
    plt.ylabel('Velocity (m/s)')
    #plt.xlabel('time (s)')
    plt.subplot(2,2,3)
    plt.plot(time,altitude)
    plt.ylabel('Altitude (m)')
    plt.xlabel('time (s)')
    plt.subplot(2,2,4)
    plt.plot(time,drag_accel)
    plt.ylabel('Acceleration due to Drag (m/s^2)')
    plt.xlabel('time (s)')


    plt.show()

def first_pass():   #set the initial values at t=0 to prepare the iterative loop to commence, goes through 2 time steps
    time.append(0)
    thrust.append(0)
    coefficient_of_drag.append(simple_cd)
    dynamic_pressure.append(0)
    drag_force.append(0)
    vehicle_mass.append(float(start_vehicle_wet_mass))
    thrust_accel.append(0)
    gravity_accel.append(float(-simple_gravity))
    drag_accel.append(0)
    net_accel.append(thrust_accel[0] + gravity_accel[0] + drag_accel[0])
    velocity.append(0)
    altitude.append(0)
    temp.append(float(temp_sea_level))
    air_pressure.append(float(air_pressure_sea_level))
    air_density.append(float(rho_sea_level))

def flight_loop(count):  #normal flight process
    global past_apogee
    global hit_ground
    global vehicle_crossectional_area
    global apogee

    time.append(count*time_step)
    if count*time_step <= burn_time:
        thrust.append(float(thrust_goal))  #super simple assumed constant thrust
        vehicle_mass.append(start_vehicle_wet_mass-(propellent_mass_flow*time[count]))  #mass dispelled by propellent
    else:
        thrust.append(0)    #post burnout
        vehicle_mass.append(vehicle_mass[count-1])  #assuming constant mass flow

    coefficient_of_drag.append(simple_cd)   #assuming a constant Cd
    dynamic_pressure.append((air_density[count-1]*velocity[count-1]**2)/2)   #calculated the dynamic pressure based on previous cycle values
    if past_apogee == 1:
        drag_force.append(.5*coefficient_of_drag[count]*vehicle_crossectional_area*air_density[count-1]*velocity[count-1]**2) #calculate the current drag force from the previous velocity
        if is_parachute == 1:
            if parachute_deployment_altitude <= apogee: #ensure that the deployment altitude is not above the apogee
                if altitude[count-1] < parachute_deployment_altitude :
                    vehicle_crossectional_area = m.pi*(parachute_diameter/2)**2
                    coefficient_of_drag[count] = cd_for_parachute_design
            else:
                print("Deployment Altitude Out of Range")
    else:
        apogee = altitude[count-1]
        drag_force.append(-.5*coefficient_of_drag[count]*vehicle_crossectional_area*air_density[count-1]*velocity[count-1]**2) #calculate the current drag force from the previous velocity
    thrust_accel.append(thrust[count]/vehicle_mass[count])
    gravity_accel.append(-0.00000000006673*(5.97219E+24/(6.378*10**6+altitude[count-1])**2))    #variable gravitational acceleration based on altitude
    drag_accel.append(drag_force[count]/vehicle_mass[count])
    net_accel.append(thrust_accel[count] + gravity_accel[count] + drag_accel[count])
    velocity.append(time_step*net_accel[count]+velocity[count-1])
    altitude.append(altitude[count-1]+(velocity[count]*time_step)+(.5*net_accel[count]*time_step**2))
    temp.append(0.0000002*(altitude[count]**2)-0.008925*altitude[count]+299.9)  #not sure where I found this equation
    air_pressure.append(air_pressure_sea_level*mp.exp((-0.0001)*altitude[count]))   #or this one either
    air_density.append(air_pressure[count]/(temp[count]*specific_gas_constant)) #Equation came from a trendline from Booker T. Balloon data
    if altitude[count]<altitude[count-1]:
        past_apogee = 1
    if altitude[count] <= 0 and past_apogee == 1 :
        hit_ground = 1

if __name__ == '__main__':
    main()