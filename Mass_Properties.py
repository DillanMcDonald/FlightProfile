#Dillan McDonald
#Mass_Properties Ver 1.0

import math as m
import Engine_Calc as ec
import Flight_Profile as fp

#volumes of the components
body_tube_volume = 0
fin_volume = 0
nose_cone_volume = 0
nozzle_volume = 0
tank_volume = 0
fuelgrain_volume = 0
liner_volume = 0
injector_volume = 0
valve_volume = 0
oxidiser_volume = 0

#densities of the materials used
body_tube_density = 0   #T6-6061
fin_density = 0 #T6-6061
nozzle_density = 0  #Graphite
nose_cone_density = 0   #Carbon Fibre, but heavy
tank_density = 0    #T6-6061
fuelgrain_density = 0   #HTPB and Papi 94
liner_density = 0   #Phenolic
injector_density = 0    #T6-6061
valve_density = 0   #this will have to be arbitrary and built from the densest valve
oxidiser_density = 0    #N20

#define CG locations in respect to the tip of the nosecone in meters
#this version will assume all CGs are coaxial
body_tube_CG = 0
fin_CG = 0
nozzle_CG = 0
nose_cone_density = 0
tank_CG = 0 #in the case of GE A1 (a monocoque design) this is the summed tank end caps
fuelgrain_CG = 0
liner_CG = 0
injector_CG = 0
valve_CG = 0
oxidiser_CG = 0 #this is not updated during flight for this version

#find all the masses
body_tube_mass = body_tube_density*body_tube_volume
fin_mass = fin_density*fin_volume
nozzle_mass = nozzle_density*nozzle_volume
tank_mass = tank_density*tank_volume
fuelgrain_mass = fuelgrain_density*fuelgrain_volume
liner_mass = liner_density*liner_volume
injector_mass = injector_density*injector_volume
valve_mass = valve_density*valve_volume
oxidiser_mass = oxidiser_density*oxidiser_volume

#find the CG of the vehicle from the nosecone
Dry_CG = ((body_tube_CG*body_tube_mass)+(fin_CG*fin_mass)+(nozzle_CG*nozzle_mass)+(tank_CG*tank_mass)+(liner_CG*liner_mass)+(injector_CG*injector_mass)+(valve_CG*valve_mass))/(body_tube_mass+fin_mass+nozzle_mass+tank_mass+liner_mass+injector_mass+valve_mass)
Wet_CG = ((fuelgrain_CG*fuelgrain_mass)+(oxidiser_CG*oxidiser_mass)+(body_tube_CG*body_tube_mass)+(fin_CG*fin_mass)+(nozzle_CG*nozzle_mass)+(tank_CG*tank_mass)+(liner_CG*liner_mass)+(injector_CG*injector_mass)+(valve_CG*valve_mass))/(fuelgrain_mass+oxidiser_mass+body_tube_mass+fin_mass+nozzle_mass+tank_mass+liner_mass+injector_mass+valve_mass)

#calculate the moments of inertia per component in respect to the CG for Dry
body_tube_moi = body_tube_mass*m.abs(body_tube_CG-Dry_CG)**2
fin_moi = fin_mass*m.abs(fin_CG-Dry_CG)**2
nozzle_moi = nozzle_mass*m.abs(nozzle_CG-Dry_CG)**2
tank_moi = tank_mass*m.abs(tank_CG-Dry_CG)**2
liner_moi = liner_mass*m.abs(liner_CG-Dry_CG)**2
injector_moi = injector_mass*m.abs(injector_CG-Dry_CG)**2
valve_moi = valve_mass*m.abs(valve_CG-Dry_CG)**2

vehicle_moi = body_tube_moi+fin_moi+nozzle_moi+tank_moi+liner_moi+injector_moi+valve_moi


def main():
    print("main")



if __name__ == '__main__':
    main()
