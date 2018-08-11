#Dillan McDonald
#Mass_Properties Ver 1.0

import math as m
import numpy as np
import scipy.integrate as integrate
import Engine_Calc as ec
import Flight_Profile as fp

#physical dimensions of the rocket
body_tube_od = 0.2094484    #m
body_tube_id = 0.2046224    #m
body_tube_length = 6.096    #m
nose_cone_length = 1.095375 #m
nose_cone_diameter = 0.2032 #m
nose_cone_surface_area = 0.9219723496  #m^2     pull from solidworks
nose_cone_thickness = .125     #figure out the average thickness of the number of layups planned for useage
fin_root = 0.762     #m
fin_span = 0.381     #m
fin_tip = 0.4572     #m
fin_thickness = 0.00635   #m
fin_placement = nose_cone_length+body_tube_length-fin_root   #m  this is to the front tip of the fin

#Run the other scripts for stored values
ec.simple_fuelgrain_calc_old(fp.fuelgrain_mass,body_tube_id-0.00635, .88)   #run the engine calc script

#2D Projection Areas of the rocket
nose_cone_projected_area_calc = integrate.quad(lambda x: ((nose_cone_diameter/m.sqrt(m.pi))*m.sqrt((m.acos(1-((2*x)/(nose_cone_length))))-(m.sin(2*(m.acos(1-((2*x)/(nose_cone_length)))))/2))),0,nose_cone_length) #m^2
nose_cone_projected_area = 2 * nose_cone_projected_area_calc[0]
body_tube_projected_area = body_tube_length * body_tube_od  #m^2 assuing a 20 ft airframe
fin_projected_area_sec1 = fin_tip * (fin_root*m.sin(60))    #the rectangular part of the fin
fin_projected_area_sec2 = ((fin_root-fin_tip) * (fin_root*m.sin(60)))/2     #the triangular part of the fin

#fin surface area as a paramter for optimization
fin_area = (fin_tip*fin_span)+(((fin_root-fin_tip)*fin_span)/2)
fin_surface_area = (fin_area*2) #assuming the thickness is negligable

#volumes of the components
body_tube_volume = ((m.pi()*(body_tube_od/2)**2)-(m.pi()*(body_tube_id/2)**2))*body_tube_length
fin_volume = fin_area*fin_thickness
nose_cone_volume = nose_cone_surface_area*nose_cone_thickness
nozzle_volume = 0
tank_volume = 2*((((4/3)*m.pi()*(body_tube_od)**2))-((4/3)*m.pi()*(body_tube_id)**2))     #again, this is just the end caps
fuelgrain_volume = ec.total_fuelgrain_volume
liner_volume = 0
injector_volume = 0
valve_volume = 0
oxidiser_volume = 0   #this is the actual oxidiser required volume

#densities of the materials used in kg/m^3
body_tube_density = 2.7   #T6-6061
fin_density = 2.7 #T6-6061
nozzle_density = 2.7  #Graphite, this is from a range of 2.3 to 2.7
nose_cone_density = 0  #Carbon Fibre, but heavy , Carbon is 1800, epoxy is 1100 we use a mixing ration of 100 to 22
tank_density = 2.7    #T6-6061
fuelgrain_density = ec.total_grain_density   #HTPB and Papi 94 *pull this from the engine calc parameters
liner_density = 1.32   #Phenolic 1.24-1.32
injector_density = 2.7    #T6-6061
valve_density = 0   #this will have to be arbitrary and built from the densest valve
oxidiser_density = 750    #N20

#define CG locations in respect to the tip of the nosecone in meters
#this version will assume all CGs are coaxial
#these values are approximations based on simplified geometric forms
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

#find the CG of the vehicle from the tip of the nosecone
Dry_CG = ((body_tube_CG*body_tube_mass)+(fin_CG*fin_mass)+(nozzle_CG*nozzle_mass)+(tank_CG*tank_mass)+(liner_CG*liner_mass)+(injector_CG*injector_mass)+(valve_CG*valve_mass))/(body_tube_mass+fin_mass+nozzle_mass+tank_mass+liner_mass+injector_mass+valve_mass)
Wet_CG = ((fuelgrain_CG*fuelgrain_mass)+(oxidiser_CG*oxidiser_mass)+(body_tube_CG*body_tube_mass)+(fin_CG*fin_mass)+(nozzle_CG*nozzle_mass)+(tank_CG*tank_mass)+(liner_CG*liner_mass)+(injector_CG*injector_mass)+(valve_CG*valve_mass))/(fuelgrain_mass+oxidiser_mass+body_tube_mass+fin_mass+nozzle_mass+tank_mass+liner_mass+injector_mass+valve_mass)

#calculate the moments of inertia per component in respect to the CG for Dry
body_tube_moi = body_tube_mass*np.abs(body_tube_CG-Dry_CG)**2
fin_moi = fin_mass*np.abs(fin_CG-Dry_CG)**2
nozzle_moi = nozzle_mass*np.abs(nozzle_CG-Dry_CG)**2
tank_moi = tank_mass*np.abs(tank_CG-Dry_CG)**2
liner_moi = liner_mass*np.abs(liner_CG-Dry_CG)**2
injector_moi = injector_mass*np.abs(injector_CG-Dry_CG)**2
valve_moi = valve_mass*np.abs(valve_CG-Dry_CG)**2

#overall vehicle moi DRY
vehicle_moi = body_tube_moi+fin_moi+nozzle_moi+tank_moi+liner_moi+injector_moi+valve_moi

#find the center of pressure from the tip of the nosecone
simple_cp = ((nose_cone_projected_area*(2/3)*nose_cone_length) + (body_tube_projected_area*((body_tube_length/2)+nose_cone_length))+(fin_placement*2*(fin_projected_area_sec1+fin_projected_area_sec2)))/(((fin_projected_area_sec2+fin_projected_area_sec1)*2)+nose_cone_projected_area_calc+body_tube_projected_area)

#find the stability margin
stability_margin = (Dry_CG-simple_cp)/body_tube_od

def main():
    print(stability_margin)



if __name__ == '__main__':
    main()
