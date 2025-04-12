from cc3d.core.PySteppables import *
import numpy as np
import pandas as pd
import os
import random
import xml.etree.ElementTree as ET
from piffGenerator import generate_pif


path = os.path.dirname(os.path.abspath(__file__))

# Generate the pif file from the input_image21 file
generate_pif(path,"input_image21")

# Path to the parameter file sup.dat
param_file_path = os.path.join(path, 'sup.dat')

# Path to the xml file
xml_file_path =os.path.join(path, 'fortran2python.xml')

# Reading parameters from file as instructed in the file
def read_parameters(file_path):
    params = {}
    with open(file_path, 'r') as file:
        params['iseed'] = int(file.readline().strip())
        params['xsize'], params['ysize'] = map(int, file.readline().strip().split())
        params['nradio1'], params['nradio2'] = map(int, file.readline().strip().split())
        params['Jij'], params['Kappa'], params['temperature'] = map(float, file.readline().strip().split())
        params['porc_rad1'], params['porc_tang1'], params['porc_die1'] = map(float, file.readline().strip().split())
        params['porc_rad2'], params['porc_tang2'], params['porc_die2'] = map(float, file.readline().strip().split())
        params['sout'], params['sout2'] = map(int, file.readline().strip().split())
    return params

params = read_parameters(param_file_path)

# CHANGE PARAMETERS FROM SUP.DAT FILE MANUALLY OVER HERE IF NEEDED 
# best to change only the parameters that you are being tested
# params = {
#     'ise ed': 123456,
#     'xsize': 1050,
#     'ysize': 1050,
#     'nradio1': 8,
#     'nradio2': 8,
#     'Jij': 0.0,
#     'Kappa': 0.0,
#     'temperature': 0.0,
#     'porc_rad1': 0.0,
#     'porc_tang1': 0.0,
#     'porc_die1': 0.0,
#     'porc_rad2': 0.0,
#     'porc_tang2': 0.0,
#     'porc_die2': 0.0,
#     'sout': 0,
#     'sout2': 0
# }

# Derived parameters following Fortran code
# From the Fortran parameter (Lx=1050,Ly=1050,listQ=7000,nnn=16,k5=8)
# listQ here is a placeholder for n_posible cells
# 0.5 before converting to an integer for rounding
listQ = 7000
params['n_rad1'] = int(params['porc_rad1'] * listQ + 0.5)   # Number of radial fusions in the first stage
params['n_tang1'] = int(params['porc_tang1'] * listQ + 0.5) # Number of tangential fusions in the first stage
params['n_die1'] = int(params['porc_die1'] * listQ + 0.5)   # Number of cell deaths in the first stage
params['n_rad2'] = int(params['porc_rad2'] * listQ + 0.5)   # Number of total radial fusions
params['n_tang2'] = int(params['porc_tang2'] * listQ + 0.5) # Number of total tangential fusions
params['n_die2'] = int(params['porc_die2'] * listQ + 0.5)   # Number of total cell deaths 
# 'nfinal1' Represents the final number of steps in the first stage
# Calculated as the product of sout2, listQ, and the sum of the first stage percentages (porc_rad1 and porc_tang1), then rounded.
params['nfinal1'] = int(params['sout2'] * listQ * (params['porc_rad1'] + params['porc_tang1']) + 0.5)
# 'nfinal2' Represents the final number of steps in the second stage end of second stage
# Calculated as the product of sout2, listQ, and the sum of the second stage percentages 
params['nfinal2'] = int(params['sout2'] * listQ * (params['porc_rad2'] + params['porc_tang2']) + 0.5) 
# 'n_sout3' Represents the interval for triggering cell shrinkage deaths in the first stage.
# Calculated as the quotient of nfinal1 and n_die1, then rounded.
params['n_sout3'] = int(params['nfinal1'] / params['n_die1'] + 0.5)
# 'n_sout4' Represents the interval for triggering cell shrinkage deaths in the second stage.
# Calculated as the quotient of the difference between nfinal2 and nfinal1 and n_die2, then rounded.
params['n_sout4'] = int(params['nfinal2'] / params['n_die2'] + 0.5)
# 'nstep' Represents the total number of steps in the simulation, set to 120% of nfinal2. Actual END OF SIMULATION after for relaxation time nfinal2!
params['nstep'] = int(params['nfinal2'] * 1.2)

# [ ] Debugging print statements
print(params)

# Initialize random number generator from random module
random.seed(params['iseed'])

# TODO Michaels suggestion on start function PottsElmnt.ElementCC3D("Temperature id = temp_elem",{}, "10")
# Parse the XML file
tree = ET.parse(xml_file_path)
root = tree.getroot()

# Update XML elements with parameters
for potts in root.findall('Potts'):
    dimensions = potts.find('Dimensions')
    if dimensions is not None:
        dimensions.set('x', str(params['xsize']))
        dimensions.set('y', str(params['ysize']))
    
    temperature = potts.find('Temperature')
    if temperature is not None:
        temperature.text = str(params['temperature'])

    random_seed = potts.find('RandomSeed')
    if random_seed is not None:
        random_seed.text = str(params['iseed'])

# Write the updated XML to a new file
tree.write(xml_file_path)
print(f"Updated XML file created successfully at {xml_file_path}")

# Plot flag for player plots on the screen
plots = False

# Time for each plot to update
time_for_each_plot = 100

# (If true) cells from wound area should be deleted from the simulation, ablations (If False) they shirink and die, apoptosis 
disapear = False

# If recurrent damage is possible, this flag should be set to True
recurrent_damage = True

# Interval for wound area to be updated (the recuperation of time of the previus wound will be sum)
wound_intervals = 1

class fortran2pythonSteppable(SteppableBasePy):

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self,frequency)
        self.xsize = params['xsize'] # TODO The width of the simulation lattice or cell field? Assumed to be cell field size
        self.ysize = params['ysize']

        self.nradio1 = params['nradio1'] # The radius of the wound
        self.nradio2 = self.nradio1 + 25 # The radius of the wound's affected zone where fusion and dying cells are produced
        self.sout = params['sout'] # Interval for recording data or performing specific actions in the simulation
        self.sout2 = params['sout2'] #Time interval for performing fusion
        # From parameter (Lx=1050,Ly=1050,listQ=7000,nnn=16,k5=8)
        self.k5 = 8 # A fixed parameter used in distance calculations. In the Fortran code, this appears to be a constant related to cell distance metrics.
        self.nn = 16  # Number of neighbors considered. In the Fortran code, it corresponds to the nnn parameter which defines the number of offsets used for neighbor calculations (likely 16 nearest neighbors).
        
        self.initialize_variables()
        self.wounder_counter = 1
        self.wait_counter = 0

    # Each part of this method corresponds directly to initialization steps in the provided Fortran code. This ensures that the simulation starts with a consistent and defined state, similar to the original Fortran implementation.
    def initialize_variables(self):
        self.deltae = 0 # Energy difference between the current and candidate states in a Monte Carlo step. Used to decide whether a cell's state should change based on the Metropolis criterion.
        self.intercalar = 0 # A flag used to determine the type of fusion to perform (radial (0) or tangential (1)) based on the current state of the simulation.
        self.ntotal_rad = 0 # Tracks total number of radial fusions performed in the simulation.
        self.ntotal_tang = 0 # Tracks total number of tangential fusions performed in the simulation.
        self.ntotal_shrink = 0 # Tracks total number of cell shrinkage and death events in the simulation.
        self.nfinal = params['nfinal1'] # The final number of steps in the simulation. This value is updated based on the stage of the simulation (first or second).
        self.n_rad = params['n_rad1'] # Number of radial fusions allowed in the first stage of the simulation.
        self.n_tang = params['n_tang1'] # Number of tangential fusions allowed in the first stage of the simulation.
        self.n_shrink = params['n_die1'] # Number of cell shrinkage and death events allowed in the first stage of the simulation.
        self.n_sout1 = params['n_sout3'] # Interval for triggering cell shrinkage and death events in the first stage of the simulation.
        # Stores the area of cells at various distances from the wound center.
        self.airedistan = np.zeros((int(self.xsize / self.sout) + 1, self.k5))
        self.airedistan2 = np.zeros((int(self.xsize / self.sout) + 1, self.k5))
        # Stores the target area of cells at various distances from the wound center.
        self.kar = np.zeros((int(self.xsize / self.sout) + 1, self.k5))
        # Stores and keep track of the current area of each cell in the simulation. 
        self.aire = np.zeros(listQ)
        # Stores the x and y-coordinate of the center of mass for each cell  
        self.cx = np.zeros(listQ)
        self.cy = np.zeros(listQ)
        # Keeps track of neighboring cells to determine possible fusions and interactions
        self.neighlist = np.zeros((listQ, 100), dtype=int)
        # Used to keep track of connections or interactions between cells, possibly for calculating adhesion or other interactions
        self.links = np.zeros((listQ, listQ), dtype=int)
        # Used as the main lattice to keep track of the state and position of each cell
        self.tableau = np.zeros((self.xsize, self.ysize), dtype=int)
        # TODO 3d, I dont understand this?
        self.texture = np.zeros((self.xsize, self.ysize, 3))
        # Keeps track of how many radial fusions have been performed.
        self.n_radial = []
        # Keeps track of how many tangential fusions have been performed.
        self.n_tangen = []
        # Keeps track of how many cell shrinkage and death events have been performed.
        self.n_die = np.zeros(1, dtype=int)
        # Keeps track of cells that are candidates for fusion or death processes
        self.n_posible = np.zeros(1000, dtype=int)

        # Calculate nxred and nyred based on the Fortran code logic
        nxred = self.xsize - 2*self.nradio2
        nyred = self.ysize - 2*self.nradio2

        # Initialize ncen_x and ncen_y
        self.ncen_x = int(random.random() * nxred * 0.9999999) + self.nradio2
        self.ncen_y = int(random.random() * nyred * 0.9999999) + self.nradio2

    def start(self):
        """
        Called before MCS=0 while building the initial simulation
        """
        self.build_wall(self.WALL)
        
        
        for cell in self.cell_list:
            cell.dict["wound"] = False
            cell.dict["marked_for_shrinkage"] = False
            cell.dict["to_die"] = False  
            cell.dict["neighbor_wound"] = False 
            cell.dict["fusion_tang_rad"] = False    
            cell.targetVolume = cell.volume
            cell.lambdaVolume = 4.0 # This changes how much cells whant to keep their target volume
        
        if plots: # New plot windows for player
            # Plotting average target and actual volume
            self.plot_vol = self.add_new_plot_window(title='Average Target and Actual Volume',
                                                    x_axis_title='MonteCarlo Step (MCS)',
                                                    y_axis_title='Variables', x_scale_type='linear', y_scale_type='linear',
                                                    grid=True, config_options={'legend':True})
            self.plot_vol.add_plot("Avg.Tar.Vol.", style='Lines', color='Cyan', size=5)
            self.plot_vol.add_plot("Avg.Act.Vol.", style='Lines', color='Blue', size=5)
            self.plot_vol.add_plot("Med.Act.Vol.", style='Lines', color='lightblue', size=5)
            self.plot_vol.add_plot("Max.Vol.", style='Lines', color='red', size=5)
            # Plotting average target and actual surface
            self.plot_sur = self.add_new_plot_window(title='Average Target and Actual Surface',
                                                    x_axis_title='MonteCarlo Step (MCS)',
                                                    y_axis_title='Variables', x_scale_type='linear', y_scale_type='linear',
                                                    grid=True, config_options={'legend':True})
            self.plot_sur.add_plot("Avg.Tar.Sur.", style='Lines', color='lightgreen', size=5)
            self.plot_sur.add_plot("Avg.Act.Sur.", style='Lines', color='green', size=5)
            self.plot_sur.add_plot("Med.Act.Sur.", style='Lines', color='Aqua', size=5)
            self.plot_sur.add_plot("Max.Sur.", style='Lines', color='red', size=5)
            # Plotting cell count
            self.plot_count = self.add_new_plot_window(title='Cell Count',
                                                    x_axis_title='MonteCarlo Step (MCS)',
                                                    y_axis_title='Variables', x_scale_type='linear', y_scale_type='linear',
                                                    grid=True, config_options={'legend':True})
            self.plot_count.add_plot("Cell Count", style='Lines', color='Red', size=5)

        # Trying new implementation for data collection 
        self.df = pd.DataFrame(columns=['MCS', 'Cell_ID', 'Cell_TarVol', 'Cell_Vol', 'Cell_TarSur', 'Cell_Sur'])
        self.data_list = []
           
    def introduce_wound(self):
         # Calculate nxred and nyred based on the Fortran code logic
        nxred = self.xsize - 2 * self.nradio2
        nyred = self.ysize - 2 * self.nradio2

        # Randomize everytime a wound is introduced
        self.ncen_x = int(random.random() * nxred * 0.9999999) + self.nradio2
        self.ncen_y = int(random.random() * nyred * 0.9999999) + self.nradio2 
        
        # Determine the bounds of the wound area
        x_min = max(self.ncen_x - self.nradio2, 0)
        x_max = min(self.ncen_x + self.nradio2, self.xsize - 1)
        y_min = max(self.ncen_y - self.nradio2, 0)
        y_max = min(self.ncen_y + self.nradio2, self.ysize - 1)
        # Iterate over the wound area and mark cells as part of the wound
    
        for x in range(x_min, x_max + 1):
            for y in range(y_min, y_max + 1):
                if (x - self.ncen_x) ** 2 + (y - self.ncen_y) ** 2 <= self.nradio1 ** 2:
                    cell = self.cell_field[x, y, 0]  # Assuming 2D lattice, z = 0
                    if cell.type != 0 :  
                              cell.dict["wound"] = True
                              cell.dict["to_die"] = True
                              cell.type = self.WOUNDED
                elif(x - self.ncen_x) ** 2 + (y - self.ncen_y) ** 2 <= self.nradio2 ** 2:
                        cell = self.cell_field[x, y, 0]  # Assuming 2D lattice, z = 0
                        #if cell in self.cell_list_by_type(self.CELL):
                        if cell.type != 0 :                                   
                             cell.dict["neighbor_wound"] = True
                            
        neighbor_wound_count = 0
        Area_target_of_neighbor= 0
        for cell in self.cell_list_by_type(self.CELL):
            if cell.dict["neighbor_wound"]:
               neighbor_wound_count += 1
               Area_target_of_neighbor= Area_target_of_neighbor+ cell.targetVolume
          
        Area_target_of_wound= 0    
        
        for cell in self.cell_list_by_type(self.WOUNDED  ) : 
            Area_target_of_wound= Area_target_of_wound+ cell.targetVolume
            cell.dict["wound"] = True
            cell.dict["to_die"] = True
        
                    
        params['n_rad1'] = int(params['porc_rad1'] * neighbor_wound_count + 0.5)   # Number of radial fusions
        params['n_tang1'] = int(params['porc_tang1'] * neighbor_wound_count + 0.5) # Number of tangential fusions
        params['n_die1'] = int(params['porc_die1'] * neighbor_wound_count + 0.5)   # Number of cell deaths
        params['n_rad2'] = int(params['porc_rad2'] * neighbor_wound_count + 0.5)   # Number of  total radial fusions
        params['n_tang2'] = int(params['porc_tang2'] * neighbor_wound_count + 0.5) # Number of total tangential fusions 
        params['n_die2'] = int(params['porc_die2'] * neighbor_wound_count + 0.5)   # Number of total cell deaths 
        # 'nfinal1' Represents the final number of steps in the first stage
        # Calculated as, then rounded.
        params['nfinal1'] = int(params['sout2'] * neighbor_wound_count * (params['porc_rad1'] + params['porc_tang1']) + 0.5)
        # 'nfinal2' Represents the final number of steps in the second stage end of second stage
        # Calculated as the product of sout2, listQ, and the sum of the second stage percentages (porc_rad2 and porc_tang2) and rounded.
        params['nfinal2'] = int(params['sout2'] * neighbor_wound_count * (params['porc_rad2'] + params['porc_tang2']) + 0.5)

        params['n_sout3'] = int(params['nfinal1']/params['n_die1']+0.5)
        params['n_sout4'] = int(params['nfinal2']/params['n_die2']+0.5)
        
        params['suma_Area_target'] = int((Area_target_of_wound+ Area_target_of_neighbor*params['porc_die2'])/(params['n_rad2']+params['n_tang2'])+0.5)
        
        self.n_rad    = params['n_rad1']
        self.n_tang   = params['n_tang1']
        self.n_shrink = params['n_die1']
        self.n_sout1  = params['n_sout3']
        self.suma_Area_target= params['suma_Area_target'] 
        self.n_final_total=params['nfinal2']
        self.wound_intervals_1 = wound_intervals+ params['nfinal2']  #self.n_final_total
        
    # Create a DataFrame with the values
        data_test = {
            'neighbor_wound_count': [neighbor_wound_count],
            'Area_target_of_wound': [Area_target_of_wound],
            'n_rad1':[self.n_rad],
            'n_tang1':[self.n_tang],
            'nfinal2':[self.n_final_total]
            }
        test_df = pd.DataFrame(data_test)

        # Save the DataFrame to a CSV file
        #path = 'path/to/save'  # Replace with your desired path
        test_df.to_csv(os.path.join(path, 'wound_data.csv'), index=False) 

    def step(self, mcs):
        """
        Called every frequency MCS while executing the simulation
        
        :param mcs: current Monte Carlo step
        """
        # Relaxation the simulation if the final MCS is reached
        #if mcs >= params['nfinal2']:
        #    if mcs == params['nstep']:                
        #        CompuCellSetup.stop_simulation()
        #    return 
        
        if mcs==50 and self.wounder_counter != 10: #repeat the wound creation process for 10 times
            self.introduce_wound()
            self.wounder_counter += 1 
        
        # TODO Assuming that the cells will maintain their target volume
        for cell in self.cell_list_by_type(self.CELL):
             
            # Set the target volume to 0 to trigger cell death
            if cell.dict["to_die"] and cell.dict["marked_for_shrinkage"]:
                cell.targetVolume = 0.0
                cell.lambdaVolume = 10.0 # This changes how fast cells want to keep their target volume
                
            NEIGHBOR_DICT = self.get_cell_neighbor_data_list(cell).neighbor_count_by_type()
            if self.WOUNDED in NEIGHBOR_DICT.keys():
                continue # TODO Do something with those cells that are neighbors to the wounded cells
            # Here if we want to implement a way to trigger those cells to fuse with other cells around it we need to change how the methods for radial and tangential fusion work
            # We should pass a list of these cells as cadidates with a higher chance of being selected for fusion with other cells
            
        # Here is where you can control how the wounded cells will behave every MCS            
        for cell in self.cell_list_by_type(self.WOUNDED) :  
            #continue
            cell.targetVolume = 0.0
            cell.lambdaVolume = 100


        # CompuCell3D will automatically handle the redistribution of the area due to contact energies

        if mcs % time_for_each_plot == 0: # Plotting data on player's screen 
            # Data collection
            for cell in self.cell_list: # Loop through all cells in the simulation
                self.data_list.append({
                'MCS': mcs,
                'Cell_ID': cell.id,
                'Cell_TarVol': cell.targetVolume,
                'Cell_Vol': cell.volume,
                'Cell_TarSur': cell.targetSurface,
                'Cell_Sur': cell.surface})
            
            new_df = pd.DataFrame(self.data_list) 

        if plots and mcs % time_for_each_plot == 0:
            self.plot_vol.add_data_point("Avg.Tar.Vol.", mcs, new_df['Cell_TarVol'].mean())
            self.plot_vol.add_data_point("Avg.Act.Vol.", mcs, new_df['Cell_Vol'].mean())
            self.plot_vol.add_data_point("Med.Act.Vol.", mcs, new_df['Cell_Vol'].median())
            self.plot_vol.add_data_point("Max.Vol.", mcs, new_df['Cell_Vol'].max())

            self.plot_sur.add_data_point("Avg.Tar.Sur.", mcs, new_df['Cell_TarSur'].mean())
            self.plot_sur.add_data_point("Avg.Act.Sur.", mcs, new_df['Cell_Sur'].mean())
            self.plot_sur.add_data_point("Med.Act.Sur.", mcs, new_df['Cell_Sur'].median())
            self.plot_sur.add_data_point("Max.Sur.", mcs, new_df['Cell_Sur'].max())

            self.plot_count.add_data_point("Cell Count", mcs, len(self.cell_list))
 
        if mcs > 50:
        # Change parameters based on mcs for relaxtion         time
            if mcs == params['nfinal1'] + 1:
                self.n_rad = params['n_rad2']
                self.n_tang = params['n_tang2']
                self.n_shrink = params['n_die2']
                self.n_sout1 = params['n_sout4']
            
            #if mcs > 0:
            # At especific mcs intervals, perform fusion if aplicable           
            if mcs % self.sout2 == 0:
                if self.intercalar == 0: #and self.ntotal_rad < self.n_rad:
                    self.perform_radial_fusion()
                elif  self.intercalar == 1: #and self.ntotal_tang < self.n_tang:
                    self.perform_tangential_fusion()

            # At especific mcs intervals, perform cell shrinkage and death
            if mcs % self.n_sout1 == 0 and self.ntotal_shrink < self.n_shrink:
                self.perform_cell_shrinkage_death()

            # if the  wound is closed recurrent_damage=true , other wound will start immediately  
            Area_of_wound= 0       
            for cell in self.cell_list_by_type(self.WOUNDED  ) : 
                Area_of_wound= Area_of_wound+ cell.volume

            # TODO there might be recurrent damage in the same area which means cells needs to be tagged for death again        
            # Determine the bounds of the wound area
            if Area_of_wound== 0: #changed the logic so that wound appears every 50 MCS
                if self.wait_counter < 50:
                    self.wait_counter += 1
                else:
                    self.wait_counter = 0
                    if recurrent_damage and mcs > 0:
                    # if mcs % self.wound_intervals_1 == 0 and recurrent_damage and mcs > 0:
                        self.ntotal_rad= 0
                        self.ntotal_tang= 0
                        self.ntotal_shrink= 0
                        
                        for cell in self.cell_list:
                            #cell.dict["wound"] = False
                            cell.dict["marked_for_shrinkage"] = False
                            cell.dict["to_die"] = False  
                            cell.dict["neighbor_wound"] = False 
                            cell.dict["fusion_tang_rad"] = False
                
                        
                        nxred = self.xsize - 2 * self.nradio2
                        nyred = self.ysize - 2 * self.nradio2
                        self.ncen_x = int(random.random() * nxred * 0.9999999) + self.nradio2
                        self.ncen_y = int(random.random() * nyred * 0.9999999) + self.nradio2 
                        
                        x_min = max(self.ncen_x - self.nradio2, 0)
                        x_max = min(self.ncen_x + self.nradio2, self.xsize - 1)
                        y_min = max(self.ncen_y - self.nradio2, 0)
                        y_max = min(self.ncen_y + self.nradio2, self.ysize - 1)

                        # Iterate over the wound area and mark cells as part of the wound
                        for x in range(x_min, x_max + 1):
                            for y in range(y_min, y_max + 1):
                                if (x - self.ncen_x) ** 2 + (y - self.ncen_y) ** 2 <= self.nradio1 ** 2:
                                    cell = self.cell_field[x, y, 0]  # Assuming 2D lattice, z = 0
                                    if cell.type != 0 :  
                                        cell.dict["wound"] = True
                                        cell.dict["to_die"] = True
                                        cell.type = self.WOUNDED
                                elif(x - self.ncen_x) ** 2 + (y - self.ncen_y) ** 2 <= self.nradio2 ** 2:
                                    cell = self.cell_field[x, y, 0]  # Assuming 2D lattice, z = 0
                                    if cell.type != 0 : 
                                        cell.dict["neighbor_wound"] = True
                        
                        neighbor_wound_count = 0
                        Area_target_of_neighbor= 0
                        for cell in self.cell_list_by_type(self.CELL):
                            if cell.dict["neighbor_wound"]:
                                neighbor_wound_count += 1
                                Area_target_of_neighbor= Area_target_of_neighbor+ cell.targetVolume
                    
                        
                        Area_target_of_wound= 0        
                        for cell in self.cell_list_by_type(self.WOUNDED  ) : 
                            Area_target_of_wound= Area_target_of_wound+ cell.targetVolume
                            cell.dict["wound"] = True
                            cell.dict["to_die"] = True
                            
                        params['n_rad1'] = int(params['porc_rad1'] * neighbor_wound_count + 0.5)   # Number of radial fusions
                        params['n_tang1'] = int(params['porc_tang1'] * neighbor_wound_count + 0.5) # Number of tangential fusions
                        params['n_die1'] = int(params['porc_die1'] * neighbor_wound_count + 0.5)   # Number of cell deaths
                        params['n_rad2'] = int(params['porc_rad2'] * neighbor_wound_count + 0.5)   # Number of  total radial fusions
                        params['n_tang2'] = int(params['porc_tang2'] * neighbor_wound_count + 0.5) # Number of total tangential fusions 
                        params['n_die2'] = int(params['porc_die2'] * neighbor_wound_count + 0.5)   # Number of total cell deaths 
                        # 'nfinal1' Represents the final number of steps in the first stage
                        # Calculated as, then rounded.
                        params['nfinal1'] = int(params['sout2'] * neighbor_wound_count * (params['porc_rad1'] + params['porc_tang1']) + 0.5)
                        # 'nfinal2' Represents the final number of steps in the second stage end of second stage
                        # Calculated as the product of sout2, listQ, and the sum of the second stage percentages (porc_rad2 and porc_tang2) and rounded.
                        params['nfinal2'] = int(params['sout2'] * neighbor_wound_count * (params['porc_rad2'] + params['porc_tang2']) + 0.5)

                        params['n_sout3'] = int(params['nfinal1']/params['n_die1']+0.5)
                        params['n_sout4'] = int(params['nfinal2']/params['n_die2']+0.5)
                    
                        if params['n_rad2'] + params['n_tang2'] != 0:
                            params['suma_Area_target'] = int((Area_target_of_wound+ Area_target_of_neighbor*params['porc_die2'])/(params['n_rad2']+params['n_tang2'])+0.5)
                        else:
                            params['suma_Area_target'] = 0
                        
                        params['nfinal1'] = params['nfinal1'] + mcs
                        params['nfinal2'] = params['nfinal2'] + mcs
                        
                        self.n_rad    = params['n_rad1']
                        self.n_tang   = params['n_tang1']
                        self.n_shrink = params['n_die1']
                        self.n_sout1  = params['n_sout3']
                        self.suma_Area_target= params['suma_Area_target'] 
                        self.n_final_total=params['nfinal2']
                        # Create a DataFrame with the values
                        data_test = {
                        'neighbor_wound_count': [neighbor_wound_count],
                        'Area_target_of_wound': [Area_target_of_wound],
                        'n_rad1':[self.n_rad],
                        'n_tang1':[self.n_tang],
                        'nfinal2':[self.n_final_total]
                        }
            
                        test_df = pd.DataFrame(data_test)
                        test_df.to_csv(os.path.join(path, 'wound_data.csv'), index=False)    

                        self.wound_intervals_1 = wound_intervals+ params['nfinal2']  #self.n_final_total
                                
        
    def perform_radial_fusion(self):
        while True:
            candidate_cells = [cell for cell in self.cell_list_by_type(self.CELL) if cell.dict["neighbor_wound"] and not cell.dict["marked_for_shrinkage"] and not cell.dict["to_die"]]
            if not candidate_cells:
                return  # No valid candidate cells found
            cell_dummy = random.choice(candidate_cells)
            break
        self.ntotal_rad += 1

        # Find the nearest cell (nrestar) to fuse with dummy
        nrestar1 = None
        nrestar2 = None
        min_dist = float('inf')
        max_dist = 0
        while nrestar1 is None and nrestar2 is None:
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell_dummy):
                if neighbor == None:
                    continue
                else:
                    if not neighbor.dict["wound"] and not neighbor.dict["marked_for_shrinkage"] and not neighbor.dict["to_die"] and not neighbor.type == self.MEDIUM:
                        nradila_R = (cell_dummy.xCOM - self.ncen_x) ** 2 + (cell_dummy.yCOM - self.ncen_y) ** 2
                        nradila = (neighbor.xCOM - self.ncen_x) ** 2 + (neighbor.yCOM - self.ncen_y) ** 2
                        dist = abs(nradila - nradila_R)
                        if dist < min_dist:
                            min_dist = dist
                            nrestar1 = neighbor
                        if dist > max_dist:
                            max_dist = dist
                            nrestar2 = neighbor
          
        nrestar = random.choice([nrestar1, nrestar2])

        # Fuse the cells
        if nrestar.type == self.CELL:        
            self.merge_cells(nrestar, cell_dummy)
             # The number that multiplies is ad-hoc (we choose 1.5)
            cell_dummy.targetVolume = (nrestar.targetVolume + cell_dummy.targetVolume)*1.2 + self.suma_Area_target*1.5
            cell_dummy.dict["fusion_tang_rad"] = True
            
        # Update intercalar and radial fusion tracking
        self.intercalar = 1
        if nrestar and nrestar.id not in self.n_radial:
            self.n_radial.append(nrestar.id)
    
    def perform_tangential_fusion(self):
        while True:
            candidate_cells = [cell for cell in self.cell_list_by_type(self.CELL) if cell.dict["neighbor_wound"] and not cell.dict["marked_for_shrinkage"] and not cell.dict["to_die"]]
            if not candidate_cells:
                return  # No valid candidate cells found
            cell_dummy = random.choice(candidate_cells)
            break

        self.ntotal_tang += 1

        # Find the nearest and farthest cell (nrestar1 and nrestar2) to fuse with dummy
        nrestar1 = None
        nrestar2 = None
        min_dist = float('inf')
        max_dist = 0
        for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell_dummy):
            if neighbor != None and not neighbor.dict["wound"] and not neighbor.dict["marked_for_shrinkage"] and not neighbor.dict["to_die"] and not neighbor.type == self.MEDIUM:
                nradila = (neighbor.xCOM - self.ncen_x) ** 2 + (neighbor.yCOM - self.ncen_y) ** 2
                if nradila <= min_dist:
                    nrestar1 = neighbor
                    min_dist = nradila
                if nradila >= max_dist:
                    nrestar2 = neighbor
                    max_dist = nradila
        
        # Randomly choose between nrestar1 and nrestar2
        nrestar = random.choice([nrestar1, nrestar2])
 
        # Fuse the cells
        if nrestar.type == self.CELL:
            self.merge_cells( nrestar, cell_dummy)
            # The number that multiplies is ad-hoc  ( we choose 1.5)
            cell_dummy.targetVolume = (nrestar.targetVolume + cell_dummy.targetVolume)*1.5+ self.suma_Area_target*1.5
            cell_dummy.dict["fusion_tang_rad"] = True

        # Update intercalar and tangential fusion tracking
        if self.ntotal_rad < self.n_rad:
            self.intercalar = 0
        if nrestar and nrestar.id not in self.n_tangen:
            self.n_tangen = np.append(self.n_tangen, nrestar.id)

    def perform_cell_shrinkage_death(self):
        while True:
            candidate_cells = [cell for cell in self.cell_list_by_type(self.CELL) if cell.dict["neighbor_wound"] and not cell.dict["marked_for_shrinkage"] and not cell.dict["to_die"] and not cell.dict["fusion_tang_rad"]]
            if not candidate_cells:
                return  # No valid candidate cells found
            cell_dummy = random.choice(candidate_cells)
            
            # Check if the selected cell is not in the lists of cells that have undergone radial or tangential fusion
            if cell_dummy.id not in self.n_radial and cell_dummy.id not in self.n_tangen:
                break   
            
        self.ntotal_shrink += 1

        # Mark the selected cell for shrinkage
        cell_dummy.dict["marked_for_shrinkage"] = True
        cell_dummy.dict["to_die"] = True

    def on_stop(self):
        """
        Called if the simulation is stopped before the last MCS
        """
       
        new_df = pd.DataFrame(self.data_list)
        self.df = pd.concat([self.df, new_df], ignore_index=True)
        
        self.df.to_csv(os.path.join(path, 'cell_data.csv'), index=False)