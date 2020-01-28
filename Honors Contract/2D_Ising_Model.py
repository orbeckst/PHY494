import numpy as np
import random as rand
import matplotlib.pyplot as plt

latt = 10
rlx_mc_runs = 50
mc_runs = 1000

def init_spin_arr(r,c):
    spin_arr_fill = np.ones((r,c))
    return spin_arr_fill

#Compares the state of the nieghbors and then creates an array of the nieghbors for each lattice site
def neighbors(spin_arr, latt, x, y):
    side1 = (x - 1, y)
    side2 = (x, y + 1 if (y + 1) < latt else 0)
    side3 = (x + 1 if (x + 1) < latt else 0, y)
    side4 = (x, y - 1)

    neighbors_arr = [spin_arr[side1[0], side1[1]],
                     spin_arr[side2[0], side2[1]],
                     spin_arr[side3[0], side3[1]],
                     spin_arr[side4[0], side4[1]]]

    return neighbors_arr

#This will be used the various energy of each lattice site and then calculate the total energy of the lattice
def energy_of_latt(spin_arr, latt, x, y):
    single = -1 * spin_arr[x, y] #Current lattice site energy
    sum_neighbors = sum(neighbors(spin_arr, latt, x, y)) #Energy of neighbors

    return single*sum_neighbors


def simulate():
    #Asking User for inputs
    latt = eval(input("Lattice Size: "))
    mc_runs = eval(input("Number of Monte Carlo runs: "))
    t_input = eval(input("Input Temperature Range Zero to: "))
    step_size = eval(input("Input number of Temperature Steps: "))

    #Definind array for storage of future energy and magnetization values
    mag_arr = []
    ene_arr = []

    #Defining a global variable temperature for the for loop and plotting
    global temperature
    temperature = np.linspace(0, t_input, step_size)

    for temp in temperature:

        #creating lattice that is of size latt x latt
        spin_arr = init_spin_arr(latt, latt)

        #Creating Matricies for energy and magnetism
        magnetism = np.zeros(mc_runs + rlx_mc_runs)
        energy = np.zeros(mc_runs + rlx_mc_runs)

        #Running Monte Carlo for temperature range
        for run in range(mc_runs + rlx_mc_runs):
            for i in range(latt):
                for j in range(latt):
                    ene = energy_of_latt(spin_arr, latt, i, j)

                    #Code deciding whether or not to flip that depends on random circumstance
                    if ene <= 0:
                        spin_arr[i, j] = spin_arr[i, j] * (-1)
                    elif np.exp(((-1) * ene)/temp) > rand.uniform(0,1):
                        spin_arr[i, j] = spin_arr[i, j] * (-1)
                    else:
                        continue

            #Adding to arrays of energy and magnetization and then figuring the energy/magnetization per lattice site
            energy[run] = energy_of_latt(spin_arr, latt, i, j)/ (latt**2)
            magnetism[run] = np.absolute(sum(sum(spin_arr)))/(latt**2)

        #The temperature and magnetism results for running log
        print(temp, sum(magnetism[rlx_mc_runs:])/ mc_runs)
        print(temp, sum(energy[rlx_mc_runs:])/ mc_runs)

        #The total magnetization and energy after summing each monte carlo run and getting the average value for each mc run
        mag_arr.append(sum(magnetism[rlx_mc_runs:])/ mc_runs)
        ene_arr.append(sum(energy[rlx_mc_runs:])/ mc_runs)

    return {'mag': mag_arr, 'ene': ene_arr}

output = simulate()

#plotting magnetization and energy vs temperature
temp_ene = plt.plot(temperature, output['mag'], color = 'green')
temp_mag = plt.plot(temperature, output['ene'], color = 'purple')

#To show plots
plt.show()
