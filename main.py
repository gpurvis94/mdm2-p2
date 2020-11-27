import sys
import numpy as np
import csv

def main(args):

    # x axis time values
    t_interval = 0.01
    t_min = 0 # This is inclusive
    t_max  = 1800 + t_interval # Not inclusive
    t = np.arange(0, t_max, t_interval)
    
    # constants values
    c = {}
    r = { 't' : t }
    initialise_variables(r, c)
    
    # Skip the first value in t as values for 0 are defined
    t_iterator = iter(t)
    next(t_iterator)
    
    if len(args) > 1:
        args_iter = iter(args)
        for k in args_iter:
            if k in r:
                r[k][0] = float(next(args_iter))
            elif k in c:
                c[k] = float(next(args_iter))
            else:
                input(f'idk {k }')
    
    for dt in t_iterator:
    
        # Each iteration represents a step from
        # t1 to t2 where t_interval = t2 - t1
        
        # First thing to do is calculate temperatures at t2 
        # using values from t1, T(t=t2) = T(t=t1) + dTdt * t_interval
        r['T_R'].append(r['T_R'][-1] + t_interval * dT_Rdt(r, c))
        r['T_F'].append(r['T_F'][-1] + t_interval * dT_Fdt(r, c))
        r['T_O'].append(r['T_O'][-1] + t_interval * dT_Odt(r, c))
        
        # Next we calculate the concentrations of A and B at t2
        # given feed rate and consumption from r at t1 to t2
        
        print("\nstart", t_interval)
        print(f"C_A1 {r['C_A'][-1]} C_B1 {r['C_B'][-1]}")
        
        # From time t1 to t2 where t2-t1=d:
        # C_A(t=t2) = C_A(t=t1) + (C from feed * dt) - (C from r * dt)
        # need to make sure consumption from r isn't > amount in the tank at t=t2
        C_A = r['C_A'][-1] + t_interval * ((c['F_R']/c['V_R']) * (c['C_Afeed'] - r['C_A'][-1]))
        C_B = r['C_B'][-1] + t_interval * ((c['F_R']/c['V_R']) * -r['C_B'][-1])
        print(f"Feed/drain: C_A2 {C_A} C_B2 {C_B}")
        # Set to 0 if negative
        C_B = max(C_B, 0)
        C_A = max(C_A, 0)
        # The rate of reaction consumes r (conc/s) of A and B,
        # cannot consume more of one than the other and cannot
        # decrease concentration to below 0
        adjusted_r = -min(C_A/t_interval, C_B/t_interval, -r['r'][-1])
        r['adjusted_r'].append(adjusted_r)
        print(f"r1 {r['r'][-1]} adj {adjusted_r}")
        # Final modifications to C_A C_B from r gives C_A C_B at t2
        r['C_A'].append(C_A + adjusted_r * t_interval)
        r['C_B'].append(C_B + adjusted_r * t_interval)
        
        # Next is to calc the rate of reaction at t2 using t2's values
        r['r'].append(calc_r(r, c))
        
        
        # Now we calculate the power exchange balances for t2 from t2 values
        
        # Power from Q_ER is the amount of r calculated above
        r['Q_ER'].append(300000 * (-adjusted_r) * c['V_R'])
        r['Q_RF'].append(calc_Q_RF(r, c))
        r['Q_FO'].append(calc_Q_FO(r, c))
        r['Q_OA'].append(c['sbc'] * c['emissivity'] * c['A_Oo']
                 * ((r['T_O'][-1] + 273.15)**2 + (c['T_A']**2 + 273.15))
                 * (r['T_O'][-1] - c['T_A']) * (546.3 + r['T_O'][-1] + c['T_A']))
        
        # Cooling energy loss is only if coolant is on if reactants temp over threshold temperature
        threshold_temperature = 100
        if r['T_R'][-1] > threshold_temperature:
            r['Q_FC'].append(c['m_F'] * c['c_F']
                             * (c['F_C']/c['V_F']) * (r['T_F'][-1] - c['T_C']))
        else:
            r['Q_FC'].append(0)
    
    write_results(r, c)
    
def initialise_variables(r, c):

    ## Universal constants
    c['R'] = 8.3144626  # Ideal gas constant
    c['sbc'] = 5.6703744 * 10**-8  # stefan boltzmann constant
    

    ## Reactor Dimensions (width, height, area, volume)
    
    # ALL  volumes and flow rates in m3 (= L / 1000)
    # Tank volume is NOT the volume of reactants due to air padding
    c['V_T'] = 10 / 1000
    # Tank (jacket) height
    c['L_J'] = 0.26
    # thickness of the jacket wall (assumed both walls equal thickness)
    c['width_J_wall'] = 0.001
    c['width_J_fluid'] = 0.01
    
    # We define the height and volume of the tank and calculate the
    # inside radius: distance to the INSIDE of the inner jacket wall (m)
    c['r_Ii'] = np.sqrt(c['V_T'] / (np.pi * c['L_J']))    
    
    # The volume of the tank is defined as 10L. Given that the reactants
    # have a padding of air, the volume of reactants is assumed 95% 10L
    c['V_R'] = 0.95 * c['V_T']
    c['L_R'] = c['V_R'] / (np.pi * c['r_Ii']**2)
    
    # 4 wall areas need to be calculated, 2 pi r * height
    # These areas are used for heat transfer so the area in contact
    # with fluid is required, so reactant height is used for inside
    c['A_Ii'] = 2 * np.pi * (c['r_Ii']) * (c['L_R'])
    
    c['r_Io'] = c['r_Ii'] + c['width_J_wall']
    c['A_Io'] = 2 * np.pi * (c['r_Io']) * (c['L_J'])
    
    c['r_Oi'] = c['r_Ii'] + c['width_J_wall'] + c['width_J_fluid']
    c['A_Oi'] = 2 * np.pi * (c['r_Oi']) * (c['L_J'])
    
    c['r_Oo'] = c['r_Ii'] + c['width_J_wall'] * 2 + c['width_J_fluid']
    c['A_Oo'] = 2 * np.pi * (c['r_Oo']) * (c['L_J'])
    
    # 3 volumes are required for specific heat capacity stuff, V_R included
    c['V_F'] = np.pi * c['L_J'] * (c['r_Oi']**2 - c['r_Io']**2)
    c['V_O'] = np.pi * c['L_J'] * (c['r_Oo']**2 - c['r_Oi']**2)
    
    
    ## Material properties (Heat coefficients, specific heats, densities, mass)
    
    # Convection and conduction transfer coefficients
    c['h_R'] = 2000 # todo
    c['k_J'] = 50 # todo
    c['h_F'] = 2000 # todo
    
    # todo A_Ii and L_J values - inside area different L value to outside area???
    c['U'] = 1 / ( (1 / (c['h_R'] * c['A_Ii']))
                 + (np.log(c['r_Io']/c['r_Ii']) / (2 * np.pi * c['L_J'] * c['k_J']))
                 + (1 / (c['h_F'] * c['A_Io'])) )
    
    # Specific heats
    c['c_R'] = 4184
    c['c_F'] = 4184
    c['c_O'] = 420
    
    # Densities of reactants = density fluid = density water, units are kg per m3
    c['rho_R'] = 1000
    c['rho_F'] = 1000
    c['rho_J'] = 8050 # Density of jacket wall material (steel)
    
    # Masses = density * volume
    c['m_R'] = c['rho_R'] * c['V_R']
    c['m_F'] = c['rho_F'] * c['V_F']
    c['m_O'] = c['rho_J'] * c['V_O'] # only outer jacket volume required
    
    # Steel properties for emissivity
    c['emissivity'] = 0.5 # todo
    
    
    ## Flow rates and concentrations
    
    # Temperature & flow rate of coolant (approximated as 10% jacket volume/s)
    c['T_C'] = 10
    c['F_C'] = 0.1 * c['V_F']
    
    # Concentrations of A in feed (F_R = feed rate in = out of tank)
    # todo: for now the feed rate is assumed to "replace" the tank over 30mins
    c['F_R'] = c['V_R'] / 1800
    c['C_Afeed'] = 1000 # todo, kg per m3
    
    
    ## Rates of reaction constants
    
    c['B'] = 10**10 # todo
    c['E_a'] = 51217 # calculated for a value of 300 kelvins
    
    
    ## Initial values for iteration
    
    # results string maps to results array
    r['T_R'] = [25]
    r['T_F'] = [10]
    r['T_O'] = [25]
    c['T_A'] = 25  # Air temperature is constant
    r['C_A'] = [0]
    r['C_B'] = [1000]  # todo, kg per m3
    r['adjusted_r'] = [0]
    r['r'] = [0]
    r['Q_ER'] = [0]
    c['Q_SR'] = 0 # Todo: Power from stirring is constant
    r['Q_RF'] = [calc_Q_RF(r, c)]
    r['Q_FO'] = [calc_Q_FO(r, c)]
    r['Q_FC'] = [0]
    r['Q_OA'] = [0]    


## Iteration calculations

def dT_Rdt(r, c):
    return (r['Q_ER'][-1] + c['Q_SR'] - r['Q_RF'][-1]) / (c['m_R'] * c['c_R'])

def dT_Fdt(r, c):
    return (r['Q_RF'][-1] - r['Q_FO'][-1] - r['Q_FC'][-1]) / (c['m_F'] * c['c_F'])

def dT_Odt(r, c):
    return (r['Q_FO'][-1] - r['Q_OA'][-1]) / (c['m_O'] * c['c_O'])
    
def calc_r(r, c):
    # r is in mols per m3
    return (-c['B'] * r['C_A'][-1] * r['C_B'][-1] * 
            np.e**(-c['E_a'] / (c['R']*(r['T_R'][-1] + 273.15))) )
    
# Power calculations

def calc_Q_RF(r, c):
    return c['U'] * (r['T_R'][-1] - r['T_F'][-1])
    
def calc_Q_FO(r, c):
    return c['h_F'] * c['A_Oi'] * (r['T_F'][-1] - r['T_O'][-1])


def write_results(r, c):

    with open("constants.txt", 'w') as f:
    
        f.write("Constants\n")
        for k, v in sorted(c.items()):
            f.write(f'{k:20} {v}\n')
    
        f.write("\nInitial values\n")
        for k, v in sorted(r.items()):
            f.write(f'{k:20} {v[0]}\n')
            
    with open("data.csv", 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(r.keys())
        for i in range(len(r['t'])):
            w.writerow(r[k][i] for k in r.keys())
            
    # also write data for graphing at a certain resolution
    with open("graphing.csv", 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(r.keys())
        for i in range(len(r['t'])):
            if i % 10 == 0:
                w.writerow(r[k][i] for k in r.keys())
        
        
if __name__ == '__main__':
    main(sys.argv[1:])
