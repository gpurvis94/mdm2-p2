import sys
import numpy as np
import csv


def main(args):

    # This is for if you're defining parameters in the command line, eg:
    # To set the air temp to 30 and conc of B to 100
    # $ python main.py T_A 30 C_B 100
    cmd_line_param = {}
    if len(args) > 1 and len(args) % 2 == 0:
        args_iter = iter(args)
        for k in args_iter:
            cmd_line_param[k] = float(next(args_iter))

    # x axis time values
    t_interval = 0.01
    t_min = 0    # This is inclusive
    t_max = 600  # Not inclusive
    t = np.arange(t_min, t_max + t_interval, t_interval)

    # Containers for values, constants and results
    # constants c maps a 'str' : value
    # results r maps a 'str' : [list of results]
    c = {}
    r = {'t': t}
    initialise_variables(r, c, t_max, cmd_line_param)

    # Skip the first value in t as values for 0 are defined
    t_iterator = iter(t)
    next(t_iterator)

    # Enter the iterative loop
    pump_active = False
    time_pump_active = 0
    for dt in t_iterator:

        # Each iteration represents a step from
        # t1 to t2 where t_interval = t2 - t1

        # First thing to do is calculate temperatures at t2 
        # using values from t1, T(t=t2) = T(t=t1) + dTdt * t_interval
        r['T_R'].append(r['T_R'][-1] + t_interval * dT_Rdt(r, c))
        r['T_F'].append(r['T_F'][-1] + t_interval * dT_Fdt(r, c))
        r['T_O'].append(r['T_O'][-1] + t_interval * dT_Odt(r, c))

        # Rate of reaction calcs. Not sure about this

        # Next we calculate the concentrations of A and B at t2
        # given feed rate and consumption from r at t1 to t2

        # From time t1 to t2 where t2-t1=dt: (remember r is negative)
        # C_A(t=t2) = C_A(t=t1) + (C from feed * dt) + (C from r * dt)
        # need to make sure consumption from r during dt isn't > amount in the tank at t=t2
        C_A = r['C_A'][-1] + t_interval * ((c['F_R']/c['V_R']) * (c['C_Afeed'] - r['C_A'][-1]))
        C_B = r['C_B'][-1] + t_interval * ((c['F_R']/c['V_R']) * -r['C_B'][-1])

        # Set to 0 if negative
        C_A = max(C_A, 0)
        C_B = max(C_B, 0)

        # Calc the rate of reaction at t2 after t1->t2 materials added
        r['r'].append(calc_r(r, c, C_A, C_B))

        # Cannot consume more reactants than there are in the tank
        adjusted_r = -min(C_A/t_interval, C_B/t_interval, -r['r'][-1])
        r['adjusted_r'].append(adjusted_r)

        # Final modifications to C_A C_B from r gives C_A C_B at t2
        r['C_A'].append(C_A + adjusted_r * t_interval)
        r['C_B'].append(C_B + adjusted_r * t_interval)

        # Now we calculate the power exchange balances for t2 from t2 values

        # Power from Q_ER is the amount of r calculated above
        r['Q_ER'].append(300000 * (-adjusted_r) * c['V_R'])
        r['Q_RF'].append(calc_Q_RF(r, c))
        r['Q_FO'].append(calc_Q_FO(r, c))
        r['Q_OA'].append(calc_Q_OA(r, c))

        # If there is no minimum uptime for the pump, then at the threshold
        # temperature the pump turns off and on nearly every t_interval
        # which is unrealistic for small values of t_interval. A minimum
        # time interval min_uptime is specified and this code
        if pump_active:
            # If the pump is running check how long for
            if time_pump_active > c['min_uptime']:
                # Stop the pump if temperature is within threshold and min_uptime met
                if r['T_R'][-1] < c['threshold']:
                    pump_active = False
                    time_pump_active = 0
                else:
                    time_pump_active += t_interval
            else:
                # If the pump hasn't been on long enough, continue running
                time_pump_active += t_interval
        else:
            if r['T_R'][-1] > c['threshold']:
                pump_active = True
                time_pump_active = t_interval

        # If the pump is active then calculate energy loss to coolant flow
        if pump_active:
            r['Q_FC'].append(c['m_F'] * c['c_F']
                             * (c['F_C']/c['V_F']) * (r['T_F'][-1] - c['T_C']))
        else:
            r['Q_FC'].append(0)

    write_results(r, c)


def initialise_variables(r, c, t_max, param):

    # Cheeky func for parsing param, returns "x" if no "var" in cmd line params
    get_param = lambda x, var: x if var not in param else param[var]

    # Universal constants ds[1] if 2 in ds else 5
    c['R'] = 8.3144626  # Ideal gas constant
    c['sbc'] = 5.6703744 * 10**-8  # stefan boltzmann constant
    
    # Reactor Dimensions (width, height, area, volume)
    
    # ALL  volumes and flow rates in m3 (= L / 1000)
    # Tank volume is NOT the volume of reactants due to air padding
    c['V_T'] = 10 / 1000
    # Tank (jacket) height
    c['L_J'] = get_param(0.26, 'L_J')
    # thickness of the jacket wall (assumed both walls equal thickness)
    c['width_J_wall'] = get_param(0.001, 'width_J_wall')
    c['width_J_fluid'] = get_param(0.01, 'width_J_fluid')
    
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
    
    # Material properties (Heat coefficients, specific heats, densities, mass)
    
    # Convection and conduction transfer coefficients
    c['h_R'] = get_param(2000, 'h_R')  # todo
    c['k_J'] = get_param(50, 'k_J')  # todo
    c['h_F'] = get_param(2000, 'h_F')  # todo

    # Overall heat transfer coefficient
    # todo A_Ii and L_J values - inside area different L value to outside area???
    c['U'] = 1 / ((1 / (c['h_R'] * c['A_Ii']))
                  + (np.log(c['r_Io']/c['r_Ii']) / (2 * np.pi * c['L_J'] * c['k_J']))
                  + (1 / (c['h_F'] * c['A_Io'])))
    
    # Specific heats
    c['c_R'] = get_param(4184, 'c_R')
    c['c_F'] = get_param(4184, 'c_F')
    c['c_O'] = get_param(420, 'c_O')
    
    # Densities of reactants = density fluid = density water, units are kg per m3
    c['rho_R'] = get_param(1000, 'rho_R')
    c['rho_F'] = get_param(1000, 'rho_F')
    c['rho_J'] = get_param(8050, 'rho_J')  # Density of jacket wall material (steel)
    
    # Masses = density * volume
    c['m_R'] = c['rho_R'] * c['V_R']
    c['m_F'] = c['rho_F'] * c['V_F']
    c['m_O'] = c['rho_J'] * c['V_O']  # only outer jacket volume required
    
    # Steel properties for emissivity
    c['emissivity'] = get_param(0.5, 'emissivity')  # todo
    
    # Flow rates and concentrations

    # Concentrations of A in feed (F_R = feed rate in = out of tank)
    # todo: for now the feed rate is assumed to "replace" the tank over 30mins
    c['F_R'] = c['V_R'] / t_max
    c['C_Afeed'] = get_param(1000, 'C_Afeed')  # todo, kg per m3
    
    # Temperature & flow rate of coolant (approximated as 10% jacket volume/s)
    c['T_C'] = get_param(10, 'T_C')
    c['F_C'] = 0.1 * c['V_F']

    # If there is no minimum up-time for the pump, then at the threshold
    # temperature the pump turns off and on nearly every t_interval
    # which is unrealistic for small values of t_interval.
    c['threshold'] = get_param(30, 'threshold')  # The threshold at which to turn cooling on
    c['min_uptime'] = get_param(10, 'min_uptime')  # The minimum amount of time cooling pump is on
    
    # Rates of reaction constants
    
    c['B'] = 10**4  # todo
    c['E_a'] = 51217  # calculated for a value of 300 kelvins
    
    # Initial values for iteration
    
    # results string maps to results array
    r['T_R'] = [get_param(20, 'T_R')]
    r['T_F'] = [get_param(10, 'T_F')]
    r['T_O'] = [get_param(20, 'T_O')]
    c['T_A'] = get_param(20, 'T_A')  # Air temperature is constant
    r['C_A'] = [0]
    r['C_B'] = [get_param(1000, 'C_B')]  # todo, kg per m3
    r['adjusted_r'] = [0]
    r['r'] = [0]
    r['Q_ER'] = [0]
    c['Q_SR'] = 0  # Todo: Power from stirring is constant
    r['Q_RF'] = [calc_Q_RF(r, c)]
    r['Q_FO'] = [calc_Q_FO(r, c)]
    r['Q_FC'] = [0]
    r['Q_OA'] = [0]


def dT_Rdt(r, c):
    return (r['Q_ER'][-1] + c['Q_SR'] - r['Q_RF'][-1]) / (c['m_R'] * c['c_R'])


def dT_Fdt(r, c):
    return (r['Q_RF'][-1] - r['Q_FO'][-1] - r['Q_FC'][-1]) / (c['m_F'] * c['c_F'])


def dT_Odt(r, c):
    return (r['Q_FO'][-1] - r['Q_OA'][-1]) / (c['m_O'] * c['c_O'])


def calc_r(r, c, C_A, C_B):
    # r is in moles per m3
    return (-c['B'] * C_A * C_B *
            np.e**(-c['E_a'] / (c['R']*(r['T_R'][-1] + 273.15))))
    

def calc_Q_RF(r, c):
    return c['U'] * (r['T_R'][-1] - r['T_F'][-1])


def calc_Q_FO(r, c):
    return c['h_F'] * c['A_Oi'] * (r['T_F'][-1] - r['T_O'][-1])


def calc_Q_OA(r, c):
    radiation = (c['sbc'] * c['emissivity'] * c['A_Oo']
                 * ((r['T_O'][-1] + 273.15)**2 + (c['T_A'] + 273.15)**2)
                 * (r['T_O'][-1] - c['T_A']) * (546.3 + r['T_O'][-1] + c['T_A']))

    convection = 50 * c['A_Oo'] * (r['T_O'][-1] - c['T_A'])

    return radiation + convection


def write_results(r, c):

    # Save all constant values and values for t=0 to constants.txt
    with open("constants.txt", 'w') as f:
    
        f.write("Constants\n")
        for k, v in sorted(c.items()):
            f.write(f'{k:20} {v}\n')
    
        f.write("\nInitial values\n")
        for k, v in sorted(r.items()):
            f.write(f'{k:20} {v[0]}\n')
            
    # Save all data at resolution t=t_interval to "data.csv"
    with open("data.csv", 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(r.keys())
        for i in range(len(r['t'])):
            w.writerow(r[k][i] for k in r.keys())
            
    # Save all data at resolution t=t_interval * 10 to "graphing.csv"
    with open("graphing.csv", 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(r.keys())
        for i in range(len(r['t'])):
            if i % 10 == 0:
                w.writerow(r[k][i] for k in r.keys())
        
        
if __name__ == '__main__':
    main(sys.argv[1:])
