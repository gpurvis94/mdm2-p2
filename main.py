import sys
import numpy as np
import csv

from graphing import multi_plot


def main(args):

    # Iterate the model
    r, c = iterate(args)

    # graphing
    write_initial_conditions(r, c)
    graph_results(r, c)
    write_results(r, c)


def iterate(args):

    # Read command line parameters for later use
    cmd_line_param = {}
    if len(args) > 0:

        # Exit a variable name was passed without a value
        if len(args) % 2 != 0:
            print("User defined args must be followed by a numerical value, e.g:\n"
                  "$ python main.py C_B 100 C_Afeed 1000")
            exit()

        # Store command params in a dictionary to parse later
        args_iter = iter(args)
        for k in args_iter:
            cmd_line_param[k] = float(next(args_iter))

    # x axis time values
    dt = 0.01
    t_min = 0    # This is inclusive
    t_feed_max = 600
    t_max = 900  # Not inclusive
    t = np.arange(t_min, t_max + dt, dt)

    # Containers for constant values and results lists.
    # c maps a 'str' : value, r maps a 'str' : [list of results]
    c = {}
    r = {'t': t}
    initialise_variables(r, c, t_feed_max, cmd_line_param)

    # Skip the first value in t as values for 0 are defined
    t_iterator = iter(t)
    next(t_iterator)

    # Enter the iterative loop
    pump_active = False
    time_pump_active = 0
    total_time_pump_active = 0
    for time in t_iterator:

        # Each iteration represents a step from t1 to t2 where dt = t2 - t1

        # Calculate temperatures at t2 using values from t1: T(t=t2) = T(t=t1) + dTdt * dt
        r['T_R'].append(r['T_R'][-1] + dt * dT_Rdt(r, c))
        r['T_F'].append(r['T_F'][-1] + dt * dT_Fdt(r, c))
        r['T_O'].append(r['T_O'][-1] + dt * dT_Odt(r, c))

        # Update properties dependent upon reactants area
        r['V_R'].append(r['V_R'][-1] + dt * dV_Rdt(r, c, time, t_feed_max))
        r['m_R'].append(c['rho_R'] * r['V_R'][-1])
        r['A_Ii'].append((2 * r['V_R'][-1]) / c['r_Ii'])
        r['U'].append(calc_U(r, c))

        # Mass of reactants in tank after flow in/out, & mass of reactants consumed
        F_Rin = 0 if time > t_feed_max else c['F_Rin']
        F_Rout = 0 if time > t_feed_max else c['F_Rout']
        m_A = r['C_A'][-1]*r['V_R'][-2] + dt*(F_Rin*c['C_Afeed'] - F_Rout*r['C_A'][-1])
        m_B = r['C_B'][-1]*r['V_R'][-2] + dt*(F_Rin*c['C_Bfeed'] - F_Rout*r['C_B'][-1])
        m_r = min(dt * -r['r'][-1] * r['V_R'][-2], m_A, m_B)

        # Concentrations of reactants after mass & volume updates
        r['C_A'].append((m_A - m_r)/r['V_R'][-1])
        r['C_B'].append((m_B - m_r)/r['V_R'][-1])

        # Calc the rate of reaction at t2 after t1->t2 materials added
        r['r'].append(calc_r(r, c))
        r['adjusted_r'].append(-m_r/(r['V_R'][-1] * dt))

        # Now we calculate the power exchange balances for t2 from t2 values

        # Power from Q_ER is the amount of r calculated above
        r['Q_ER'].append(300000 * m_r/dt)
        r['Q_RF'].append(calc_Q_RF(r, c))
        r['Q_FO'].append(calc_Q_FO(r, c))
        r['Q_OA'].append(calc_Q_OA(r, c))

        # Calculate heat loss depending on whether or not the pump is running
        if pump_active:
            # If the pump is running check how long for
            if time_pump_active > c['min_pump_uptime']:
                # Stop the pump if temperature is within threshold and min_pump_uptime met
                if r['T_R'][-1] < c['T_thresh']:
                    pump_active = False
                    time_pump_active = 0
                else:
                    time_pump_active += dt
                    c['total_time_pump_active'] += dt
            else:
                # If the pump hasn't been on long enough, continue running
                time_pump_active += dt
                c['total_time_pump_active'] += dt
        else:
            if r['T_R'][-1] > c['T_thresh']:
                pump_active = True
                time_pump_active = dt
                c['total_time_pump_active'] += dt

        # If the pump is active then calculate energy loss to coolant flow
        r['Q_FC'].append(calc_Q_FC(r, c, pump_active))

    return r, c


def initialise_variables(r, c, t_feed_max, param):

    # Checks if user defined custom value, returns "default" if no "var_name" in cmd line params
    get_user_value = lambda default, var_name: default if var_name not in param else param[var_name]

    # Universal constants ds[1] if 2 in ds else 5
    c['R'] = 8.3144626  # Ideal gas constant
    c['sbc'] = 5.6703744 * 10**-8  # stefan boltzmann constant
    
    # Reactor Dimensions (width, height, area, volume)
    # ALL  volumes and flow rates in m3 (= L / 1000)

    # Tank volume is NOT the volume of reactants due to air padding
    c['V_T'] = 10 / 1000
    c['L_J'] = get_user_value(0.26, 'L_J')                      # Tank (jacket) height
    c['width_J_wall'] = get_user_value(0.001, 'width_J_wall')   # Jacket wall thickness
    c['width_J_fluid'] = get_user_value(0.01, 'width_J_fluid')  # Jacket fluid thickness
    
    # Inside radius: distance to the INSIDE of the inner jacket wall (m)
    c['r_Ii'] = np.sqrt(c['V_T'] / (np.pi * c['L_J']))
    
    # The volume of the tank is defined as 10L. Given that the reactants
    # have a padding of air, the volume of reactants is assumed 95% 10L
    r['V_R'] = [0.5 * c['V_T']]
    # L_R is not constant but there's no point tracking results for it
    c['L_R'] = [r['V_R'][-1] / (np.pi * c['r_Ii']**2)]
    
    # These areas are used for heat transfer so the area in contact
    # with fluid is required, so reactant height is used for inside
    r['A_Ii'] = [2 * np.pi * (c['r_Ii']) * (c['L_R'][-1])]
    
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
    c['h_R'] = get_user_value(2000, 'h_R')  # todo
    c['k_J'] = get_user_value(50, 'k_J')  # todo
    c['h_F'] = get_user_value(2000, 'h_F')  # todo

    # Overall heat transfer coefficient
    r['U'] = [calc_U(r, c)]

    # Specific heats
    c['c_R'] = get_user_value(4184, 'c_R')
    c['c_F'] = get_user_value(4184, 'c_F')
    c['c_O'] = get_user_value(420, 'c_O')
    
    # Densities of reactants = density fluid = density water, units are kg per m3
    c['rho_R'] = get_user_value(1000, 'rho_R')
    c['rho_F'] = get_user_value(1000, 'rho_F')
    c['rho_J'] = get_user_value(8050, 'rho_J')  # Density of jacket wall material (steel)
    
    # Masses = density * volume
    r['m_R'] = [c['rho_R'] * r['V_R'][-1]]
    c['m_F'] = c['rho_F'] * c['V_F']
    c['m_O'] = c['rho_J'] * c['V_O']  # only outer jacket volume required

    # Steel properties for emissivity
    c['emissivity'] = get_user_value(0.5, 'emissivity')  # todo
    
    # Flow rates and concentrations

    # Concentrations of A in feed (F_R = feed rate in = out of tank)
    # todo: for now the feed rate is assumed to "fill" the tank over t_max
    #c['F_Rin'] = 10*((0.005) / 600)
    c['F_Rin'] = get_user_value((c['V_T'] - r['V_R'][-1]) / t_feed_max, 'F_Rin')
    c['F_Rout'] = 0
    c['C_Afeed'] = get_user_value(1000, 'C_Afeed')  # todo, kg per m3
    c['C_Bfeed'] = get_user_value(0, 'C_Bfeed')  # todo, kg per m3

    # Temperature & flow rate of coolant (approximated as 10% jacket volume/s)
    c['T_C'] = get_user_value(10, 'T_C')
    c['F_C'] = 0.1 * c['V_F']

    # Pump temperature thresholds and minimum run time
    c['T_thresh'] = get_user_value(35, 'T_thresh')
    c['min_pump_uptime'] = get_user_value(10, 'min_pump_uptime')
    
    # Rates of reaction constants
    #c['B'] = 0.7 * 10**4  # USE THIS VALUE FOR INITIAL TEST
    c['B'] = 1 * 10**5  # USE THIS VALUE FOR VARIABLE VARIATION
    c['E_a'] = 51217  # calculated for a value of 300 kelvins
    
    # Initial values for iteration
    
    # results string maps to results array
    r['T_R'] = [get_user_value(20, 'T_R')]
    r['T_F'] = [get_user_value(10, 'T_F')]
    r['T_O'] = [get_user_value(20, 'T_O')]
    c['T_A'] = get_user_value(20, 'T_A')  # Air temperature is constant
    r['C_A'] = [0]
    r['C_B'] = [get_user_value(1000, 'C_B')]  # todo, kg per m3
    r['adjusted_r'] = [0]
    r['r'] = [0]
    r['Q_ER'] = [0]
    c['Q_SR'] = 0  # Todo: Power from stirring is constant
    r['Q_RF'] = [calc_Q_RF(r, c)]
    r['Q_FO'] = [calc_Q_FO(r, c)]
    r['Q_FC'] = [0]
    r['Q_OA'] = [0]

    # Extra stuff
    c['total_time_pump_active'] = 0


def dT_Rdt(r, c):
    return (r['Q_ER'][-1] + c['Q_SR'] - r['Q_RF'][-1]) / (r['m_R'][-1] * c['c_R'])


def dT_Fdt(r, c):
    return (r['Q_RF'][-1] - r['Q_FO'][-1] - r['Q_FC'][-1]) / (c['m_F'] * c['c_F'])


def dT_Odt(r, c):
    return (r['Q_FO'][-1] - r['Q_OA'][-1]) / (c['m_O'] * c['c_O'])


def dV_Rdt(r, c, time, t_feed_max):
    return get_F_Rin(r, c, time, t_feed_max) - get_F_Rout(r, c, time, t_feed_max)


def get_F_Rin(r, c, time, t_feed_max):
    return 0 if time > t_feed_max else c['F_Rin']


def get_F_Rout(r, c, time, t_feed_max):
    return 0 if time > t_feed_max else c['F_Rout']


def calc_U(r, c):
    return 1 / ((1 / (c['h_R'] * r['A_Ii'][-1]))
                   + (np.log(c['r_Io']/c['r_Ii']) / (2 * np.pi * c['L_J'] * c['k_J']))
                   + (1 / (c['h_F'] * c['A_Io'])))


def calc_r(r, c):
    return (-c['B'] * r['C_A'][-1] * r['C_B'][-1] *
            np.e**(-c['E_a'] / (c['R']*(r['T_R'][-1] + 273.15))))
    

def calc_Q_RF(r, c):
    return r['U'][-1] * (r['T_R'][-1] - r['T_F'][-1])


def calc_Q_FO(r, c):
    return c['h_F'] * c['A_Oi'] * (r['T_F'][-1] - r['T_O'][-1])


def calc_Q_FC(r, c, pump_active):
    if pump_active:
        return c['m_F'] * c['c_F'] * (c['F_C']/c['V_F']) * (r['T_F'][-1] - c['T_C'])
    else:
        return 0


def calc_Q_OA(r, c):
    radiation = (c['sbc'] * c['emissivity'] * c['A_Oo']
                 * ((r['T_O'][-1] + 273.15)**2 + (c['T_A'] + 273.15)**2)
                 * (r['T_O'][-1] - c['T_A']) * (546.3 + r['T_O'][-1] + c['T_A']))

    convection = 50 * c['A_Oo'] * (r['T_O'][-1] - c['T_A'])

    return radiation + convection


def graph_results(r, c):

    multi_plot(r)


def write_initial_conditions(r, c):
    # Save all constant values and values for t=0 to constants.txt
    with open("constants.txt", 'w') as f:

        f.write("Constants\n")
        for k, v in sorted(c.items()):
            f.write(f'{k:20} {v}\n')

        f.write("\nInitial values\n")
        for k, v in sorted(r.items()):
            f.write(f'{k:20} {v[0]}\n')


def write_results(r, c):
    # Save all data at resolution t=dt * 10 to "graphing.csv"
    with open("graphing.csv", 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(r.keys())
        for i in range(len(r['t'])):
            if i % 10 == 0:
                w.writerow(r[k][i] for k in r.keys())
        
        
if __name__ == '__main__':
    main(sys.argv[1:])
