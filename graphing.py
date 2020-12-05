import sys
import csv
import numpy as np
import matplotlib.pyplot as plt


def main():
    with open("graphing.csv", 'r') as f:
        reader = csv.reader(f)
        r_iter = iter(reader)
        header = next(r_iter)
        
        data = {k: [] for k in header}
        for i, row in enumerate(r_iter):
            if i % 10 != 0:
                continue
            for k, value in zip(header, row):
                data[k].append(float(value))     

    for k in data.keys():
        print(f'{k:15} min {min(data[k]):15.6f} max {max(data[k]):15.6f}')   

    multi_plot(data)


def multi_plot(data):
    fig, axs = plt.subplots(2, 2)

    # x axis
    x = [d/60 for d in data['t']]
    y = [np.sin(x ** 2) for x in x] # placeholder


    # Temperature
    T_R = data['T_R']
    T_F = data['T_F']
    T_O = data['T_O']

    # Plot the data
    axs[0, 0].plot(x, T_R, color='red', label='T_R')
    axs[0, 0].plot(x, T_O, color='green', label='T_O')
    axs[0, 0].plot(x, T_F, color='blue', label='T_F')

    axs[0, 0].set_ylim(bottom=0)


    # Concentration / Rate

    # Left y axis
    C_A = data['C_A']
    C_B = data['C_B']

    axs[0, 1].plot(x, C_A, color='red', label='C_A')
    axs[0, 1].plot(x, C_B, color='blue', label='C_B')

    axs[0, 1].set_ylim(bottom=0)

    # right y axis
    r = data['r']
    adj_r = data['adjusted_r']

    ax01_2 = axs[0, 1].twinx()
    ax01_2.plot(x, r, color='black', label='r')
    ax01_2.plot(x, adj_r, color='grey', label='r_adj')

    ax01_2.set_ylim(top=0)
    ax01_2.legend(loc='upper center', frameon=False, fontsize='x-small')


    # Power

    powers = ['Q_ER', 'Q_RF', 'Q_FO', 'Q_FC', 'Q_OA']
    y_powers = [data[s] for s in powers]

    for y, s in zip(y_powers, powers):
        axs[1, 0].plot(x, y, label=s)

    #y axis needs to be logarithmic
    axs[1, 0].set_yscale('log')
    axs[1, 0].set_ylim((0, 10**5))


    # Last plot unused for now

    var = data['V_R']
    #var2 = data['var2']
    axs[1, 1].plot(x, var, color='red', label='V_R')
    #axs[1, 1].plot(x, var2, color='red', label='var')


    # Set xlim and format axis
    for i in axs:
        for j in i:
            j.set_xlim((min(x), max(x)))
            j.legend(frameon=False, fontsize='x-small', loc='best')

    plt.show()

    
if __name__ == '__main__':
    main()
