import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

def main():
    with open("graphing.csv", 'r') as f:
        reader = csv.reader(f)
        r_iter = iter(reader)
        header = next(r_iter)
        
        data = { k : [] for k in header}
        for i, row in enumerate(r_iter):
            if i % 10 != 0:
                continue
            for k, value in zip(header, row):
                data[k].append(float(value))     

    for k in data.keys():
        print(f'{k:15} min {min(data[k]):15.6f} max {max(data[k]):15.6f}')   
        
    plot_q(data)
    plot_rate(data)
    plot_temp(data)
    
def plot_q(data):
    # Prepare the data
    x = [d/60 for d in data['t']]
    
    # Plot the normalized values
    
    plt.plot(x, data['Q_ER'], label='Q_ER')
    #plt.plot(x, data['Q_RF'], label='Q_RF')
    #plt.plot(x, data['Q_FO'], label='Q_FO')
    #plt.plot(x, data['Q_OA'], label='Q_OA')
    #plt.plot(x, data['Q_FC'], label='Q_FC')
    
    plt.xlim((0, 30))
    plt.xticks(np.arange(0, 30, 1))

    plt.legend()

    plt.show()

def plot_temp(data):
    # Prepare the data
    x = [d/60 for d in data['t']]
    y1 = data['T_R']
    y2 = data['T_F']
    y3 = data['T_O']

    # Plot the data
    plt.plot(x, y1, label='T_R')
    plt.plot(x, y2, label='T_F')
    plt.plot(x, y3, label='T_O')
    plt.xlim((0, 30))
    plt.xticks(np.arange(0, 30, 2))
    plt.ylim((0, 30))
    plt.yticks(np.arange(0, 70,5))

    plt.legend()

    plt.show()
    
def plot_rate(data):
    # Prepare the data
    x = [d/60 for d in data['t']]
    y1 = data['C_A']
    y2 = data['C_B']
    
    fig, ax = plt.subplots()
    ax.plot(x, y1, color='red', label='C_A')
    ax.plot(x, y2, color='blue', label='C_B')
    
    y3 = data['r']
    y4 = data['adjusted_r']
    
    
    ax2 = ax.twinx()
    ax2.plot(x, y3, color='green', label='r')
    ax2.plot(x, y4, color='yellow', label='adj_r')
    
    plt.xlim((0, 30))
    plt.xticks(np.arange(0, 31, 1))

    ax.legend()
    ax2.legend()

    plt.show()
    
    
if __name__ == '__main__':
    main()