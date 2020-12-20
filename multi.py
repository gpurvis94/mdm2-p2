import time

from main import iterate


# Progress bar for long iterations
def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd="\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()


def dimensions():

    # Don't want to deal with debugging rounding errors from np linspace
    # with decimal step sizes so just a static definition for domain
    steps = [0.1, 0.2, 0.3, 0.4, 0.5]
    trials = [['L_J', x] for x in steps]
    data = {}

    printProgressBar(0, len(trials), prefix='Progress:', suffix='Complete', length=50)
    for i, trial in enumerate(trials):
        data[trial[1]] = iterate(trial)
        printProgressBar(i + 1, len(trials), prefix='Progress:', suffix='Complete', length=50)

    max_temps = [max(data[step][0]['T_R']) for step in steps]
    total_time_pump_active = [data[step][1]['total_time_pump_active'] for step in steps]

    print("")
    print("Dimensions", steps)
    print(f"Max temperatures")
    print(max_temps)
    print(f"Total time pump active")
    print(total_time_pump_active)
    print("")


def feed_rate():

    # Don't want to deal with debugging rounding errors from np linspace
    # with decimal step sizes so just a static definition for domain
    base = 0.005 / 600
    steps = [base * 1/2, base * 2/3, base, base * 1.5, base * 2]
    trials = [['F_Rin', x] for x in steps]
    data = {}

    printProgressBar(0, len(trials), prefix='Progress:', suffix='Complete', length=50)
    for i, trial in enumerate(trials):
        data[trial[1]] = iterate(trial)
        printProgressBar(i + 1, len(trials), prefix='Progress:', suffix='Complete', length=50)

    max_temps = [max(data[step][0]['T_R']) for step in steps]
    max_rate = [min(data[step][0]['r']) for step in steps]
    avg_rate = [average(data[step][0]['r']) for step in steps]

    print("")
    print("Feed rate", steps)
    print([i * 1000000 for i in steps])
    print(f"Max temperatures")
    print(max_temps)
    print(f"Max rate of reaction")
    print(max_rate)
    print(f"Avg rate of reaction")
    print(avg_rate)
    print("")


def average(lst):
    return sum(lst) / len(lst)


dimensions()

feed_rate()
