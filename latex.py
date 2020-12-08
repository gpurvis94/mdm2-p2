import csv

with open("graphing.csv", 'r') as f:
    reader = csv.reader(f)
    r_iter = iter(reader)
    header = next(r_iter)

    data = {k: [] for k in header}
    for i, row in enumerate(r_iter):
        if i % 10 != 0:
            continue
        for k, value in zip(header, row):
            if k == 't':
                data[k].append(round(float(value)/60, 3))
            else:
                data[k].append(round(float(value), 3))

# Specify which rows of data you want to write
keys = ['t', 'T_R', 'T_F', 'C_A', 'C_B', 'r', 'Q_ER', 'Q_RF', 'Q_FC', 'Q_FO', 'Q_OA']

# Write data into format for latex graph
with open("latex.dat", 'w', newline='') as f:
    # latex requires header be fieldnames
    f.write(' '.join(keys) + '\n')
    # space separated values for each field name
    for i in range(len(data[keys[0]])):
        f.write(' '.join(str(data[k][i]) for k in keys) + ' \n')
