# mdm2-p2

## Usage
Firstly, everything will still be a buggy mess, pls lmk if something breaks 
or doesn't work the way you expect it to. Also I haven't included any sanity
checks so if you use it wrong (e.g. passing `T_R not_a_number` as a command 
line parameter) it'll break. 
### main.py
Use this file to generate data. For example:

`$ python main.py`

This runs the program with the constants and initial values defined in the 
code. This generates three files.
 - `constants.txt` lists all initial and constant values used.
 - `data.csv` a csv containing every recorded data point. The default time 
 interval is 0.01 seconds so this file is large.
 - `graphing.csv` lists every tenth data point, a resolution of 0.1 
 seconds, which makes the file less 
 time consuming to parse when reading for graphing purposes.

If you want to run the code with different constant or initial values, 
you can either edit the code or use the command line, like so.

`$ python main.py min_uptime 20 T_A 30`

This sets the air temperature to 30 celsius and minimum cooling pump uptime to
0 seconds. Not all variables can be assigned to this way, e.g. all reactor
dimensions are calculated from `L_J` and the tank volume being 10 litres. 
To see which variables can be assigned to using the command line search 
`get_param` in the code.
 
 ### graphing.py
 Run this script without any command line params.
 
 `$ python grpahing.py`
 
 This file reads `graphing.csv` and generates graphs displaying 
 temperatures, concentrations & rates of reaction, and powers. The 4th 
 subplot is left blank if you want to add anything else.
 
 The minimum and maximum values of each results array is printed to the 
 console.
