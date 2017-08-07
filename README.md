# fort36_analyzer
Does several analysis steps on the data file fort.36 (a product of a successful run of MCCE)
### Usage
The main.py file uses the ArgumentParses module in python, meaning that arguments must be given within the command line itself. Here is what the help command returns:
```
usage: main.py [-h] -f FORT36_LOCATION -out OUTPUT_LOCATION
               [-stdev STANDARD_DEVIATION]

Analyze values in fort36, calculates mean of energies at each pH by default.
Can process more values by if required by user

optional arguments:
  -h, --help            show this help message and exit

Required arguments
                please enter absolute paths where possible:
  -f FORT36_LOCATION, --fort36_location FORT36_LOCATION
                        Location of fort36 file to be analyzed
  -out OUTPUT_LOCATION, --output_location OUTPUT_LOCATION
                        Location of FOLDER containing output csv files

Optional arguments:
  -stdev STANDARD_DEVIATION, --standard_deviation STANDARD_DEVIATION
                        Flag indicating whether stdev should be calculated.
                        False by default
```

### What it does
All the data for analysis lies within the fort.36 file
1. Runs the pandas module to create a dataframe that only consists of the Ave. Energies at each pH (their are multiple such occurences. Finds mean of all the Ave. Energy values for each pH (ranging from 0 through 14). Also, includes a column with standard deviation for the same data (all Ave. Energies at each pH)
2. For each conformer listed in the fort.36 file, finds standard deivation of all occupancy values at each pH. A gigantic dataframe is created.
3. Both dataframes (one in step 1 and one in step 2) are written into files and saved in a directory defined by the user as one of the arguments.

### Example
The command given to produce the `analysis_protein` folder in this repository...
```
python main.py -f fort.36_protein -out analysis_protein
```
