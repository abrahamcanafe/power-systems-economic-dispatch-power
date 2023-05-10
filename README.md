# Optimal Economic Dispatch Linear Program
Mixed-integer linear program that simulates the optimal economic dispatch over a 24-hour horizon for a local power system with hydroelectric generators. The final project written 2022 for my Power Systems Economics class at my undergraduate university.

## Built With

* MATLAB
  * YALMIP toolbox

<!-- ABOUT THE PROJECT -->
## Project Details

Information from the power system model diagrammed in `Course Project.xlsx` was used to write a linear program (LP) for optimal economic dispatch.

The following data were considered when writing the LP:
* Hydroelectric unit capacity limits: minimum and maximum power, startup costs, duty factors (cfs/MWh)
* Reservoir elevation limits
* Locational marginal prices over a 24-hour horizon

## Power System Layout
![image](https://github.com/abrahamcanafe/power-systems-optimal-economic-dispatch/blob/main/hydro_power_system.png)


<!-- GETTING STARTED -->
## Sample output
![image](https://github.com/abrahamcanafe/power-systems-optimal-economic-dispatch/blob/main/EEE259_Final_Project_Output.png)


