# Vapor-Liquid-Equilibrium-VLE-Calculations-for-Multicomponent-Mixtures-Using-MATLAB
This project involves creating a MATLAB program for Vapor-Liquid Equilibrium (VLE) calculations of multicomponent mixtures. It predicts phase compositions using Raoult's Law and the Antoine equation, providing key insights for designing and optimizing separation units in chemical engineering.

# Description
This MATLAB script focuses on solving several thermodynamic problems related to bubble point and dew point calculations using iterative numerical methods, specifically Newton's method. The script includes:

# 1. Initialization:
The script starts by clearing the workspace and defining key constants, including pressure, liquid compositions, and parameters for the Antoine equation, which are crucial for calculating the saturated vapor pressures and temperatures of the chemical components.
# 2. Problem 2a:
Bubble Point Calculation:
The script calculates the bubble point temperature and the corresponding vapor composition for a given liquid mixture. This involves estimating the temperature iteratively using Newton's method until convergence is achieved, based on the Antoine equation's parameters.
# 3. Problem 2b:
Dew Point Calculation:
Here, the script calculates the dew point temperature for a given vapor mixture. Similar to the bubble point, it uses Newton's method to iteratively estimate the temperature and determine the liquid composition.
# 4. Problem 3a:
Advanced Bubble Point Calculation:
This section extends the bubble point calculation by incorporating non-idealities using the UNIFAC model to compute activity coefficients and the virial equation to account for non-ideal gas behavior. The process involves multiple nested iterations to account for the interactions between the different species in the mixture.
# 5. Problem 3b:
Advanced Dew Point Calculation:
The script performs a more complex dew point calculation, similar to Problem 3a, but for the dew point instead of the bubble point. It iterates over temperature and composition until the conditions for the dew point are satisfied, considering non-ideal interactions.
# 6. Problem 4:
Pressure Calculations:
This problem involves calculating the dew and bubble pressures for a given temperature and composition, followed by checking if the mixture is within the two-phase region. If so, it computes the vapor and liquid phase compositions and their respective mole fractions.
# 7. Supporting Functions:
The script defines several custom functions (CalcTsats, CalcPsats, CalcdPsats, TbubbleNewton, TdewNewton, UNIFAC, etc.) to handle specific calculations like determining saturated pressures, applying Newton's method, and calculating activity coefficients using the UNIFAC model.
The script prints the results of each calculation, including the bubble and dew temperatures, vapor and liquid compositions, and phase mole fractions, helping solve complex thermodynamic problems for multi-component systems.
