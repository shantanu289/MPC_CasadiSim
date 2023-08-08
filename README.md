# MPC_CasadiSim
MATLAB-based MPC execution using Casdadi in co-sim

Necessary GIT repos - 
1. Casadi (MATLAB) v3.5.5 : https://web.casadi.org/get/
2. Clothoid fitting (Enrico Bertolazzi et. al) : https://github.com/ebertolazzi/G1fitting
Reference GIT repos -
1. MPC MHE Workshop (Md. Mehrez et. al): https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi/tree/master/workshop_github
2. TIAGo hamr navigation (for obstacle avoidance examples in MATLAB simulations - kinematic model) : https://github.com/rejzi/TIAGo-hamr-navigation/tree/master

Folders at the same level - 
1. Casadi Bolcks --> MATLAB system blocks invoking casadi-based optimization (solver : ipopt by default)
2. Simulink models --> Casadi_Custom Driver and standalone simulink models for each maneuver (for testing purpose)
3. Casadi repo
4. Clothoid fitting repo
5. Test codegen zip --> To test casadi blog example with different solvers for code genereation
6. MPC Run results --> for all different maneuvers 
7. Road data for reference trajectories (build-up purpose)
8. Trial Code Blocks --> MATLAB simulations (early steps) 
