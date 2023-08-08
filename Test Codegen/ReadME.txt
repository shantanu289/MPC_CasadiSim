1. Download the casadi package for MATLAB 2016 or later from - https://web.casadi.org/get/ --> already present in this zip
2. Download 'casadi_block' and 'mpc_demo' from - https://web.casadi.org/blog/mpc-simulink/ --> save with a different name as updated codes are present in the zip
3. Put in the same folder as the Simulink model and casadi_block 
4. Open the Simulink model and casadi_block 
5. Change the continuous integrator to discrete integrator
6. Ensure the sample times for all ZOH and Noise blocks are the same 
5. in casadi_block : check for lines 147-153
   for 'ipopt' solver --> line 153 : solver = solver = casadi.nlpsol('solver', 'ipopt', prob, opts);
   for 'qpoases' solver --> line 153 : solver = casadi.qpsol('solver', 'qpoases', prob); 
6. Try code generation 
7. Error : "addpath function not supported for code generation" --> remove the add path statement from the code and paste and run it in the terminal 
8. Next Error : "import statements not supported for codegen" --> change all the functions from casadi to 'casadi.Function' 
   Note --> steps 7 and 8 are already done in the code in present in this zip folder.
9. Try code generation with the new code --> Error shown in the 'error.png' flie in this zip 