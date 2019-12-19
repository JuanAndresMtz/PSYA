# PSYA
Polar Space Yield Analysis

This repository is for the Elementary Mode Analysis approach presented in Martinez et al. (Submited 2019)
This repository contains the PSYA code on MATLAB for the calculation of the Polar Solution Space for a given EM distribution. On this work we used efmtool (Terzer and Stelling,) to compute the elementary flux modes (EFMs) from the metabolic network (MetabolicNetwork.txt). The EM matrix output was then reduced by constraints by the EMConstraints.m program (available on this repository) as described in Martinez et al.(Submited 2019).

For upload space reasons on this repositores Matrices needed for runing PSYA.mat are available at: https://drive.google.com/drive/folders/1hdTPjNkeGHQT-1gBfgc52-vyDucY3IRs?usp=sharing

The google drive repository contains 4 initial matrices (EM_Matrix,EM_Matrix_4, EM_Matrix_6 and EM_Matrix_8) which correspond to the efmtools output and the reduced EM matrices outputs by 4, 6 and 8 constraints as described in Martinez et al.(Submited 2019).

The google drive repository also contains a matrix of experimental data of a CHO-S fermentation (Exp_Data), which is conformed by columns of Time, Biomass, Glucose, Glutamine, Lactate and Glutamate data.

The google drive repository also contains an array of modeled biomass rates on 0.1 h time steps, calculated from extracellular mathematical modeling approaches as described in Martinez et al.(Submited 2019).

This repository last update was performed on: Dec 18 2009.

For more information please contact: andres.amdg@gmail.com

efttool can be found at: https://csb.ethz.ch/tools/software/efmtool.html

 ***************************************************************************************************************************
 **** Polar Space Yield Analysis (PSYA) program solution example ***********************************************************

 ***************************************************************************************************************************
 
                                                                              by J.A.Martinez, L.A.Palomares and O.T.Ramirez
 
                                                                                               contact:andres.amdg@gmail.com 
 
 ************************************************************************************************************************** 
 
THIS IS AN EXAMPLE CODE ON MATLAB FOR THE USE OF THE PSYA APPROACH. IT WAS MADE FOR SOLVING THE SUBMITED PRESENTED DATA.  
IT IS IMPLEMENTED FOR CCM NETWORK CHO-S CELL ELEMENTARY MODE SELECTION. ITS USE IS HEREIN AUTHORIZED WITH OR WITHOUT 
MODIFICATIONS EITHER ON THE SAME NETWORK, OTHER SYSTEMS OR OTHERS INTENDED METABOLIC MODELING NEEDS, 
PLEASE CITE THE IF TOTAL OR FRACTIONAL SECTIONS USED OR MODIFIED. IF YOU NEED MORE DATA PLEASE CONTACT US.

**************************************************************************************************************************
 
%%%%%%%%%% EmSoftwaretool was used for EM calculation                             
EM Matrix Loading from Efmtools software is performed on this script with its output diverted to the following variables:                   

efmMatrix.int_met = names of intracellular metabolites (array)       
efmMatrix.ext_met = names of extracellular metabolites (array)        
efmMatrix.react_name = reaction names (array)                        
efmMatrix.st = stoichiometric matrix (rows correspond to internal metabolites, columns to reactions)                    
efmMatrix.ext = same structure as st, but rows correspond to external metabolites
efmMatrix.ems = Elementary modes data (EMs), (rows correspond to reactions, columns to EM)                            
efmMatrix.ext_out = EMs output matrix (efmMatrix.ext*efmMatrix.ems)    


