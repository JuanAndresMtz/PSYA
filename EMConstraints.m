%************************************************************************
%**** EM Number Reduction by Simple Constraints PSYA Complement  ********
%************************************************************************
%                          by J.A.Martinez, L.A.Palomares and O.T.Ramirez 
%                                           contact:andres.amdg@gmail.com 
%************************************************************************
% THIS IS A COMPLEMENT FOR THE EXAMPLE CODE PSYA.m 
% IT WAS MADE FOR SOLVING THE SUBMITED PRESENTED DATA BY MARTINEZ ET. AL 
% IT IS IMPLEMENTED FOR CCM NETWORK CHO-S CELL ELEMENTARY MODE SELECTION
% ITS USE IS HEREIN AUTHORIZED WITH OR WITHOUT MODIFICATIONS 
% EITHER ON THE SAME NETWORK, OTHER SYSTEMS OR OTHERS INTENDED METABOLIC
% MODELING NEEDS, PLEASE CITE THE IF TOTAL OR FRACTIONAL SECTIONS 
% USED OR MODIFIED. IF YOU NEED MORE DATA PLEASE CONTACT US.
%************************************************************************
%************************************************************************

 
%%%%%%%%%%%%%% Initial Data allocation and normalization %%%%%%%%%%%%%%%%%
% EmSoftwaretool was used for EM calculation
% EM Matrix Loading from Efmtools software is performed on this script 
% with its output diverted to the following variables:
%  efmMatrix.int_met = names of intracellular metabolites (array)
%  efmMatrix.ext_met = names of extracellular metabolites (array)
%  efmMatrix.react_name = reaction names (array)
%  efmMatrix.st = stoichiometric matrix (rows correspond to internal 
%                 metabolites, columns to reactions)
%  efmMatrix.ext = same structure as st, but rows correspond to 
%                  external metabolites
%  efmMatrix.ems = Elementary modes data (EMs), (rows correspond to 
%                  reactions, columns to EM)
%  efmMatrix.ext_out = EMs output matrix (efmMatrix.ext*efmMatrix.ems)    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 4 Constraints EM number reduction
clear variables %clear previous data
close all
clc
load('EM_Matrix.mat') %%% All EMs efmMatrix struct as described previously. 

% al restrictions can be changed, number increased or reduced 
Restriction_1=double((efmMatrix.ext_out(25,:)>0));%BIOMAS constraint
Restriction_2=double((efmMatrix.ext_out(2,:)>0)); %LACTATE constraint
Restriction_3=double((efmMatrix.ext_out(7,:)>0)); %GLUTAMATE constraint
Restriction_4=double((efmMatrix.ext_out(12,:)<0));%GLUTAMINE constraint

%%% Sum previously boleean (1,0) restriction data for all EMs 
ResctrictionSum=Restriction_1+Restriction_2+Restriction_3+Restriction_4;
Restriction_Pos=find(ResctrictionSum==4); %Find positions with all contraints

%%% Find EM data
for l =1:numel(Restriction_Pos)
   Temporary.ems(:,l)=efmMatrix.ems(:,Restriction_Pos(l));
   Temporary.ext_out(:,l)=efmMatrix.ext_out(:,Restriction_Pos(l));
end
%%% Transpose EM data
efmMatrix.ems=Temporary.ems;
efmMatrix.ext_out=Temporary.ext_out;

%clear all not needed variables
clear Temporary Restriction_1 Restriction_2 Restriction_3 Restriction_4 ResctrictionSum Restriction_Pos l

%Save contraint matrix
save('EM_Matrix_4.mat','-v7.3')


%%% 6 Constraints EM number reduction
clear variables %clear previous data
close all
clc
load('EM_Matrix.mat') %%% All EMs efmMatrix struct as described previously. 

% al restrictions can be changed, number increased or reduced 
Restriction_1=double((efmMatrix.ext_out(25,:)>0));%BIOMAS constraint
Restriction_2=double((efmMatrix.ext_out(2,:)>0)); %LACTATE constraint
Restriction_3=double((efmMatrix.ext_out(7,:)>0)); %GLUTAMATE constraint
Restriction_4=double((efmMatrix.ext_out(12,:)<0));%GLUTAMINE constraint
Restriction_5=double((efmMatrix.ems(24,:)>0)); % V39 PYR[c]toPYR[m] constraint
Restriction_6=double((efmMatrix.ems(79,:)>0)); % V53r SUCCOA[m] = SUCC[m] + GTP constraint

%%% Sum previously boleean (1,0) restriction data for all EMs
ResctrictionSum=Restriction_1+Restriction_2+Restriction_3+Restriction_4+Restriction_5+Restriction_6;
Restriction_Pos=find(ResctrictionSum==6); %Find positions with all contraints

%%% Find EM data
for l =1:numel(Restriction_Pos)
   Temporary.ems(:,l)=efmMatrix.ems(:,Restriction_Pos(l));
   Temporary.ext_out(:,l)=efmMatrix.ext_out(:,Restriction_Pos(l));
end
%%% Transpose EM data
efmMatrix.ems=Temporary.ems;
efmMatrix.ext_out=Temporary.ext_out;

%clear all not needed variables
clear Temporary Restriction_1 Restriction_2 Restriction_3 Restriction_4 Restriction_5 Restriction_6 ResctrictionSum Restriction_Pos  l
%Save contraint matrix
save('EM_Matrix_6.mat','-v7.3')


%%% 8 Constraints EM number reduction
clear variables %clear previous data
close all
clc
load('EM_Matrix.mat') %%% All EMs efmMatrix struct as described previously. 

% al restrictions can be changed, number increased or reduced 
Restriction_1=double((efmMatrix.ext_out(25,:)>0)); %BIOMAS constraint
Restriction_2=double((efmMatrix.ext_out(2,:)>0));%sxz2>0 %LACTATE constraint
Restriction_3=double((efmMatrix.ext_out(7,:)>0)); %GLUTAMATE constraint
Restriction_4=double((efmMatrix.ext_out(12,:)<0));%GLUTAMINE constraint
Restriction_5=double((efmMatrix.ems(24,:)>0)); % v39 : PYR[c]toPYR[m] constraint
Restriction_6=double((efmMatrix.ems(79,:)>0)); % v53r : SUCCOA[m] = SUCC[m] + GTP constraint
Restriction_7=double((efmMatrix.ems(26,:)>0)); % v41 : GLN[c] = GLN[m] constraint
Restriction_8=double((efmMatrix.ems(22,:)>0)); % v30 : PEP[c] = PYR[c] + ATP constraint

%%% Sum previously boleean (1,0) restriction data for all EMs 
ResctrictionSum=Restriction_1+Restriction_2+Restriction_3+Restriction_4+Restriction_5+Restriction_6+Restriction_7+Restriction_8;
Restriction_Pos=find(ResctrictionSum==8); %Find positions with all contraints

%%% Find EM data
for l =1:numel(Restriction_Pos)
   Temporary.ems(:,l)=efmMatrix.ems(:,Restriction_Pos(l));
   Temporary.ext_out(:,l)=efmMatrix.ext_out(:,Restriction_Pos(l));
end
%%% Transpose EM data
efmMatrix.ems=Temporary.ems;
efmMatrix.ext_out=Temporary.ext_out;

%clear all not needed variables
clear Temporary Restriction_1 Restriction_2 Restriction_3 Restriction_4 Restriction_5 Restriction_6 Restriction_7 Restriction_8 ResctrictionSum Restriction_Pos  l

%Save contraint matrix
save('EM_Matrix_8.mat','-v7.3')

clear variables %clear previous data
close all
clc
