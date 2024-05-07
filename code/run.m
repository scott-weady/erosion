%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run script for phase-field erosion simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Re = 100; %Reynolds number
sigma = 100; %erodibility ratio

% Domain size
Lx = 8; %length
Ly = 4; %width
Lz = 4; %height

Nz = 64; %vertical resolution

C_clay = 0.01; %clay erodibility
C_inc = C_clay/sigma; %inclusion erodibility

tf = 1000; %final time
tplt = 1; %plotting time (set to inf to not plot)
dt = 0.01; %time step
tsave = 5*tplt; %save time

main %run