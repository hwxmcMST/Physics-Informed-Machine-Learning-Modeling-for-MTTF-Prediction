clc;clear
smodel = createpde("structural","static-solid");
importGeometry(smodel,"Blade.stl");
msh = generateMesh(smodel,"Hmax",0.01);
E = 227E9; % in Pa
CTE = 12.7E-6; % in 1/K
nu = 0.27; 
structuralProperties(smodel,"YoungsModulus",E, ...
                            "PoissonsRatio",nu, ...
                            "CTE",CTE);
structuralBC(smodel,"Face",3,"Constraint","fixed");
p1 = 5e5; %in Pa
p2 = 5e5; %in Pa

structuralBoundaryLoad(smodel,"Face",11,"Pressure",p1); % Pressure side
structuralBoundaryLoad(smodel,"Face",10,"Pressure",p2);  % Suction side 
Rs = solve(smodel);
find(Rs.VonMisesStress==max(Rs.VonMisesStress))
find(Rs.Displacement.uz==max(Rs.Displacement.uz))
max(Rs.Displacement.uz)
max(Rs.VonMisesStress)