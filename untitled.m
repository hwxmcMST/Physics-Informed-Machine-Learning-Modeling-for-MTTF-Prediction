smodel = createpde("structural","static-solid");
importGeometry(smodel,"Blade.stl");
figure
pdegplot(smodel,"FaceLabels","on","FaceAlpha",0.5)
msh = generateMesh(smodel,"Hmax",0.01);
pdemesh(smodel)