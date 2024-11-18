function [SpS,El] = myload(filename,species)
load(filename);
isp = myfind({Sp.Name},species);
SpS = Sp(isp);
for i=1:length(El)
    El(i).Mass = El(i).Mass*1e-3;
end