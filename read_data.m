% Read Ksn data Green's functions and lithologic units

clear all

if 0
    data = load('ksnData_empiricalGreensFCNs.csv');
elseif 0
    load('ksnData_empiricalGreensFCNs_preAGUtest.mat');
    data = greens_ksn;
elseif 0
    load('ksnData_empiricalGreensFCNs_preAGUtest_03_coarse.mat');
    data = greens_ksn;
elseif 1
    load ksnData_iBemUzGreens_rstrain_cl1
    data = greens_ksn;
end

%%

x = data(:,1);
y = data(:,2);
ksn = data(:,5);
geo = data(:,6);
G = data(:,[7 8 9]);

if 0
    save ksnData_ibem_cl1.mat x y ksn geo G
end

%% Test plotting

figure
scatter(x,y,[],ksn,'filled')
colorbar
caxis([min(ksn) 300])


figure
scatter(x,y,[],geo,'filled')
colormap(hsv(20))
colorbar
% caxis([min(ksn) 300])

% 
% -- 
% Andreas Petros Mavrommatis
% PhD Candidate, Geophysics
% Stanford University
% https://sites.google.com/site/amavrommatis/