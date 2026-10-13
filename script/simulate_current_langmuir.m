% Monte Carlo Langmuir simulation for the CURRENT double-ellipse experiment.
% Run this script manually. Requires Statistics and Machine Learning Toolbox
% for exact binomial draws. Python and MATLAB seeds do not produce identical
% records, but implement the same distribution and update rule.
clear; clc;
nRecords = 10000;
logCMin = -10.5; logCMax = -8;
dt = 60; timeStop = 7200;
kon = 55800; koff = 1.79e-5;
diameterNm = 6.52;
nSites = 5852; % Current geometry: plane 2654 + top 738 + side 2460.
% Recompute sites if geometry or diameter changes; Python contains geometry code.
rng(7, 'twister');
package = fileparts(fileparts(mfilename('fullpath')));
output = fullfile(package, 'work', 'binding_matlab');
mkdir(output);
time = 0:dt:timeStop;
concentration = 10.^(logCMin + (logCMax-logCMin)*rand(nRecords,1));
pOn = min(kon*concentration*dt, 1);
pOff = min(koff*dt, 1);
occupied = zeros(nRecords,1);
thickness = zeros(nRecords,numel(time),'single');
for j = 2:numel(time)
    % Both events use occupancy at the START of this time step.
    bind = binornd(nSites-occupied, pOn);
    unbind = binornd(occupied, pOff);
    occupied = occupied + bind - unbind;
    thickness(:,j) = (pi*diameterNm/6)*(occupied/nSites);
end
metadata = table((1:nRecords)', (1:nRecords)', ones(nRecords,1), concentration, ...
    'VariableNames', {'recordId','concentrationIndex','monteCarloIndex','concentrationM'});
names = cellstr(compose('thicknessNm_t%04ds',time));
writetable([metadata array2table(thickness,'VariableNames',names)], ...
    fullfile(output, sprintf('protein_total_thickness_0_to_%04ds.csv',timeStop)));
metadata.finalThicknessNm = thickness(:,end);
metadata.finalOccupiedSites = occupied;
writetable(metadata, fullfile(output,'protein_langmuir_records_metadata.csv'));
summary = table(10^logCMin,10^logCMax,dt,kon,koff,nSites,diameterNm, ...
    'VariableNames',{'concentrationMinM','concentrationMaxM','dtSeconds', ...
    'kon_Minv_sinv','koff_sinv','nBindingSites','proteinDiameterNm'});
writetable(summary,fullfile(output,'protein_langmuir_dataset_summary.csv'));
fprintf('Saved %d traces, %d time points to %s\n',nRecords,numel(time),output);
