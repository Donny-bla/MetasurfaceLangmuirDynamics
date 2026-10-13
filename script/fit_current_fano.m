function fit_current_fano()
% Refit the correct raw 0 nm and 1-6 nm spectra. Requires Curve Fitting Toolbox.
% This updates the bundled fit tables consumed by python/01_fit_fano_forward.py.
here = fileparts(mfilename('fullpath'));
addpath(here);
package = fileparts(here);
fitOne(fullfile(package,'inputs','raw','DoubleEllipis_S843_theta205.txt'), ...
    fullfile(package,'inputs','fit_zero'),0,200,here);
fitOne(fullfile(package,'inputs','raw','DoubleEllipis_S843_theta205_thickness1-6nm.txt'), ...
    fullfile(package,'inputs','fit_coated'),1:6,100,here);
end

function fitOne(dataFile,outputFolder,thicknessNm,pointsPerSpectrum,here)
makePlots = true;
% Run in a function workspace because the fitting script clears local variables.
run(fullfile(here,'analyzeFanoThicknessSpectra.m'));
end
