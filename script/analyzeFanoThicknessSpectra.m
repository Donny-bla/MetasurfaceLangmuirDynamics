% Fit Fano resonances for spectra measured at different input thicknesses.
%
% Expected TXT formats:
%   1) Two columns:   wavelength_nm   transmission
%      Rows 1:50 are spectrum 1, rows 51:100 are spectrum 2, etc.
%      In this case set thicknessNm manually below.
%
%   2) Three columns: thickness_nm   wavelength_nm   transmission
%      Thickness is read from the first column of each 50-row block.
%
% Output:
%   fano_thickness_fit_results.csv contains one fitted row per spectrum.
%   fano_parameter_changes.csv contains changes relative to the first spectrum.

clearvars -except dataFile outputFolder thicknessNm pointsPerSpectrum ...
    makePlots fitGridPoints fwhmGridPoints fitWavelengthStartNm ...
    fitWavelengthStopNm; clc;

if ~exist('dataFile', 'var')
    dataFile = fullfile(pwd, 'smalldataset', 'thickness005-060nm.txt');
end
if ~exist('outputFolder', 'var')
    outputFolder = fullfile(pwd, 'fano_thickness_analysis');
end
if ~exist('pointsPerSpectrum', 'var')
    pointsPerSpectrum = 100;
end

% Use this only for two-column input. Leave [] if the TXT has thickness in
% column 1 or if thickness should be inferred from a filename. The number of
% entries must equal the number of spectra.
if ~exist('thicknessNm', 'var')
    thicknessNm = 1:1:6;
end
% Example manual override:
% thicknessNm = [1 2 3 4 5 6];

if ~exist('makePlots', 'var')
    makePlots = true;
end
if ~exist('fitGridPoints', 'var')
    fitGridPoints = 5000;
end
if ~exist('fwhmGridPoints', 'var')
    fwhmGridPoints = 20000;
end
if ~exist('fitWavelengthStartNm', 'var')
    fitWavelengthStartNm = [];
end
if ~exist('fitWavelengthStopNm', 'var')
    fitWavelengthStopNm = [];
end

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

spectra = readThicknessSpectra(dataFile, pointsPerSpectrum, thicknessNm);
nSpectra = numel(spectra);

fitRows = table();

for i = 1:nSpectra
    lambda = spectra(i).lambda;
    transmission = spectra(i).transmission;
    [lambda, transmission] = selectWavelengthWindow(lambda, transmission, ...
        fitWavelengthStartNm, fitWavelengthStopNm);

        [fitResult, stats] = fitITFFano([lambda.'; transmission.'], false);
        m = fittedDipMetrics(fitResult, lambda, fitGridPoints, fwhmGridPoints);

        fitRows = [fitRows; table(i, spectra(i).thicknessNm, true, "", ...
            m.lambdaMin, m.fwhm, fitResult.q, m.visibility, ...
            fitResult.centre, fitResult.gamma, fitResult.offset, fitResult.amp, ...
            m.tMin, m.tMax, stats.gof.rsquare, stats.sse, ...
            'VariableNames', {'spectrumIndex', 'thicknessNm', 'fitOk', ...
            'fitMessage', 'resonanceWavelengthNm', 'fwhmNm', 'fanoQ', ...
            'visibility', 'fanoCentreNm', 'fanoGammaNm', 'fanoOffset', ...
            'fanoAmplitude', 'fittedMinTransmission', ...
            'fittedMaxTransmission', 'rSquared', 'sse'})]; %#ok<AGROW>

        if makePlots
            plotOneFit(lambda, transmission, fitResult, spectra(i).thicknessNm, ...
                outputFolder, i);
        end
end

changeRows = makeChangeTable(fitRows);

writetable(fitRows, fullfile(outputFolder, 'fano_thickness_fit_results.csv'));
writetable(changeRows, fullfile(outputFolder, 'fano_parameter_changes.csv'));

if makePlots
    plotParameterChanges(changeRows, outputFolder);
end

disp(fitRows);
disp(changeRows);

function [lambda, transmission] = selectWavelengthWindow(lambda, transmission, lowerNm, upperNm)
    if isempty(lowerNm) && isempty(upperNm)
        return;
    end

    keep = true(size(lambda));
    if ~isempty(lowerNm)
        keep = keep & lambda >= lowerNm;
    end
    if ~isempty(upperNm)
        keep = keep & lambda <= upperNm;
    end

    lambda = lambda(keep);
    transmission = transmission(keep);
end

function spectra = readThicknessSpectra(dataFile, pointsPerSpectrum, thicknessNm)

    data = readmatrix(dataFile, 'FileType', 'text', 'CommentStyle', '%');
    data = data(all(isfinite(data), 2), :);


    nSpectra = size(data, 1) / pointsPerSpectrum;

    if size(data, 2) >= 3
        rawThickness = data(:, 1);
        wavelength = data(:, 2);
        transmission = data(:, 3);
        thicknessNm = zeros(nSpectra, 1);
        for i = 1:nSpectra
            rows = blockRows(i, pointsPerSpectrum);
            thicknessNm(i) = mean(rawThickness(rows), 'omitnan');
        end
    else
        wavelength = data(:, 1);
        transmission = data(:, 2);
        if isempty(thicknessNm)
            thicknessNm = inferThicknessFromFileName(dataFile, nSpectra);
        end
        thicknessNm = thicknessNm(:);
    end


    spectra = struct('thicknessNm', num2cell(thicknessNm), ...
        'lambda', [], 'transmission', []);

    for i = 1:nSpectra
        rows = blockRows(i, pointsPerSpectrum);
        lambdaBlock = wavelength(rows);
        transmissionBlock = transmission(rows);
        [lambdaBlock, order] = sort(lambdaBlock);

        spectra(i).lambda = lambdaBlock(:);
        spectra(i).transmission = transmissionBlock(order);
    end
end

function rows = blockRows(blockIndex, pointsPerSpectrum)
    rows = (blockIndex - 1) * pointsPerSpectrum + (1:pointsPerSpectrum);
end

function thicknessNm = inferThicknessFromFileName(dataFile, nSpectra)
    [~, baseName] = fileparts(dataFile);
    token = regexp(baseName, 'thickness(?<start>\d+)-(?<stop>\d+)nm', ...
        'names', 'once');

    if isempty(token)
        thicknessNm = (0:nSpectra-1).';
        return;
    end

    startNm = parseCompactThickness(token.start);
    stopNm = parseCompactThickness(token.stop);
    thicknessNm = linspace(startNm, stopNm, nSpectra).';
end

function thicknessNm = parseCompactThickness(token)
    token = string(token);
    if strlength(token) == 3
        thicknessNm = str2double(token) / 100;
    else
        thicknessNm = str2double(token);
    end
end

function changeRows = makeChangeTable(fitRows)
    goodRows = fitRows(fitRows.fitOk, :);
    if isempty(goodRows)
        changeRows = fitRows(:, {'spectrumIndex', 'thicknessNm', ...
            'resonanceWavelengthNm', 'fwhmNm', 'visibility', 'fanoGammaNm'});
        return;
    end

    reference = goodRows(1, :);
    changeRows = goodRows(:, {'spectrumIndex', 'thicknessNm', ...
        'resonanceWavelengthNm', 'fwhmNm', 'visibility', 'fanoGammaNm'});

    changeRows.deltaThicknessNm = changeRows.thicknessNm - reference.thicknessNm;
    changeRows.deltaResonanceWavelengthNm = changeRows.resonanceWavelengthNm - ...
        reference.resonanceWavelengthNm;
    changeRows.deltaFwhmNm = changeRows.fwhmNm - reference.fwhmNm;
    changeRows.deltaVisibility = changeRows.visibility - reference.visibility;
    changeRows.deltaFanoGammaNm = changeRows.fanoGammaNm - reference.fanoGammaNm;

    if height(goodRows) >= 2
        p = polyfit(goodRows.thicknessNm, goodRows.resonanceWavelengthNm, 1);
        changeRows.resonanceSensitivityNmPerNm = repmat(p(1), height(changeRows), 1);
    else
        changeRows.resonanceSensitivityNmPerNm = NaN(height(changeRows), 1);
    end
end

function m = fittedDipMetrics(fitResult, lambda, dipGridPoints, fwhmGridPoints)
    gamma = abs(fitResult.gamma);
    span = max(lambda) - min(lambda);
    halfWindow = max([0.75*span, 20*gamma, 25]);
    lowerFit = min(min(lambda), fitResult.centre - halfWindow);
    upperFit = max(max(lambda), fitResult.centre + halfWindow);

    lambdaDipGrid = linspace(lowerFit, upperFit, dipGridPoints).';
    fitTransmission = fitResult(lambdaDipGrid);
    [~, idx] = min(fitTransmission);

    leftIdx = max(1, idx - 1);
    rightIdx = min(numel(lambdaDipGrid), idx + 1);
    lower = lambdaDipGrid(leftIdx);
    upper = lambdaDipGrid(rightIdx);
    if lower == upper
        lower = lowerFit;
        upper = upperFit;
    end

    [lambdaMin, tMin] = fminbnd(@(x) fitResult(x), lower, upper);
    m.lambdaMin = lambdaMin;
    m.tMin = tMin;
    m.fwhm = fittedFwhm(fitResult, lowerFit, upperFit, lambdaMin, tMin, ...
        fwhmGridPoints);
    m.tMax = max(fitTransmission);
    m.visibility = (m.tMax - m.tMin) / (m.tMax + m.tMin);
end

function fwhm = fittedFwhm(fitResult, lowerFit, upperFit, lambdaMin, tMin, gridPoints)
    lambdaFit = linspace(lowerFit, upperFit, gridPoints).';
    transmissionFit = fitResult(lambdaFit);
    baseline = max(transmissionFit);
    halfLevel = (baseline + tMin) / 2;
    [~, idxMin] = min(abs(lambdaFit - lambdaMin));

    leftIdx = find(transmissionFit(1:idxMin-1) > halfLevel & ...
        transmissionFit(2:idxMin) <= halfLevel, 1, 'last');
    rightIdx = idxMin - 1 + find(transmissionFit(idxMin:end-1) <= halfLevel & ...
        transmissionFit(idxMin+1:end) > halfLevel, 1, 'first');

    if isempty(leftIdx) || isempty(rightIdx)
        fwhm = NaN;
        return;
    end

    leftHalf = interp1(transmissionFit(leftIdx:leftIdx+1), ...
        lambdaFit(leftIdx:leftIdx+1), halfLevel);
    rightHalf = interp1(transmissionFit(rightIdx:rightIdx+1), ...
        lambdaFit(rightIdx:rightIdx+1), halfLevel);
    fwhm = rightHalf - leftHalf;
end

function plotOneFit(lambda, transmission, fitResult, thicknessNm, outputFolder, index)
    lambdaFit = linspace(min(lambda), max(lambda), 1000).';

    figure('Color', 'w', 'Visible', 'off');
    plot(lambda, transmission, 'k.', 'DisplayName', 'Data'); hold on;
    plot(lambdaFit, fitResult(lambdaFit), 'b-', 'LineWidth', 1.5, ...
        'DisplayName', 'Fano fit');
    grid on;
    xlabel('Wavelength (nm)');
    ylabel('Transmission');
    title(sprintf('Spectrum %d, thickness %.4g nm', index, thicknessNm));
    legend('Location', 'best');
    saveas(gcf, fullfile(outputFolder, sprintf('fano_fit_%03d.svg', index)));
    close(gcf);
end

function plotParameterChanges(changeRows, outputFolder)
    if isempty(changeRows)
        return;
    end

    figure('Color', 'w', 'Visible', 'off');
    tiledlayout(2, 2, 'TileSpacing', 'compact');

    fitRows = table();

    nexttile;
    row = plotQuadraticChange(changeRows.thicknessNm, ...
        changeRows.deltaResonanceWavelengthNm, ...
        "deltaResonanceWavelengthNm", "\Delta resonance wavelength (nm)");
    fitRows = [fitRows; row];

    nexttile;
    row = plotQuadraticChange(changeRows.thicknessNm, changeRows.deltaFwhmNm, ...
        "deltaFwhmNm", "\Delta FWHM (nm)");
    fitRows = [fitRows; row];

    nexttile;
    row = plotQuadraticChange(changeRows.thicknessNm, changeRows.deltaVisibility, ...
        "deltaVisibility", "\Delta visibility");
    fitRows = [fitRows; row];

    nexttile;
    row = plotQuadraticChange(changeRows.thicknessNm, ...
        changeRows.deltaFanoGammaNm, "deltaFanoGammaNm", ...
        "\Delta Fano gamma (nm)");
    fitRows = [fitRows; row];

    saveas(gcf, fullfile(outputFolder, 'fano_parameter_changes.svg'));
    savefig(gcf, fullfile(outputFolder, 'fano_parameter_changes.fig'));
    writetable(fitRows, fullfile(outputFolder, ...
        'fano_parameter_quadratic_fits.csv'));
    save(fullfile(outputFolder, 'fano_parameter_quadratic_fits.mat'), ...
        'fitRows');
    close(gcf);
end

function fitRow = plotQuadraticChange(x, y, parameterName, yLabelText)
    x = x(:);
    y = y(:);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    plot(x, y, 'o', 'LineWidth', 1.4, 'DisplayName', 'Fitted data');
    hold on;

    if numel(x) >= 3
        p = polyfit(x, y, 2);
        xFit = linspace(min(x), max(x), 300).';
        yFit = polyval(p, xFit);
        equation = sprintf('y = %.6g x^2 %+.6g x %+.6g', ...
            p(1), p(2), p(3));
        matlabFunction = sprintf('@(x) %.16g*x.^2 %+.16g*x %+.16g', ...
            p(1), p(2), p(3));

        plot(xFit, yFit, '-', 'LineWidth', 1.6, ...
            'DisplayName', sprintf('2nd order fit: a = %.4g, b = %.4g', ...
            p(1), p(2)));
        text(0.03, 0.94, equation, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', 9, ...
            'BackgroundColor', 'w', 'EdgeColor', [0.75 0.75 0.75], ...
            'Margin', 3);
    else
        p = [NaN, NaN, NaN];
        equation = "Not enough points for quadratic fit";
        matlabFunction = "";
    end

    grid on;
    xlabel('Thickness (nm)');
    ylabel(yLabelText);
    legend('Location', 'best');

    fitRow = table(parameterName, p(1), p(2), p(3), string(equation), ...
        string(matlabFunction), ...
        'VariableNames', {'parameterName', 'quadraticA', 'linearB', ...
        'constantC', 'equation', 'matlabFunction'});
end
