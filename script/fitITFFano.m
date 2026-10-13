function [f, stats] = fitITFFano(TF, showPlot, fitConfig)
% Simple Fano fit for one downward transmission resonance.
% TF(1,:) = wavelength in nm
% TF(2,:) = transmission

if nargin < 2
    showPlot = true;
end
if nargin < 3
    fitConfig = struct();
end

lambda = TF(1, :).';
T = TF(2, :).';
lambda = lambda(:);
T = T(:);

[lambda, order] = sort(lambda);
T = T(order);

[Tmin, idx] = min(T);
Tmax = max(T);
lambda0 = lambda(idx);
span = max(lambda) - min(lambda);
fwhm0 = estimateFWHM(lambda, T, idx);
gamma0 = max(fwhm0/2, span/1000);
depth0 = max(Tmax - Tmin, eps);

% Downward Fano dip:
% T = 1 - amp * (q + eps)^2 / ((1 + q^2) * (1 + eps^2))
% eps = (lambda - centre) / gamma
model = fittype( ...
    @(offset, amp, q, centre, gamma, lambda) ...
    offset - amp .* ((q + (lambda-centre)./gamma).^2 ./ ...
    ((1 + q.^2) .* (1 + ((lambda-centre)./gamma).^2))), ...
    'independent', 'lambda', ...
    'coefficients', {'offset', 'amp', 'q', 'centre', 'gamma'});

centreWindow = max([2*fwhm0, 2, 0.1*span]);
gammaLower = max(span/10000, eps);
gammaUpper = max(span, 10*gamma0);
offsetUpper = max([1.2*Tmax, Tmax + depth0, 1]);
ampUpper = max([2*depth0, Tmax - Tmin, 1]);

qStarts = [-50, -20, -10, -5, -2, 2, 5, 10, 20, 50];
gammaStarts = unique(max(gammaLower, [0.5, 1, 2, 4] * gamma0));
bestSse = Inf;
f = [];
bestGof = [];
bestOutput = [];

if isfield(fitConfig, 'StartPoint')
    startPoint = fitConfig.StartPoint(:).';
    localCentreWindow = max([centreWindow, abs(startPoint(5))*8, 1]);
    opts = fitoptions('Method', 'NonlinearLeastSquares', ...
        'StartPoint', startPoint, ...
        'Lower',      [0, 0, -200, startPoint(4)-localCentreWindow, gammaLower], ...
        'Upper',      [offsetUpper, ampUpper, 200, startPoint(4)+localCentreWindow, gammaUpper], ...
        'MaxIter', 800, ...
        'MaxFunEvals', 1600);

    [f, bestGof, bestOutput] = fit(lambda, T, model, opts);
    residual = T - f(lambda);
    bestSse = sum(residual.^2);
else
    for q0 = qStarts
        for g0 = gammaStarts
            % In this parameterisation the dip is near centre + gamma/q.
            centre0 = lambda0 - g0/q0;
            centre0 = min(max(centre0, lambda0 - centreWindow), lambda0 + centreWindow);

            opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                'StartPoint', [Tmax, depth0, q0, centre0, g0], ...
                'Lower',      [0,    0,     -200, lambda0-centreWindow, gammaLower], ...
                'Upper',      [offsetUpper, ampUpper, 200, lambda0+centreWindow, gammaUpper], ...
                'MaxIter', 10000, ...
                'MaxFunEvals', 20000);

                [candidate, candidateGof, candidateOutput] = fit(lambda, T, model, opts);
                residual = T - candidate(lambda);
                sse = sum(residual.^2);
                if sse < bestSse
                    bestSse = sse;
                    f = candidate;
                    bestGof = candidateGof;
                    bestOutput = candidateOutput;
                end
        end
    end
end


if nargout > 1
    stats = makeFitStats(f, bestGof, bestOutput, bestSse);
end

if showPlot
    lambdaFit = linspace(min(lambda), max(lambda), 2000).';
    figure('Color', 'w');
    plot(lambda, T, 'k.', 'DisplayName', 'Data'); hold on;
    plot(lambdaFit, f(lambdaFit), 'b-', 'LineWidth', 1.5, ...
        'DisplayName', 'Fano fit');
    grid on;
    xlabel('Wavelength (nm)');
    ylabel('Transmission');
    legend('Location', 'best');
    title(sprintf('centre = %.4f nm, gamma = %.4f nm, q = %.3f', ...
        f.centre, f.gamma, f.q));
end
end

function stats = makeFitStats(f, gof, output, sse)
stats = struct();
stats.gof = gof;
stats.output = output;
stats.sse = sse;
stats.coeffNames = coeffnames(f).';
stats.coeffValues = coeffvalues(f);

stats.conf95 = confint(f, 0.95);
stats.paramSigma = (stats.conf95(2, :) - stats.conf95(1, :)) / (2 * 1.96);
end

function fwhm = estimateFWHM(lambda, T, idx)
baseline = max(T);
halfLevel = (baseline + T(idx)) / 2;

leftIdx = find(T(1:idx-1) > halfLevel & T(2:idx) <= halfLevel, 1, 'last');
rightIdx = idx - 1 + find(T(idx:end-1) <= halfLevel & T(idx+1:end) > halfLevel, 1, 'first');

if isempty(leftIdx) || isempty(rightIdx)
    fwhm = 0.05 * (max(lambda) - min(lambda));
    return;
end

leftHalf = interp1(T(leftIdx:leftIdx+1), lambda(leftIdx:leftIdx+1), halfLevel);
rightHalf = interp1(T(rightIdx:rightIdx+1), lambda(rightIdx:rightIdx+1), halfLevel);
fwhm = rightHalf - leftHalf;
end
