classdef SLM
% Class definition for SLM calculation and analysis functions.
% SLM = single level model
%===============================================================================
% Authors:  Artur Erbe, Bernd Briechle [Octave Script & Functions]
%           Matthias Wieser [MATLAB Class]
% Version:  0.31 [2013-02-28]
%===============================================================================
  methods (Static)



    function [FitModel,FitGOF,FitInfo] = fit (V,I,StartPoints,Weights,MaxIter,MaxFunEvals,TolFun,TolX)
        %GOF = goodness-of-fit
        % sse:          Sum of squares due to error
        % rsquare:      Coefficient of determination
        % dfe:          Degrees of freedom
        % adjrsquare:   Degree-of-freedom adjusted coefficient of determination
        % rmse:         Root mean squared error (standard error)

        %GOF DETAILS:
        % sse:          A value closer to 0 indicates that the model has a smaller random error component, and that the fit will be more useful for prediction.
        % rsquare:      R-square can take on any value between 0 and 1, with a value closer to 1 indicating that a greater proportion of variance is accounted for by the model.
        % adjrsquare:   The adjusted R-square statistic can take on any value less than or equal to 1, with a value closer to 1 indicating a better fit. Negative values can occur when the model contains terms that do not help to predict the response.
        % rmse:         Also known as the fit standard error and the standard error of the regression. A value closer to 0 indicates a fit that is more useful for prediction.

        [FitModel,FitGOF,FitInfo] = fit(V, I, fittype('SLM.current(E0, Gamma1, Gamma2, V)','independent','V','coefficients',{'E0','Gamma1','Gamma2'}), 'Algorithm', 'Levenberg-Marquardt', 'MaxIter', MaxIter, 'MaxFunEvals', MaxFunEvals, 'TolFun', TolFun, 'TolX', TolX, 'StartPoint', StartPoints, 'Weights', Weights);
    end



	function [Ifit,IfitPred] = Ifit (FitModel, V)
        Ifit = feval(FitModel, V);
        IfitPred = predint(FitModel, V, 0.95, 'observation', 'off');
    end



    function Transmission = transmission (E, E0, Gamma1, Gamma2, V)
        Transmission = 4 .* Gamma1 .* Gamma2 ./ ((E - (E0 + ((Gamma1 - Gamma2) ./ (Gamma1 + Gamma2)) * V./2)).^2 + (Gamma1 + Gamma2).^2);
    end



    function Fermi = fermi (E, mu, b)
        Fermi = 1 ./ (1 + exp((E - mu) .* b));
    end



    function Current = current (E0, Gamma1, Gamma2, V)
        % b = (kB*T)^-1 = 40
        % kB*T=25meV, INFO: 1/(40 * 8.6173324*10^-5) = 290.11 K
        b = 40;

        Current=struja(E0, Gamma1, Gamma2, V);
    end



  end %methods
end %classdef


function Struja= struja(E0, Gamma1, Gamma2, V) 
        Struja = zeros(size(V)); %Preallocation
        b = 40;
       %2do: GET "Parallel Computing Toolbox" AND USE "matlabpool" + "parfor" (SHOULD BE FASTER!)
            
            
            Struja =gather(arrayfun(@(Z) (7.74809174e-5 .* quadgk(@(E) SLM.transmission(E,E0,Gamma1,Gamma2,Z) .* ( SLM.fermi(E, Z./2, b) - SLM.fermi(E, -Z./2, b) ),-2,2)),V));
        
end