classdef GammaVariateStrategy < HPKinetics.VIFStrategy
    %VIFSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function paramsOut = parseParams(self,paramsIn)
            % parseParams: a function to fill default param values if they are
            % not defined
            default = struct('gammaPDFA',2.8,'gammaPDFB',4.5);
            tmpNames = fieldnames(default);
            paramsOut = paramsIn;
            for i = 1:numel(tmpNames)
                if ~isfield(paramsOut,tmpNames{i})
                    paramsOut.(tmpNames{i}) = default.(tmpNames{i});
                end
            end
                        paramsOut.VIF = @(t) gampdf(t,paramsOut.gammaPDFA,paramsOut.gammaPDFB);
%             paramsOut.VIF = @(t) t.^(paramsOut.gammaPDFA-1).*exp(-t./paramsOut.gammaPDFB);
            paramsOut = parseParams@HPKinetics.VIFStrategy(self,paramsOut);
        end
    end
end

