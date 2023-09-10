classdef VIFStrategy
    %VIFSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function paramsOut = parseParams(~,paramsIn)
            % parseParams: a function to fill default param values if they are
            % not defined
            default = struct('VIF',@(t)0,'Extravasation',[0,0],...
                'bloodVolumeFraction',[],'VIFScale',1,'t0',0,'M0',[]);
            tmpNames = fieldnames(default);
            paramsOut = paramsIn;
            for i = 1:numel(tmpNames)
                if ~isfield(paramsOut,tmpNames{i})
                    paramsOut.(tmpNames{i}) = default.(tmpNames{i});
                end
            end
            paramsOut.VIF = paramsOut.VIF;
            %% Fill vascular voluem
            if isempty(paramsOut.bloodVolumeFraction)
                paramsOut.bloodVolumeFraction = 1-sum(paramsOut.volumeFractions);
            end
            if(paramsOut.bloodVolumeFraction <= 0)
                error('the volume fraction of the vascular space must be greater than zero. Current value: %f/n',...
                    paramsOut.bloodVolumeFraction)
            end
        end
        function [TRList,Mxy,Mz] = compile(self,TRList,FaList,M0,A,params)
            % COMPILE: runs the model based on some input parameters
            if ~isempty(params.M0)
                M0 = params.M0;
            end
            tmpVIF = params.VIF;
            bloodVolumeFraction = params.bloodVolumeFraction;
            tmpExtravasation = params.Extravasation;
            % Fill extravesation and VIF for compartments not in contact
            % with the vascular space
            zPadVIF = size(A,2)-length(tmpVIF(params.t0));
            if(length(params.VIFScale) == 1)
                VIF = @(t)params.VIFScale.*[tmpVIF(t-params.t0);zeros(zPadVIF,1)];
            else
                VIF = @(t)[params.VIFScale;zeros(zPadVIF,1)].*...
                    [tmpVIF(t-params.t0);zeros(zPadVIF,1)];
            end
            Extravasation = zeros(size(A,2),2);
            Extravasation(1:size(tmpExtravasation,1),:)=tmpExtravasation;
            [TRList, Mxy, Mz] = self.evaluate(TRList,FaList,M0,A,...
                Extravasation,bloodVolumeFraction,VIF);
        end
    end
    methods(Access = private)
        function [TRList, Mxy, Mz] = evaluate(~,TRList,FaList,M0,A,...
                Extravasation,bloodVolumeFraction,VIF)
            % EVALUATE: runs the model based on some input parameters
            A = A - diag(squeeze(Extravasation(:,2)./(1-bloodVolumeFraction)));
            fun = @(t,y)(A*y+(Extravasation(:,1)./(1-bloodVolumeFraction)).*VIF(t));
            Mz = zeros(size(FaList));
            Mxy = zeros(size(FaList));
            Mz(:,1) = M0.*cos(FaList(:,1));
            Mxy(:,1) = (M0.*sin(FaList(:,1)));
            for i = 2:length(TRList)
                % the transpose on Mz dose not matter as matlab
                % automatically converts the vector to match matrix
                % multiplication conventions but is done for clairity to
                % match with fun as it is declared above
                [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(:,i-1));
                Mz(:,i) = Y(end,:).';
                Mxy(:,i) = sin(FaList(:,i)).*((1-bloodVolumeFraction).*Mz(:,i)+...
                    bloodVolumeFraction.*VIF(TRList(i)));
                Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
            end
        end
    end
    
end

