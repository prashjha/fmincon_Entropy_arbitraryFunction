classdef AnalyticStrategy 
    %ANALYTICSTRATEGY the first pool is fed from a measured signal that
    %needs to be excitation angle corrected
    %   Detailed explanation goes here
    
    properties
    end
    methods (Static)
        function paramsOut = parseParams(paramsIn)
            % parseParams: a function to fill default param values if they are
            % not defined
            default = struct('M0',[],'inputTime',[0,1],...
                'inputSignal',[0,0]);
            tmpNames = fieldnames(default);
            paramsOut = paramsIn;
            for i = 1:numel(tmpNames)
                if ~isfield(paramsOut,tmpNames{i})
                    paramsOut.(tmpNames{i}) = default.(tmpNames{i});
                end
            end
        end
    end
    methods
        function [TRList,Mxy,Mz] = compile(self,TRList,FaList,M0,A,params)
            % COMPILE: runs the model based on some input parameters
            if ~isempty(params.M0)
                M0 = params.M0;
            end
            inputTime = params.inputTime;
            inputSignal = params.inputSignal;
            [TRList, Mxy, Mz] = self.evaluate(TRList,FaList,M0,A,inputTime,inputSignal);
        end
    end
    methods(Access = private)
        function [TRList, Mxy, Mz] = evaluate(~,TRList,FaList,M0,A,inputTime,inputSignal)
            % EVALUATE: runs the model based on some input parameters     
            Mz = zeros(size(FaList));
            Mxy = zeros(size(FaList));
            M0(1) = inputSignal(1); % Set first timepoint to input signal
            Mz(:,1) = M0.*cos(FaList(:,1));
            Mxy(:,1) = (M0.*sin(FaList(:,1)));
            function dydt = tmpOdeFun(t,y)
                y(1) = interp1(inputTime,inputSignal,t,[],'extrap');
                dydt = A*y;
            end
            for i = 2:length(TRList)
                fun = @(t,y)tmpOdeFun(t,y);
                [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(:,i-1));
                Mz(:,i) = Y(end,:).';
                Mz(1,i) = interp1(inputTime,inputSignal,TRList(i),[],'extrap');
                Mxy(:,i) = sin(FaList(:,i)).*Mz(:,i);
                Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
            end
        end
    end 
end

