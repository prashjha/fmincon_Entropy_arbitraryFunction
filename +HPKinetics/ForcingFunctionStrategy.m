classdef ForcingFunctionStrategy
    %CLOSEDSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function paramsOut = parseParams(paramsIn)
            % parseParams: a function to fill default param values if they are
            % not defined
            default = struct('M0',[],'forcingFunction',@(t)0,...
                'scale',1);
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
            b = @(t)params.scale.*params.forcingFunction(t);
            [TRList, Mxy, Mz] = self.evaluate(TRList,FaList,M0,A,b);
        end
    end
    methods(Access = private)
        function [TRList, Mxy, Mz] = evaluate(~,TRList,FaList,M0,A,b)
            % EVALUATE: runs the model based on some input parameters
            Mz = zeros(size(FaList));
            Mxy = zeros(size(FaList));
            Mz(:,1) = M0.*cos(FaList(:,1));
            Mxy(:,1) = (M0.*sin(FaList(:,1)));
            fun = @(t,y)A*y+b(t);
            for i = 2:length(TRList)
                [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(:,i-1));
                Mz(:,i) = Y(end,:).';
                Mxy(:,i) = sin(FaList(:,i)).*Mz(:,i);
                Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
            end
        end
    end
    
end

