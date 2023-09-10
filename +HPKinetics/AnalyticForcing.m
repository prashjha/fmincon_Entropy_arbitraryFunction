classdef AnalyticForcing < ForcingFunctionStrategy
    %CLOSEDSTRATEGY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    methods(Access = private)
        function [TRList, Mxy, Mz] = evaluate(~,TRList,FaList,M0,A,b)
            % EVALUATE: runs the model based on some input parameters     
            A11 = A(1,1); % Conver T1 for Pool A
            Mz = zeros(size(FaList));
            Mxy = zeros(size(FaList));
            M0(1) = b(1)/sin(FaList(1,1)); % Set first timepoint to input signal
            Mz(:,1) = M0.*cos(FaList(:,1));
            Mxy(:,1) = (M0.*sin(FaList(:,1)));
            input = zeros(size(Mz(:,1)));
            for i = 2:length(TRList)
                Mz(1,i) = b(i)/sin(FaList(1,i)); %% FA correct input signal
                TR = TRList(i)-TRList(i-1);
                input(1) = (Mz(1,i)-Mz(1,i-1)*exp((A11)*TR))*(-A11)/(1-exp((A11)*TR));
                fun = @(t,y)A*y+input;
                [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(:,i-1));
                Mz(:,i) = Y(end,:).';
                Mxy(:,i) = sin(FaList(:,i)).*Mz(:,i);
                Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
            end
        end
    end
    
end

