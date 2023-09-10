classdef NewMultiPoolTofftsGammaVIF < HPKinetics.MultiPoolToffts
    %MULTIPOOLTOFFTSGAMMAVIF A chemical exchange model assuming two pooled
    %Tofts model of perfusion
    %   parameters Values
    %* ExchangeTerms - A Matrix defining chemical Exchange. Defalt: 0
    %* T1s - A row vector of T1 decay terms. Default: 100
    %* FaList - A matrix of excitation angles. Default: 0
    %* TRList - A matrix of excitation times. Default: 0
    %* t0 - A row vector for delivery delay of each metabolite. Default: 0
    %* gammaPdfA - A row vector for shape term alpha of each metabolite. Default: 2.8
    %* gammaPdfB - A row vector for shape term beta of each metabolite. Default: 4.5
    %* ScaleFactor - A row vector for each metabolite's VIF scale factor. Default: 1
    %* fitOptions - A matlab fit option structure. Default: optimset(''lsqcurvefit'')
    %* PerfusionTerms - A row vector for each metabolite's extravisation rate. Default: 0
    %* volumeFractions - A row vector for each metabolite's volume fraction. Default: 1
    %   There is NO imput validation for the parameters passed in, for more
    %   detail on the assumed data structur of these parameters use the
    %   defaults function
    properties
    end
    methods
        function defaults(self)
            % DEFAULTS explains the default values for each parameter
            names = {'t0','gammaPdfA','gammaPdfB','scaleFactor'};
            discriptions = {'A  Row vector of time delays for each metabolite'...
                ' A  Row vector of shape term Alpha, set this to zero to have no VIF for a chemical pool'...
                ' A  Row vector of shape term Beta, this cannot be zero and will be set to 1e-40 if zero is used'...
                ' A  Row vector of Scale Factor to be applied to the VIF'};
            defaultsVals = {'0','2.8','4.5','1'};
            fprintf('*Note* all terms must be a vector of size 1 x N where N is the number of chemical Pools\n')
            for i = 1:numel(names)
                fprintf('''%s'': %s\n Default Vaule: %s\n',...
                    names{i},discriptions{i},defaultsVals{i});
            end
            defaults@HPKinetics.MultiPoolToffts(self);
        end
        function paramsOut = parseParams(self,paramsIn)
            % parseParams: Parses the input shape terms of a gamma variate
            % for the VIF, each term should be a vector with shape terms
            % for each chemical species
            
            % Fill Default Values
            default = struct('t0',0,'gammaPdfA',2.8,...
                'gammaPdfB',4.5,'scaleFactor',1);
            tmpNames = fieldnames(default);
            paramsOut = paramsIn;
            for i = 1:numel(tmpNames)
                if ~isfield(paramsOut,tmpNames{i})
                    paramsOut.(tmpNames{i}) = default.(tmpNames{i});
                end
            end
            % Build VIF
            paramsOut.VIF = @(t)paramsIn.scaleFactor.*...
                gampdf(t-paramsIn.t0,paramsIn.gammaPdfA,paramsIn.gammaPdfB);
            paramsOut.gammaPdfA = paramsIn.gammaPdfA;
            paramsOut.gammaPdfB = paramsIn.gammaPdfB;
            % Fill in parent Class defaults
            paramsOut = parseParams@HPKinetics.MultiPoolToffts(self,paramsOut);
        end
        function [A] = getA(self,params)
            params = self.parseParams(params);
            A = params.A;
        end
        function [VIF] = getVIF(self,params)
            params = self.parseParams(params);
            b = params.b;
            FaList = params.FaList;
            TRList = params.TRList;
            N = size(TRList,2);
            
            VIF = zeros(size(FaList));
            for i = 1:N
                VIF(:,i) = b(TRList(i));
            end
        end
        function [TRList,Mxy,Mz] = compile(self,M0,params)
        % EVALUATE: runs the model based on some input parameters
        %params = self.parseParams(params);
        % EVALUATE: runs the model based on some input parameters
            params = self.parseParams(params);
            A = params.A;
            b = params.b;
            FaList = params.FaList;
            TRList = params.TRList;
            [TRList, Mxy, Mz] = self.evaluate(...
                TRList,FaList,M0,A,b,params);
            %[TRList,Mxy,Mz] = compile@HPKinetics.MultiPoolToffts(self,M0,params);
        end
        
        function [TRList,Mxy,Mz,G,dGdTR,dGdFA] = compile_der(self,M0,params)
        % EVALUATE: runs the model based on some input parameters and 
        % computes derivative of signal function wrt TR and FA 
        %params = self.parseParams(params);
        % EVALUATE: runs the model based on some input parameters
            params = self.parseParams(params);
            A = params.A;
            b = params.b;
            FaList = params.FaList;
            TRList = params.TRList;
            [TRList, Mxy, Mz, G, dGdTR, dGdFA] = self.evaluate_der(...
                TRList,FaList,M0,A,b,params);
        end
    end
        
    methods (Access = private)
         function [TRList, Mxy, Mz] = evaluate(self,TRList,FaList,M0,A,b,params)
            % EVALUATE: runs the model based on some input parameters 
            kve = params.kve;
            ve = params.ve;
            N = size(TRList,2);
            
            Mz = zeros(size(FaList));
            Mxy = zeros(size(FaList));
            
            % first step
            Mz(:,1) = M0.*cos(FaList(:,1));
            Mxy(:,1) = (params.ve*M0+(1-params.ve)*params.b(TRList(1))).*sin(FaList(:,1));
            
            walkermethod = true; % <-- ode45 is faster than using integral in second method
            if walkermethod
                fun = @(t,y)A*y+(kve.'/ve).*b(t);
                for i = 2:N
                    [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(:,i-1));
                    Mz(:,i) = Y(end,:)';
                    Mxy(:,i) = sin(FaList(:,i)).*(params.ve.*Mz(:,i)+...
                        (1-params.ve).*b(TRList(i)));
                    Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
                end
            else
                % compute TRs from TRList
                TR = zeros(N,1);
                TR(1) = TRList(1);
                for i = 2:N
                    TR(i) = TRList(i) - TRList(i-1);
                end
                
                fun = @(t,tk) expm((tk - t)*A)*b(t);
                for i = 2:N
                    ti = TRList(i);
                    tim = TRList(i-1);
                    
                    Mz(:,i) = expm(TR(i)*A)*Mz(:,i-1) + (kve.'/ve).*integral(@(x) fun(x,ti), tim, ti,'ArrayValued',true,'RelTol',1e-6);
                    Mxy(:,i) = sin(FaList(:,i)).*(params.ve.*Mz(:,i)+...
                        (1-params.ve).*b(TRList(i)));
                    Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
                end
            end
         end
         
         function [TRList, Mxy, Mz, G, dGdTR, dGdFA] = evaluate_der(self,TRList,FaList,M0,A,b,params)
            % EVALUATE: runs the model based on some input parameters and 
            % computes derivative of signal function wrt TR and FA
            kve = params.kve;
            ve = params.ve;
            N = size(TRList,2);
            gpdfa = params.gammaPdfA(1); % for pyruvate
            gpdfb = params.gammaPdfB(1); % for pyruvate
            
            % Note that here TRList[i] is taken as time t_i so it
            % should be sum_{k=1}^i TR_i
            % compute TRs from TRList
            TR = zeros(N,1);
            TR(1) = TRList(1);
            for i = 2:N
                TR(i) = TRList(i) - TRList(i-1);
            end
            
            Mz = zeros(size(FaList));
            Mz_nocos = zeros(size(FaList));  % without post cosine multiplication for derivative calc
            Mxy = zeros(size(FaList));
            
            % first step
            Mz(:,1) = M0.*cos(FaList(:,1));
            Mz_nocos(:,1) = M0;
            Mxy(:,1) = (ve*M0+(1-ve)*b(TRList(1))).*sin(FaList(:,1));
            
            walkermethod = true; % <--- ode45 is faster than using integral in second method
            if walkermethod
                % @cmwalker - matrix transpose nightmare here - please review - or is there a later version ? 
                fun = @(t,y)A*y+(kve.'/ve).*b(t);
                for i = 2:N
                    [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(:,i-1));
                    Mz(:,i) = Y(end,:)';
                    % @cmwalker - this doesn't match your disseration or 2019 paper?
                    % @cmwalker - measured signal lis the intracellular component +  any signal in the extracellular/extravascular space ?
                    Mxy(:,i) = sin(FaList(:,i)).*(ve.*Mz(:,i)+(1-ve).*b(TRList(i)));
                    Mz_nocos(:,i) = Mz(:,i);
                    Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
                end
            else
                fun = @(t,tk) expm((tk - t)*A)*b(t);
                for i = 2:N
                    ti = TRList(i);
                    tim = TRList(i-1);
                    
                    Mz(:,i) = expm(TR(i)*A)*Mz(:,i-1) + (kve.'/ve).*integral(@(x) fun(x,ti), tim, ti,'ArrayValued',true);
                    Mxy(:,i) = sin(FaList(:,i)).*(ve.*Mz(:,i)+(1-ve).*b(TRList(i)));
                    Mz_nocos(:,i) = Mz(:,i);
                    Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
                end
            end
            
            % adjoint step (adjoint is vector of size N+1, where last 
            % element of vector is zero
            Q = zeros(size(FaList,1), N+1); % adjoint
            Q(:,N+1) = zeros(2,1); % init adjoint
            for i = N:-1:1
                fi = FaList(:,i);
                TRip = 0;
                if i < N
                    TRip = TR(i+1);
                end
                Q(:,i) = cos(fi).*(expm(TRip*A).'*Q(:,i+1)) + ve*sin(fi);
            end
            
            % signal function and derivatives
            G =  sum(Mxy(:));
            dGdTR = zeros(size(TRList));
            dGdFA = zeros(size(FaList));
            
            % loop for derivatives
            for i = 1:N
                ti = TRList(i);
                TRi = TR(i);
                TRip = 0;
                if i < N
                    TRip = TR(i+1);
                end
                Pi = FaList(1,i);
                Li = FaList(2,i);
                Atri = expm(TRi*A);
                Atrip = expm(TRip*A);
                
                % contribution from der of sin term
                xx = (1-ve)*b(ti) + ve*Mz_nocos(:,i);
                
                ds = zeros(2,1);
                ds(1) = cos(Pi);
                dGdFA(1,i) = dot(ds,xx);
                
                ds(1) = 0;
                ds(2) = cos(Li);
                dGdFA(2,i) = dot(ds,xx);
                
                % contribution from der of VIF term
                if i > 1
                    for j=i:N
                        tj = TRList(j);
                        ds = sin(FaList(:,j));
                        ds(2) = 0; % so that VIF for lactate is ignored
                        % derivative of gamma fn: dt g(t,a,b) =
                        % g(t,a,b)[(a-1)/t - 1/b]
                        dGdTR(:,i) = dGdTR(:,i) + (1-ve)*((gpdfa - 1)/tj - 1/gpdfb)*dot(ds,b(tj));
                    end
                end
                
                % contribution from der of fwd model to dGdFA
                if i < N
                    ds = zeros(2,1);
                    ds(1) = -sin(Pi);
                    dGdFA(1,i) = dGdFA(1,i) + dot(Q(:,i+1), Atrip*(ds.*Mz_nocos(:,i)));
                    
                    ds(1) = 0;
                    ds(2) = -sin(Li);
                    dGdFA(2,i) = dGdFA(2,i) + dot(Q(:,i+1), Atrip*(ds.*Mz_nocos(:,i)));
                end
                
                % contribution to dGdTR
                if i > 1
                    ds = cos(FaList(:,i-1)).*Mz_nocos(:,i-1);
                    dGdTR(:,i) = dGdTR(:,i) + dot(Q(:,i), A*Atri*ds);
               

                    % lastly, add contribution to dGdTR from der VIF term
                    for j=i:N
                        % compute derivative of \int_{t_{i-1}}^{t_i} \exp[A(t_i
                        % - \tau)] VIF[\tau] d\tau
                        tj = TRList(j);
                        tjm = TRList(j-1);

                        % term d t_j/d TR_i (from derivative of upper 
                        % integral limit in Leibniz rule)
                        db = (kve.'/ve).*b(tj);

                        % term d t_{j-1} / d TR_i
                        if j-1 >= i
                            db = db - expm((tj-tjm)*A)*((kve.'/ve).*b(tjm));
                        end

                        % from integral
                        % compute integral from solution instead of integrating again
                        bj = Mz_nocos(:,j) - expm(TR(j)*A)*(cos(FaList(:,j-1)).*Mz_nocos(:,j-1));
                        % or compute by integration
                        %bj = (kve.'/ve).*integral(@(x) fun(x,tj), tjm, tj,'ArrayValued',true);
                        db = db + A*bj;

                        dGdTR(:,i) = dGdTR(:,i) + dot(Q(:,j), db);
                    end
                end
            end
         end
    end
end
    