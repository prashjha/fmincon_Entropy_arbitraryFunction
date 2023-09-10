classdef SliceDecomp < handle
    %FITMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fitType
        model
    end
    %% Passthrough functions
    methods
        function defaults(self)
            self.model.defaults()
        end
        function [TRList,Mxy,Mz] = compile(self,M0,params)
            % COMPILE: runs the model based on some input parameters
            [TRList,Mxy,Mz] = self.model.compile(M0,params);
        end
        function M0out = splitM0(self,M0in,params)
            M0out = self.model.splitM0(M0in,params);
        end
    end
    methods (Access = private)
        function [xNames, xIndex,x0] = fillGuess(self,guess)
            [xNames, xIndex,x0] = self.model.fillGuess(guess);
        end
        function [resultParams,allParams] = packResults(...
                self,resultParams,allParams,xNames,xIndex,x,linkNames,linkIndex,link)
            [resultParams,allParams] = self.model.packResults(...
                resultParams,allParams,xNames,xIndex,x,linkNames,linkIndex,link);
        end
    end
    %% Functions that account for sclie profile
    methods
        function self = SliceDecomp(model)
            self.model = model;
        end
        function DataCompare(self,params,xdata,ydata,varargin)
            p = inputParser();
            p.addParameter('axis',[])
            p.parse(varargin{:})
            input_axis = p.Results.axis;
            nDecomp = params.nDecomp;
            flipAngleProfiles = params.flipAngleProfiles;
            [decompFlipAngles,tmpM0] = self.decompSlice(ydata(:,1),nDecomp,flipAngleProfiles);
            tmpM0 = self.splitM0(tmpM0,params);
            for m = 1:nDecomp
                params.FaList = decompFlipAngles(:,m);
                [TRList,Mxy(m,:,:),~] = self.compile(tmpM0,params);
            end
            Mxy = nDecomp*squeeze(sum(Mxy,1));
            % Sometimes the squeze will transpose the matrix
            if(size(Mxy,2)==1)
                Mxy = Mxy.';
            end
            legendVals = cell(size(Mxy,1),1);
            if isempty(input_axis)
                figure
            else
                axis(input_axis);
            end
            for i = 1:size(Mxy,1)
                tmpLine = plot(TRList,Mxy(i,:));
                hold on
                plot(xdata,ydata(i,:),'o','MarkerEdgeColor',tmpLine.Color);
                legendVals{2*i-1} = ['Fit Pool ',char(i+'A'-1)];
                legendVals{2*i} = ['Data Pool ',char(i+'A'-1)];
            end
            hold off
            xlabel('Time (sec)')
            ylabel('Signal (arb)')
            legend(legendVals)
        end
        function [x,resultParams,allParams,resnorm,residual,exitflag,output,lambda,jacobian]...
                = fitData(self,params,guess,xdata,ydata,varargin)
            %FITDATA: Fits some set of guess parameters to input data
            %following the given model
            p = inputParser();
            p.addParameter('lb',[])
            p.addParameter('ub',[])
            p.addParameter('linker',[])
            p.parse(varargin{:})
            linker = p.Results.linker;
            % Fill fit parameters from the guesses
            [xNames, xIndex,x0] = self.model.fillGuess(guess);
            if ~isempty(linker)
                % Fill link parameters
                [linkNames, linkIndex,link] = self.model.fillGuess(linker);
            else
                linkNames = [];
                linkIndex = [];
                link = [];
            end
            if ~isempty(linker)
                fun = @(x,xdata)self.fitFunction(...
                    params,x,xNames,xIndex,ydata(:,1),link,linkNames,linkIndex);
            else
                fun = @(x,xdata)self.fitFunction(params,x,xNames,xIndex,ydata(:,1));
            end
            opts = params.fitOptions;
            [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                lsqcurvefit(fun,x0,xdata,ydata,...
                [p.Results.lb],[p.Results.ub],opts);
            resultParams = guess;
            allParams = params;
            % Pack up Fit and link Parameters
            [resultParams,allParams] = self.model.packResults(...
                resultParams,allParams,xNames,xIndex,x,linkNames,linkIndex,link);
        end
        function Y = fitFunction(self,params,x,xNames,xIndex,Y0,varargin)
            % fitFunction packs the parameter in params and x up and evaluates
            % using the evaluate funnction over some time (tSpan) with some
            % initial value (Y0)
            % Parse input
            p = inputParser();
            p.addOptional('link',[])
            p.addOptional('linkNames',[])
            p.addOptional('linkIndex',[])
            p.parse(varargin{:})
            link = p.Results.link;
            linkNames = p.Results.linkNames;
            linkIndex = p.Results.linkIndex;
            nDecomp = params.nDecomp;
            flipAngleProfiles = params.flipAngleProfiles;
            [decompFlipAngles,tmpY0] = self.decompSlice(Y0,nDecomp,flipAngleProfiles);
            
            tmpY0 = self.splitM0(tmpY0,params);
            for m = 1:nDecomp
                j=1;
                params.FaList = decompFlipAngles(:,m);
                % Fill fit variables
                for i = 1:numel(xNames)
                    % fill in the rest of the fit variables
                    for k = 1:numel(xIndex{i})
                        params.(xNames{i})(xIndex{i}(k)) = x(j);
                        % Check if fitting flip angle (there mus be a better
                        % way to do this
                        if strcmp(xNames{i}, 'FaList')
                            params.(xNames{i}) =...
                                repmat(x(j),size(params.(xNames{i})));
                        end
                        j = j+1;
                    end
                end
                if ~isempty(link)
                    j=1;
                    % Link multiple varibles to fit variables
                    for i = 1:numel(linkNames)
                        % fill in the rest of the fit variables
                        for k = 1:numel(linkIndex{i})
                            % DO NOT LINK FLIP ANGLE THIS PROBABLY WONT WORK
                            params.(linkNames{i})(linkIndex{i}(k)) = x(link(j));
                            j = j+1;
                        end
                    end
                end
                % Split out the multiple Physical compartments
                params = self.parseParams(params);
                [~, Y(m,:,:), ~] = self.compile(tmpY0,params);
            end
            Y = nDecomp*squeeze(sum(Y,1)).';
            if size(Y,1)~=1
                Y= Y.';
            end
        end
        function paramsOut = parseParams(self,paramsIn)
            % parseParams: a function to fill default param values if they are
            % not defined
            default = struct('nDecomp',64);
            tmpNames = fieldnames(default);
            paramsOut = paramsIn;
            for i = 1:numel(tmpNames)
                if ~isfield(paramsOut,tmpNames{i})
                    paramsOut.(tmpNames{i}) = default.(tmpNames{i});
                end
            end
            paramsOut = self.model.parseParams(paramsOut);
        end
    end
    methods (Access = private)
        function [decompFlipAngles,tmpY0] = decompSlice(~,Y0,nDecomp,flipAngleProfiles)
            for i = 1: size(flipAngleProfiles,1)
                decompFlipAngles(i,:) = interp1(...
                    1:length(flipAngleProfiles(i,:)),flipAngleProfiles(i,:),...
                    linspace(1,length(flipAngleProfiles(i,:)),nDecomp));
                tmpY0(i,1) = Y0(i)/sum(sin(decompFlipAngles(i,:)))/nDecomp;
            end
        end
    end
    
end

