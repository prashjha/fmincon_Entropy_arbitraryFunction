classdef MultiSliceDecomp < HPKinetics.SliceDecomp
    %FITMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    %% Passthrough functions
    methods
        function defaults(self)
            self.model.defaults()
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
        function self = MultiSliceDecomp(model)
            self = self@HPKinetics.SliceDecomp(model);
        end
        function DataCompare(self,params,xdata,ydata,varargin)
            p = inputParser();
            p.addParameter('axis',[])
            p.parse(varargin{:})
            input_axis = p.Results.axis;
            [TRList,Mxy,~] = self.compile(ydata(:,1),params);
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
        function [TRList,Mxy,Mz] = compile(self,M0,params)
            params = self.parseParams(params);
            if params.nSlices > 1
                % Fill flip angle list
                params.FaList = repmat(params.FaList.',[1,length(params.TRList)]);
                % Fill TR list with dummy timepoints for the other slices
                fullTRList = repmat(params.TRList,[params.nSlices,1])+((((1:params.nSlices).')-1)*params.SliceTR);
                fullTRList = reshape(fullTRList,[1,numel(fullTRList)]);
                params.TRList =  fullTRList;
            end
            % COMPILE: runs the model based on some input parameters
            nDecomp = params.nDecomp;
            sliceProfile = params.sliceProfile;
            [decompSliceProfile] = self.decompSlice(nDecomp,sliceProfile,params.nSlices,params.nChemical);
            if ~isempty(params.M0),M0 = params.M0;end
            tmpY0 = self.decompY0(M0,squeeze(decompSliceProfile(params.SliceNo,:,:)),params.FaList);
            tmpY0 = self.splitM0(tmpY0,params);
            tmpFaList = params.FaList;
            for m = 1:nDecomp
                tmpParams = params;
                params.FaList = repmat(squeeze(decompSliceProfile(:,:,m)),[1,size(tmpFaList,2)/size(squeeze(decompSliceProfile(:,:,m)),2)]).*tmpFaList; % Scale the flip angles
                params = self.parseParams(params);
                params = self.model.parseParams(params);
                [TRList,Mxy(m,:,:),Mz(m,:,:)] = self.model.compile(tmpY0,params);
                params = tmpParams;
            end
            Mxy = squeeze(sum(Mxy,1));
            Mz = squeeze(sum(Mz,1));
            % Sometimes the squeze will transpose the matrix
            if(size(Mxy,2)==1)
                Mxy = Mxy.';
            end
            if(size(Mz,2)==1)
                Mz = Mz.';
            end
            if params.nSlices > 1
            % Remove dummy timepoints
            TRList = TRList(params.SliceNo:params.nSlices:end);
            Mxy = Mxy(:,params.SliceNo:params.nSlices:end);
            Mz = Mz(:,params.SliceNo:params.nSlices:end);
            end
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
            params = self.parseParams(params);
            if nDecomp < size(params.sliceProfile,1),params.nDecomp = size(params.sliceProfile,1); end % Make sure there is at least 1 sub-slice for each slice;
                j=1;
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
                [~, Y, ~] = self.compile(Y0,params);
        end
        function paramsOut = parseParams(self,paramsIn)
            % parseParams: a function to fill default param values if they are
            % not defined
            default = struct('nDecomp',64,'SliceTR',0.1,'SliceNo',1);
            tmpNames = fieldnames(default);
            paramsOut = paramsIn;
            for i = 1:numel(tmpNames)
                if ~isfield(paramsOut,tmpNames{i})
                    paramsOut.(tmpNames{i}) = default.(tmpNames{i});
                end
            end
            if ndims(paramsOut.FaList)>3
            nSlices = size(paramsOut.FaList,1); % The number of slices must be reflected in the slice profiles
            else
                nSlices = 1;
            end
            paramsOut.nSlices = nSlices;
            paramsOut.nChemical = size(paramsIn.ExchangeTerms,2);
            sliceTR = paramsOut.SliceTR;
        end
    end
    methods (Access = private)
        function [decompSliceProfile] = decompSlice(~,nDecomp,sliceProfile,nSlices,nChemical)
            if size(sliceProfile,2)~=nChemical % Assume no chemical shift was passed in
                if ndims(sliceProfile)>2, warning('The Dimensionality of the slice profile in in an unexpected form. Only using the first 2 Dimensions/n'); end
                sliceProfile = repmat(squeeze(sliceProfile(:,:,1)),[1,1,nChemical]);
                sliceProfile = permute(sliceProfile,[1,3,2]);
            end
            for j = 1:nSlices % First dimension will be the slices
            for i = 1:nChemical % Second dimension will be the chemial species
                % 1 D interp the input slice profiles to match the
                % requested number of sub-slices
                decompSliceProfile(j,i,:) = interp1(...
                    1:length(sliceProfile(j,i,:)),squeeze(sliceProfile(j,i,:)),...
                    linspace(1,length(sliceProfile(j,i,:)),nDecomp));
            end
            end
        end
        function [decompY0] = decompY0(~,Y0,decompSliceProfile,FaList)
            % Devide the initial magnetization over the sub slices. Not
            % sure if this is the correct method
            decompY0 = Y0./sum(sin(decompSliceProfile.*FaList(:,1)),2);
        end
    end
end

