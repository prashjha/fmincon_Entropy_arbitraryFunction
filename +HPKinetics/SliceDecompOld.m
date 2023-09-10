classdef SliceDecompOld < HPKinetics.AbstractPools
    %SLICEDECOMP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
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
            Y = nDecomp*squeeze(sum(Y,1));
            if size(Y,1)~=1
                Y= Y.';
            end
        end
                function DataCompare(self,params,xdata,ydata)
            % DATACOMPARE: a function for comparing some data with a set of
            % parameters and the model, the inital condition for the model
            % is taken from the data.
            nDecomp = params.nDecomp;
            flipAngleProfiles = params.flipAngleProfiles;
            [decompFlipAngles,tmpM0] = self.decompSlice(ydata(:,1),nDecomp,flipAngleProfiles);
            tmpM0 = self.splitM0(tmpM0,params);
            for m = 1:nDecomp
                params.FaList = decompFlipAngles(:,m);
            [TRList,Mxy(m,:,:),~] = self.compile(tmpM0,params);
            end
            Mxy = nDecomp*squeeze(sum(Mxy,1));
            legendVals = cell(size(Mxy,1),1);
            figure
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

