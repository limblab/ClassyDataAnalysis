function [varargout] = fitCovariates(binned,varargin)
    %this is a method function of the binnedData class and should be saved
    %in the @binnedData folder with the class definition and other methods
    %files
    %
    %bd.fitCovariates uses the configuration in the bd.covariateConfig field to compute
    %GLM weights for each unit and stores the result in the
    %bd.covariateData field

    % set up function as either helper function or setting function
    if( nargin == 0 )
        covarCfg = binned.covariateConfig;
    elseif( nargin == 1 )
        covarCfg = varargin{1};
    else
        error('fitCovariates:nargin','Too many input arguments to fitCovariates')
    end

    %get our list of units
    if isempty(covarCfg.units)
        %find all our units and make a cell array containing the whole list
        unitMask=~cellfun(@(x)isempty(strfind(x,'CH')),binned.data.Properties.VariableNames) & ~cellfun(@(x)isempty(strfind(x,'ID')),binned.data.Properties.VariableNames);
        uList=binned.data.Properties.VariableNames(unitMask);
    else
        %use the list the user supplied
        uList=covarCfg.units;
        if ~iscellstring(uList)
            error('fitCovariates:unitListNotCellString','the list of units in binnedData.covariateConfig.units must be a cell array of strings, where each string is the name of a unit column in binnedData.data')
        end
    end

    %get the mask for the rows of binned data to fit
    if isempty(covarCfg.windows)
        rowMask=true(size(binned.data.t));
    else
        rowMask=windows2mask(binned.data.t,covarCfg.windows);
    end

    %% set up parallel processing
%     opt=setUpParallelProcessing(covarCfg.useParallel);

    %% Set up parameters for bootstrap and GLM
    disp('starting covariate fitting')
    % set boot function by checking stats toolbox version number
    if(verLessThan('stats','10'))
        error('COMPUTE_TUNING requires Statistics Toolbox version 10.0(R2015a) or higher');
    end
    noiseModel=covarCfg.glmNoiseModel;%if you don't abstract the noise model into a variable, then bootstrp will create copies of the whole binned object at each iteration.
    %% set up the modelSpec string that contains the wilkinson notation describing the model to fit
    fullInput={'x+y','vx+vy','fx+fy','speed'};
    inputMask=[covarCfg.pos,covarCfg.vel,covarCfg.force,covarCfg.speed];
    if covarCfg.speed
        speedTable=table(sqrt(binned.data.vx(rowMask).^2+binned.data.vy(rowMask).^2),'VariableNames',{'speed'});
    end
    inputList=[];
    if covarCfg.pos
        inputList=[inputList,{'x','y'}];
    end
    if covarCfg.vel
        inputList=[inputList,{'vx','vy'}];
    end
    if covarCfg.force
        inputList=[inputList,{'fx','fy'}];
    end
    inputSpec=strjoin(fullInput(inputMask),'+');
    %% build the data structure we will use to store all the PDs before merging into one big table
    pdType={'pos','vel','force'};
    for i=1:numel(pdType)
        type=pdType{i};
        if(covarCfg.(pdType{i}))
            data.(type).allPDs=zeros(numel(uList),1);
            data.(type).allPDCI=zeros(numel(uList),2);
            data.(type).allModdepth=zeros(numel(uList),1);
            data.(type).allModdepthCI=zeros(numel(uList),2);
            data.(type).allIstuned=zeros(numel(uList),1);
        end
    end
    %loop through each unit and compute tuning:
    tic
    for i=1:numel(uList)
        fprintf([uList{i},':','getting data subset(ET=',num2str(toc),'s).'])
        %% set up a mask for the columns we will use for this unit
%         colMask=false(1,numel(binned.data.Properties.VariableNames));
%         for j=1:numel(inputList)
%             colMask=colMask | strcmp(inputList{j},binned.data.Properties.VariableNames);
%         end
%         %% get subset of the data we will use for fitting:
%         colMask=colMask | strcmp(uList{i},binned.data.Properties.VariableNames);
        colMask=list2tableMask(binned.data,[inputList,uList(i)]);
        %if you don't make a sub-table, then bootstrp will include a copy of the WHOLE binned.data table in the output for EVERY iteration
        if covarCfg.speed
            dataTable=[binned.data(rowMask,colMask),speedTable];
        else
            dataTable=binned.data(rowMask,colMask);
        end
        %% run GLM
        fprintf(['  Bootstrapping GLM PD computation(ET=',num2str(toc),'s).'])
        %bootstrap for firing rates to get output parameters
        modelSpec=[uList{i},'~',inputSpec];
        bootfunc = @(data) fitglme(data,modelSpec,'Distribution',noiseModel);
        try
            bootTuning = bootstrp(covarCfg.bootstrapReps,@(data) {bootfunc(data)}, dataTable);
            bootCoef = cell2mat(cellfun(@(x) x.Coefficients.Estimate',bootTuning,'uniformoutput',false));
            bootPValues=cell2mat(cellfun(@(x) x.Coefficients.pValue',bootTuning,'uniformoutput',false));
            %compute the full GLM so we can drop terms later and
            %compute term significance
            fullModel=fitglme(dataTable,modelSpec,'Distribution',noiseModel);

        catch ME
            warning('fitCovariates:errorFittingGLM',['failed to fit glm for unit', uList{i}])
            disp('failed with error:')
            disp(ME.identifier)
            disp(ME.message)
            disp(ME.stack)
            disp('inserting NaN values and continuing')
            for j=1:numel(pdType)
                if(covarCfg.(pdType{j}))
                    data.(pdType{j}).allPDs(i)=nan;
                    data.(pdType{j}).allPDCIs(i,:)=[nan nan];
                    data.(pdType{j}).allModdepth(i)=nan;
                    data.(pdType{j}).allModdepthCI(i,:)=[nan nan];
                    data.(pdType{j}).allIstuned(i)=nan;
                end
            end
        end
    end
    %now compose table for the full set of tuning data:
    fprintf(['  Inserting PD data into binned.pdTable(',num2str(toc),').'])
    %get our columns describing the units in the output table:
    for i=1:numel(uList)
        %get position of 'CH' in the name string:
        CHLoc=strfind(uList{i},'CH');
        IDLoc=strfind(uList{i},'ID');
        arrayName(i)={uList{i}(1:CHLoc-1)};
        chan(i)=str2num(uList{i}(CHLoc+2:IDLoc-1));
        ID(i)=str2num(uList{i}(IDLoc+2:end));
    end
    for i=1:numel(pdType)
        type=pdType{i};
        if(binned.pdConfig.(type))
            vNames={[type,'Dir'],[type,'DirCI'],[type,'Moddepth'],[type,'ModdepthCI'],[type,'IsTuned']};
            pdTable=[pdTable,table(data.(type).allPDs,data.(type).allPDCIs,data.(type).allModdepth,data.(type).allModdepthCI,logical(data.(type).allIstuned),'VariableNames',vNames)];
        end
    end

    %% Either return data or put into binnedData object and log event
    if( nargout == 0 )
        set(binned,'covariateData',covariateTable)
    elseif( nargout == 1 )
        % return bootstrapped stuff
    else
        error('fitCovariates:nargin','Too many input arguments to fitCovariates')
    end
    evntData=loggingListenerEventData('fitCovariates',covarCfg);
    notify(binned,'ranCovariateFit',evntData)
    disp('done computing GLM')
