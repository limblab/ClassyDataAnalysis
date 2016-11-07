function [varargout] = fitCovariates(binned)
    %this is a method function of the binnedData class and should be saved
    %in the @binnedData folder with the class definition and other methods
    %files
    %
    %bd.fitCovariates uses the configuration in the bd.covariateConfig field to compute
    %GLM weights for each unit and stores the result in the
    %bd.covariateData field

    % To use as helper function for PD calculation, make new binnedData
    % object and run fitCovariates() on it.
    covarCfg = binned.covariateConfig;

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
    
    % set verbose flag
    verbose = ( ~isfield(covarCfg,'verbose') || covarCfg.verbose );

    %% set up parallel processing
%     opt=setUpParallelProcessing(covarCfg.useParallel);

    %% Set up parameters for bootstrap and GLM
    if(verbose)
        disp('starting covariate fitting')
    end
    % set boot function by checking stats toolbox version number
    if(verLessThan('stats','10'))
        error('COMPUTE_TUNING requires Statistics Toolbox version 10.0(R2015a) or higher');
    end
    noiseModel=covarCfg.glmNoiseModel;%if you don't abstract the noise model into a variable, then bootstrp will create copies of the whole binned object at each iteration.
    %% set up the modelSpec string that contains the wilkinson notation describing the model to fit
    inputSpec=strjoin(covarCfg.covarFilter,'+');
    %% build the data structure we will use to store all the PDs before merging into one big table
    models = cell(numel(uList),1);
    
    %% loop through each unit and compute tuning:
    if(verbose)
        tic
    end
    for i=1:numel(uList)
        if(verbose)
            fprintf([uList{i},':','getting data subset(ET=',num2str(toc),'s).'])
        end
        %% set up a mask for the columns we will use for this unit
        colMask=list2tableMask(binned.data,[covarCfg.covarFilter,uList(i)]);
        %if you don't make a sub-table, then bootstrp will include a copy of the WHOLE binned.data table in the output for EVERY iteration
        dataTable=binned.data(rowMask,colMask);
        
        %% run GLM
        if(verbose)
            fprintf(['  Bootstrapping GLM PD computation(ET=',num2str(toc),'s).'])
        end
        %bootstrap for firing rates to get output parameters
        modelSpec=[uList{i},'~',inputSpec];
%         bootfunc = @(data) fitglm(data,modelSpec,'Distribution',noiseModel);
        try
            %compute the full GLM so we can drop terms later and
            %compute term significance
            fullModel=fitglm(dataTable,modelSpec,'Distribution',noiseModel);
            
            if(~isfield(binned.pdConfig,'bootstrap') || binned.pdConfig.bootstrap)
%                 bootTuning = bootstrp(covarCfg.bootstrapReps,@(data) {bootfunc(data)}, dataTable);
%                 bootCoef = cell2mat(cellfun(@(x) x.Coefficients.Estimate',bootTuning,'uniformoutput',false));
%                 bootPValues=cell2mat(cellfun(@(x) x.Coefficients.pValue',bootTuning,'uniformoutput',false));
            else
%                 bootTuning = {fullModel};
%                 bootCoef = fullModel.Coefficients.Estimate';
%                 bootPValues = fullModel.Coefficients.pValue';
            end
            
            if(isfield(covarCfg,'outputModelClass') && covarCfg.outputModelClass)
                models{i} = fullModel;
            end
            

        catch ME
            warning('fitCovariates:errorFittingGLM',['failed to fit glm for unit', uList{i}])
            disp('failed with error:')
            disp(ME.identifier)
            disp(ME.message)
            disp(ME.stack)
            disp('inserting NaN values and continuing')
            if(isfield(covarCfg,'outputModelClass') && covarCfg.outputModelClass)
                models{i} = NaN;
            end
        end
    end
    %now compose table for the full set of tuning data:
    if(verbose)
        fprintf(['  Inserting GLM data into binned.covariateData(',num2str(toc),').'])
    end
    %get our columns describing the units in the output table:
    for i=1:numel(uList)
        %get position of 'CH' in the name string:
        CHLoc=strfind(uList{i},'CH');
        IDLoc=strfind(uList{i},'ID');
        arrayName(i)={uList{i}(1:CHLoc-1)};
        chan(i)=str2num(uList{i}(CHLoc+2:IDLoc-1));
        ID(i)=str2num(uList{i}(IDLoc+2:end));
    end
    covariateData=table(arrayName',chan',ID',models,'VariableNames',{'array','chan','ID','model'});
    % Add in covariate coefficients later

    %% Either return data or put into binnedData object and log event
    set(binned,'covariateData',covariateData)
    evntData=loggingListenerEventData('fitCovariates',covarCfg);
    notify(binned,'ranCovariateFit',evntData)
    if(verbose)
        disp('done computing GLM')
    end
