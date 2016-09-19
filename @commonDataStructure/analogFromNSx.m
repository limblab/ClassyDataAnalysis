function analogFromNSx(cds)
%takes a cds handle and an NEVNSx object and populates the analog cell
%array of the cds. Does not return anything
    %establish lists for force, emg, lfp
    forceList = [find(~cellfun('isempty',strfind(lower(cds.NSxInfo.NSx_labels),'force_'))),...
                find(~cellfun('isempty',strfind(cds.NSxInfo.NSx_labels,'ForceHandle')))];
    lfpList=find(~cellfun('isempty',strfind(lower(cds.NSxInfo.NSx_labels),'elec')));
    emgList=find(~cellfun('isempty',strfind(lower(cds.NSxInfo.NSx_labels),'emg_')));
    %get lists of the analog data for each frequency & remove those we have already handled
    analogList=setxor(1:length(cds.NSxInfo.NSx_labels),emgList);
    analogList=setxor(analogList,forceList);
    analogList=setxor(analogList,lfpList);
    if ~isempty(analogList)
        %get our list of frequencies:
        frequencies=unique(cds.NSxInfo.NSx_sampling(analogList));
        if ~isempty(cds.analog)
            %get a list of the frequencies already in the cds so that we
            %can merge or add the fields from the cds to the new analog
            %data
            cdsFrequencies=zeros(1,length(cds.analog));
            for i=1:length(cds.analog)
                cdsFrequencies(i)=round(1/mode(diff(cds.analog{i}.t)));
            end
        else
            cdsFrequencies=[];
        end
        for i=1:length(frequencies)
            %find the channels that area actually at this frequency:
            subset=find(cds.NSxInfo.NSx_sampling(analogList)==frequencies(i));
            analogIdx=analogList(subset);
            %append data in the subset into a single matrix:
            a=[];
            for c=1:numel(subset)
                switch frequencies(i)
                    case 500
                    nsLabel='NS1';
                    case 1000
                        nsLabel='NS2';
                    case 2000
                        nsLabel='NS3';
                    case 10000
                        nsLabel='NS4';
                    case 30000
                        nsLabel='NS5';
                    otherwise
                        error('analogFromNSx:unexpectedFrequency',['this function is not set up to handle data with collection frequency: ',num2str(frequencies(i))])
                end
                numPts=numel(cds.(nsLabel).Data(cds.NSxInfo.NSx_idx(analogIdx(c)),:));
                a{c+1}=double(cds.(nsLabel).Data(cds.NSxInfo.NSx_idx(analogIdx(c)),:))';
            end
            %get a time vector t for this sampling frequency
            a{1} = ([0:length(a{2})-1]' / frequencies(i));
            %convert the matrix of data into a table:
            match=find(cdsFrequencies==frequencies(i),1);
            if ~isempty(match)
                temp=table(a{:},'VariableNames',[{'t'};cds.NSxInfo.NSx_labels(analogList(subset))]);
                temp.Properties.VariableDescriptions=[{'time'},repmat({'analog data'},1,numel(subset))];
                temp.Properties.Description=['table of analog data with collection frequency of: ', num2str(frequencies(i))];
                analogData{i}=mergeTables(a,cds.analog{match});
            else
                analogData{i}=table(a{:},'VariableNames',[{'t'},reshape(cds.NSxInfo.NSx_labels(analogList(subset)),1,numel(cds.NSxInfo.NSx_labels(analogList(subset))))]);
                analogData{i}.Properties.VariableDescriptions=[{'time'},repmat({'analog data'},1,numel(subset))];
                analogData{i}.Properties.Description=['table of analog data with collection frequency of: ', num2str(frequencies(i))];
            end
        end
    else
        analogData={};
    end
    %push the cell array of analog data tables into the cds:
    %cds.setField('analog',analogData)
    if ~isempty(analogData)
        if ~isempty(cds.analog)
            %find any frequencies that were in the cds but not in the new data and add them:
            for i=1:length(cds.analog)
                if isempty(find(frequencies==cdsFrequencies(i),1))
                    analogData=[analogData,cds.analog{i}];
                end
            end
        else
            set(cds,'analog',analogData)
        end
        
        evntData=loggingListenerEventData('analogFromNSx',[]);
        notify(cds,'ranOperation',evntData)
    end
end
