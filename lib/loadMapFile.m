function arrayMap=loadMapFile(fName)
    %function to open and parse a map file. Returns the map file as a table
    %with columns for the:
    %hardware channel   (1:96)
    %pin on the bank    (1:32)
    %bank               A,B,C or D
    %row                row of the electrode in the 10x10 grid
    %col                column of the electrode in the 10x10 grid
    %label              typically elec#
    arrayMap={};
    fid=fopen(fName,'r');
    inheader=true;
    while ~feof(fid)
        tline=fgets(fid);
        if isempty(deblank(tline))
            continue
        end
        %skip header lines:
        if inheader && ~isempty(strfind(tline,'Cerebus mapping'))
            %this is our last header line, set the flag false and
            %continue to the next line:
            inheader=false;
            continue
        elseif inheader
            continue
        end
        if strcmp(tline(1:2),'//')
            %this should be our column labels
            colLabels=textscan(tline(3:end),'%s')';%strsplit kind of works but leaves an empty string at the end that screws up the rest of processing
            colLabels=colLabels{1};
            rowNumCol=strcmp(colLabels,'row');
            colNumCol=strcmp(colLabels,'col');
            bankCol=strcmp(colLabels,'bank');
            pinCol=strcmp(colLabels,'elec');
            labelCol=strcmp(colLabels,'label');
            continue
        end
        %if we got to this point we are on an actual data line:
        tmp=textscan(tline,'%s');
        tmp=tmp{1};
        if numel(tmp)~=5
            error('unitsFromNEV:unexpectedMapFormat',['Was expecting to find exactly 5 columns in the mapfile, instead found: ',num2str(numel(tmp))])
        end
        if exist('colLabels','var')
            %we have a newer map file with actual column labels
            %use the labels to assign values
            rowNum=str2num(tmp{rowNumCol})+1;
            colNum=str2num(tmp{colNumCol})+1;
            pin=str2num(tmp{pinCol});
            bank=char(tmp{bankCol});
            label=char(tmp{labelCol});
        else
            %we have an older file with no labels
            warning('unitsFromNEV:oldMapfile','This mapfile does not have column labels. This is typical of older mapfiles. The mapfile data has been processed assuming the older Blackrock format, but this may cause issues for custom mapfiles or other oddities')
            rowNum=str2num(tmp{1})+1;
            colNum=str2num(tmp{2})+1;
            pin=str2num(tmp{4});
            bank=char(tmp{3});
            label=char(tmp{5});
        end
        switch bank
            case 'A'
                chan=pin;
            case 'B'
                chan=pin+32;
            case 'C'
                chan=pin+64;
            case 'D'
                chan=pin+96;
            otherwise
                error('unitsFromNEV:badBankLabel',['unitsFromNEV is not configured to handle arrays with bank label: ',bank])
        end
        arrayMap=[arrayMap;{chan,pin,rowNum,colNum,bank,label}];
    end
    fclose(fid);
    arrayMap=cell2table(arrayMap,'VariableNames',{'chan','pin','row','col','bank','label'});
end