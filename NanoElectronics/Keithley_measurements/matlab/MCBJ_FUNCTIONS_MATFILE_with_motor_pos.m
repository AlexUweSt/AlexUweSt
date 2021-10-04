function f = MCBJ_FUNCTIONS_MATFILE
% MCBJ: MAT file functions
%==========================================================================================
% Authors:  Matthias Wieser
% Version:  0.50 [2014-10-14]
%==========================================================================================

    f.countHeaderlines = @countHeaderlines;
    f.MATFILE_Creator = @MATFILE_Creator;
end



%% Count header lines in text file
function n = countHeaderlines(FilePath)
    n = 0;
    fid=fopen(FilePath, 'r');
    txtline = fgetl(fid);
    while txtline(1)=='#'
        n = n+1;
        txtline = fgetl(fid);
    end;
    fclose(fid);
end



%% MAT file creation
function MATFILE_Creator(CurveType, MeasurementPath, ts)
    recycle('off'); % Don't use recycle bin when deleting files
    global temp;
    isOpening = false;
    isClosing = false;
    isSweep = false;
    isMotor = false;
    if strcmp(CurveType, 'Opening')
        isOpening = true;
        FilePath = [MeasurementPath '\R_opening'];
        if ~exist(FilePath,'dir'), FilePath = [MeasurementPath '\opening']; end; % old folder name (LabVIEW MCBJ)
        if ~exist(FilePath,'dir'), fprintf(2, '   WARNING: No opening curves directory found!\n'); return; end;
    elseif strcmp(CurveType, 'Closing')
        isClosing = true;
        FilePath = [MeasurementPath '\R_closing'];
        if ~exist(FilePath,'dir'), FilePath = [MeasurementPath '\closing']; end; % old folder name (LabVIEW MCBJ)
        if ~exist(FilePath,'dir'), fprintf(2, '   WARNING: No closing curves directory found!\n'); return; end;
    elseif strcmp(CurveType, 'IVsweeps')
        isSweep = true;
        FilePath = [MeasurementPath '\IV-Sweeps'];
        if ~exist(FilePath,'dir'), FilePath = [MeasurementPath '\IV-Sweep']; end; % old folder name (LabVIEW MCBJ)
        if ~exist(FilePath,'dir'), fprintf(2, '   WARNING: No I-V sweeps directory found!\n'); return; end;
    elseif strcmp(CurveType, 'MotorPos')
        isMotor = true;
        FilePath = [MeasurementPath '\Motor'];
        if ~exist(FilePath,'dir'), FilePath = [MeasurementPath '\MotorPos']; end; % old folder name (LabVIEW MCBJ)
        if ~exist(FilePath,'dir'), fprintf(2, '   WARNING: No motor positions directory found!\n'); return; end;        
    else
        error('MATFILE_Creator:Import:UnknownCurveType', ['Unknown CurveType: ' CurveType]);
    end

    if (isOpening || isClosing)
        % Try to find GZ or TXT files
        FileSearchGZ = '*.txt.gz';
        FileSearchTXT = '*.txt';
        FileInfoGZ = dir([FilePath '\' FileSearchGZ]);
        FileCountGZ = length(FileInfoGZ);
        if FileCountGZ > 0
            fprintf(1, '   INFO: Found %i *.txt.gz files! Processing...\n', FileCountGZ);
            FileType = 'gz'; % 'gz' = file source: MCBJ ACS
        else
            fprintf(2, '   WARNING: No *.txt.gz files found! Now trying *.txt files...\n');
            FileInfo = dir([FilePath '\' FileSearchTXT]);
            FileCount = length(FileInfo);
            if FileCount > 0
                fprintf(1, '   INFO: Found %i *.txt files! Processing...\n', FileCount);
                FileType = 'txt'; % 'txt' = file source: LabVIEW MCBJ
            else
                error('MATFILE_Creator:Import:NoTXT', 'No *.txt files found!')
            end
        end
        % ==========================================================================================
        % Scan files created by MCBJ LabVIEW program:
        if strcmp(FileType, 'txt')
            FileHeader = cell(1,FileCount); %Preallocation
            V = cell(1,FileCount); %Preallocation
            I = cell(1,FileCount); %Preallocation
            R = cell(1,FileCount); %Preallocation
            TS = cell(1,FileCount); %Preallocation
            TS_K6430 = cell(1,FileCount); %Preallocation
            for i=1:FileCount
                [temp, ~, nheaderlines] = importdata( [FilePath '\' FileInfo(i).name], '\t' );
                if nheaderlines == 0
                    V{i} = temp(:,6); % WARNING! OLD FILE FORMAT!!! (# Motorposition / Speed (set) / Speed (actual) / Timestamp (Command sent) / Timestamp (VIR Data received) / Voltage (V) / Current (A) / Resistance (Ohm) / K6430-Timestamp (s) [adjusted] / Status word (24bit))
                    I{i} = temp(:,7);
                    R{i} = temp(:,8);
                    TS{i} = temp(:,5);
                    TS_K6430{i} = temp(:,9);
                else
                    FileHeader{i} = temp.textdata;
                    V{i} = temp.data(:,6); % WARNING! OLD FILE FORMAT!!! (# Motorposition / Speed (set) / Speed (actual) / Timestamp (Command sent) / Timestamp (VIR Data received) / Voltage (V) / Current (A) / Resistance (Ohm) / K6430-Timestamp (s) [adjusted] / Status word (24bit))
                    I{i} = temp.data(:,7);
                    R{i} = temp.data(:,8);
                    TS{i} = temp.data(:,5);
                    TS_K6430{i} = temp.data(:,9);
                end
            end
            clear temp;
        % ==========================================================================================
        % Scan files created by MCBJ ACS Python program:
        elseif strcmp(FileType, 'gz')
            GZtmpPath = [FilePath '\' '___gunzip_temp'];
            if ~exist(GZtmpPath,'dir')
                mkdir(GZtmpPath)
            else
                % This case happens when the folder has not been deleted during a previous function execution!
                error('MATFILE_Creator:Import', 'Temporary path already exists!')
            end
            FileCount = FileCountGZ;
            FileHeader = cell(1,FileCount); %Preallocation
            V = cell(1,FileCount); %Preallocation
            I = cell(1,FileCount); %Preallocation
            R = cell(1,FileCount); %Preallocation
            TS = cell(1,FileCount); %Preallocation
            TS_K6430 = cell(1,FileCount); %Preallocation
            for i=1:FileCountGZ
                filenames = gunzip([FilePath '\' FileInfoGZ(i).name], GZtmpPath);
                if size(filenames) ~= 1
                    % This case should never happen!
                    error('MATFILE_Creator:Import', 'Error after gz extraction! Unknown content!')
                end
                FileInfo(i,1) = dir(filenames{1}); % get infos from extracted *.txt files
                nheaderlines = countHeaderlines( [GZtmpPath '\' FileInfo(i).name] );
                %[temp, ~, nheaderlines] = importdata( [GZtmpPath '\' FileInfo(i).name], '\t', 21 ); %2do: remove "21" when ACS uses [space]s instead of [tab] in comment lines!
                %[temp, ~, nheaderlines] = importdata( [GZtmpPath '\' FileInfo(i).name], '\t', 24 ); %2do: remove "24" when ACS uses [space]s instead of [tab] in comment lines!
                [temp, ~, nheaderlines] = importdata( [GZtmpPath '\' FileInfo(i).name], '\t', nheaderlines );
                if nheaderlines ~= 0
                    FileHeader{i} = temp.textdata;
                    V{i} = temp.data(:,1); % FILE FORMAT: # Voltage (V) / Current (A) / Device-Timestamp (s) [not adjusted]
                    I{i} = temp.data(:,2);
                    TS_K6430{i} = temp.data(:,3);
                    % calculate additional data from file data:
                    R{i} = V{i} ./ I{i};
                    x = 0; %2do!!! (scrape header for "# DeviceTimeReset:		1378057146.019")
                    TS{i} = x + TS_K6430{i};
                else
                    % This case should never happen!
                    error('MATFILE_Creator:Import', 'The *.txt.gz file has no header information!')
                end
                % Delete temporary unzipped file
                delete(filenames{1});
            end
            clear temp;
            % Delete temporary folder and its files
            if exist(GZtmpPath,'dir')
                rmdir(GZtmpPath,'s');
            end
        end
        
    elseif (isSweep)
        % I-V sweeps: Try to find GZ or TXT files
        FileSearchGZ = '*.txt.gz';
        FileSearchTXT = '*.txt';
        FileInfoGZ = dir([FilePath '\' FileSearchGZ]);
        FileCountGZ = length(FileInfoGZ);
        if FileCountGZ > 0
            fprintf(1, '   INFO: Found %i *.txt.gz files! Processing...\n', FileCountGZ);
            FileType = 'gz'; % 'gz' = file source: MCBJ ACS
        else
            fprintf(2, '   WARNING: No *.txt.gz files found! Now trying *.txt files...\n');
            FileInfo = dir([FilePath '\' FileSearchTXT]);
            FileCount = length(FileInfo);
            if FileCount > 0
                fprintf(1, '   INFO: Found %i *.txt files! Processing...\n', FileCount);
                FileType = 'txt'; % 'txt' = file source: LabVIEW MCBJ
            else
                error('MATFILE_Creator:Import:NoTXT', 'No *.txt files found!')
            end
        end
       
        % ==========================================================================================
        % I-V sweeps: Scan files created by MCBJ LabVIEW program:
        if strcmp(FileType, 'txt')
            FileHeader = cell(1,FileCount); %Preallocation
            V = cell(1,FileCount); %Preallocation
            I = cell(1,FileCount); %Preallocation
            TS = cell(1,FileCount); %Preallocation
            TS_K6430 = cell(1,FileCount); %Preallocation
            for i=1:FileCount
                [temp, ~, nheaderlines] = importdata( [FilePath '\' FileInfo(i).name], '\t' );
                if nheaderlines == 0
                    V{i} = temp(:,1); % WARNING! OLD FILE FORMAT!!! (# Voltage (V) / Current (A))
                    I{i} = temp(:,2);
                else
                    FileHeader{i} = temp.textdata;
                    V{i} = temp.data(:,1); % WARNING! OLD FILE FORMAT!!! (# Voltage (V) / Current (A))
                    I{i} = temp.data(:,2);
                end
            end
            clear temp;
        % ==========================================================================================
        % I-V sweeps: Scan files created by MCBJ ACS Python program:
        elseif strcmp(FileType, 'gz')
            GZtmpPath = [FilePath '\' '___gunzip_temp'];
            if ~exist(GZtmpPath,'dir')
                mkdir(GZtmpPath)
            else
                % This case happens when the folder has not been deleted during a previous function execution!
                error('MATFILE_Creator:Import', 'Temporary path already exists!')
            end
            FileCount = FileCountGZ;
            FileHeader = cell(1,FileCount); %Preallocation
            V = cell(1,FileCount); %Preallocation
            I = cell(1,FileCount); %Preallocation
            TS = cell(1,FileCount); %Preallocation
            TS_K6430 = cell(1,FileCount); %Preallocation
            for i=1:FileCountGZ
                filenames = gunzip([FilePath '\' FileInfoGZ(i).name], GZtmpPath);
                if size(filenames) ~= 1
                    % This case should never happen!
                    error('MATFILE_Creator:Import', 'Error after gz extraction! Unknown content!')
                end
                FileInfo(i,1) = dir(filenames{1}); % get infos from extracted *.txt files
                nheaderlines = countHeaderlines( [GZtmpPath '\' FileInfo(i).name] );
                %[temp, ~, nheaderlines] = importdata( [GZtmpPath '\' FileInfo(i).name], '\t', 21 ); %2do: remove "21" when ACS uses [space]s instead of [tab] in comment lines!
                %[temp, ~, nheaderlines] = importdata( [GZtmpPath '\' FileInfo(i).name], '\t', 24 ); %2do: remove "24" when ACS uses [space]s instead of [tab] in comment lines!
                [temp, ~, nheaderlines] = importdata( [GZtmpPath '\' FileInfo(i).name], '\t', nheaderlines );
                if nheaderlines ~= 0
                    FileHeader{i} = temp.textdata;
                    V{i} = temp.data(:,1); % FILE FORMAT: # Voltage (V) / Current (A) / Device-Timestamp (s) [not adjusted]
                    I{i} = temp.data(:,2);
                    TS_K6430{i} = temp.data(:,3);
                    % calculate additional data from file data:
                    x = 0; %2do!!! (scrape header for "# DeviceTimeReset:		1378057146.019")
                    TS{i} = x + TS_K6430{i};
                else
                    % This case should never happen!
                    error('MATFILE_Creator:Import', 'The *.txt.gz file has no header information!')
                end
                % Delete temporary unzipped file
                delete(filenames{1});
            end
            clear temp;
            % Delete temporary folder and its files
            if exist(GZtmpPath,'dir')
                rmdir(GZtmpPath,'s');
            end
        end
        
         % ==========================================================================================
        %Motor position calculation
    elseif (isMotor)
        % I-V sweeps: Try to find GZ or TXT files
        FileSearchGZ = '*.txt.gz';
        FileSearchTXT = '*.txt';
        FileInfoGZ = dir([FilePath '\' FileSearchGZ]);
        FileCountGZ = length(FileInfoGZ);
        if FileCountGZ > 0
            fprintf(1, '   INFO: Found %i *.txt.gz files! Processing...\n', FileCountGZ);
            FileType = 'gz'; % 'gz' = file source: MCBJ ACS
        else
            fprintf(2, '   WARNING: No *.txt.gz files found! Now trying *.txt files...\n');
            FileInfo = dir([FilePath '\' FileSearchTXT]);
            FileCount = length(FileInfo);
            if FileCount > 0
                fprintf(1, '   INFO: Found %i *.txt files! Processing...\n', FileCount);
                FileType = 'txt'; % 'txt' = file source: LabVIEW MCBJ
            else
                error('MATFILE_Creator:Import:NoTXT', 'No *.txt files found!')
            end
        end
        % ==========================================================================================
        % Motor position: Scan files created by MCBJ ACS Python program:
        if strcmp(FileType, 'txt')
            FileHeader = cell(1,FileCount); %Preallocation
            V = cell(1,FileCount); %Preallocation
            I = cell(1,FileCount); %Preallocation
            TS = cell(1,FileCount); %Preallocation
            TS_K6430 = cell(1,FileCount); %Preallocation
            for i=1:FileCount
                [temp, ~, nheaderlines] = importdata( [FilePath '\' FileInfo(i).name], '\t' );
                if nheaderlines == 0
                    V{i} = temp(:,1); % WARNING! OLD FILE FORMAT!!! (# Voltage (V) / Current (A))
                    I{i} = temp(:,2);
                else
                    FileHeader{i} = temp.textdata;
                    V{i} = temp.data(:,1); % WARNING! OLD FILE FORMAT!!! (# Voltage (V) / Current (A))
                    I{i} = temp.data(:,2);
                end
            end
            clear temp;
        % ==========================================================================================
        % motor position: Scan files created by MCBJ ACS Python program:
        elseif strcmp(FileType, 'gz')
            GZtmpPath = [FilePath '\' '___gunzip_temp'];
            if ~exist(GZtmpPath,'dir')
                mkdir(GZtmpPath)
            else
                % This case happens when the folder has not been deleted during a previous function execution!
                error('MATFILE_Creator:Import', 'Temporary path already exists!')
            end
            FileCount = FileCountGZ;
            FileHeader = cell(1,FileCount); %Preallocation
            SS = cell(1,FileCount); %Preallocation
            SA = cell(1,FileCount); %Preallocation
            MP = cell(1,FileCount); %Preallocation
            TS = cell(1,FileCount); %Preallocation
            TS_K6430 = cell(1,FileCount); %Preallocation
            for i=1:FileCountGZ
                filenames = gunzip([FilePath '\' FileInfoGZ(i).name], GZtmpPath);
                if size(filenames) ~= 1
                    % This case should never happen!
                    error('MATFILE_Creator:Import', 'Error after gz extraction! Unknown content!')
                end
                FileInfo(i,1) = dir(filenames{1}); % get infos from extracted *.txt files
                nheaderlines = countHeaderlines( [GZtmpPath '\' FileInfo(i).name] );
                %[temp, ~, nheaderlines] = importdata( [GZtmpPath '\' FileInfo(i).name], '\t', 21 ); %2do: remove "21" when ACS uses [space]s instead of [tab] in comment lines!
                %[temp, ~, nheaderlines] = importdata( [GZtmpPath '\' FileInfo(i).name], '\t', 24 ); %2do: remove "24" when ACS uses [space]s instead of [tab] in comment lines!
                [temp, ~, nheaderlines] = importdata( [GZtmpPath '\' FileInfo(i).name], '\t', nheaderlines );
                if nheaderlines ~= 0
                    FileHeader{i} = temp.textdata;
                    TS_K6430{i} = temp.data(:,1); % FILE FORMAT: # Timestamp / Motorposition / Speed (set) / Speed (actual)
                    MP{i} = temp.data(:,2); 
                    SS{i} = temp.data(:,3);
                    SA{i} = temp.data(:,4);
                    % calculate additional data from file data:
                    x = 0; %2do!!! (scrape header for "# DeviceTimeReset:		1378057146.019")
                    %TS{i} = x + TS_K6430{i};
                    TS{i}=TS_K6430{i}-temp.data(1,1);
                else
                    % This case should never happen!
                    error('MATFILE_Creator:Import', 'The *.txt.gz file has no header information!')
                end
                % Delete temporary unzipped file
                delete(filenames{1});
            end
            clear temp;
            % Delete temporary folder and its files
            if exist(GZtmpPath,'dir')
                rmdir(GZtmpPath,'s');
            end
        end
        
    end

    % Prepare content for MAT file
    temp = struct( 'Type', CurveType, 'Version', '0.5', 'Date', datestr(ts,'yyyy-mm-dd HH:MM:SS'), 'Timestamp', ts, 'Path', MeasurementPath );
    eval(['MATINFO_' CurveType ' = temp;']); % rename variable
    clear temp;
    if (isOpening || isClosing)
        MATfileVars = {['MATINFO_' CurveType], 'FileInfoGZ', 'FileInfo', 'FileCount', 'FileHeader', 'V', 'I', 'R', 'TS', 'TS_K6430'};
        MATfileTypeTXT = ['_R_' CurveType '_[FILES]_'];
    elseif (isSweep)
        MATfileVars = {['MATINFO_' CurveType], 'FileInfoGZ', 'FileInfo', 'FileCount', 'FileHeader', 'V', 'I', 'TS', 'TS_K6430'};
        MATfileTypeTXT = '_IV_Sweeps_[FILES]_';
    elseif (isMotor)
        MATfileVars = {['MATINFO_' CurveType], 'FileInfoGZ', 'FileInfo', 'FileCount', 'FileHeader', 'MP', 'SS', 'SA', 'TS_K6430','TS'};
        MATfileTypeTXT = '_Motor_Pos_[FILES]_';
    end
    MATfilePath = [MeasurementPath '\' MATfileTypeTXT datestr(ts,'yyyy-mm-dd_HHMM') '.mat'];
    % Delete old MAT file if it exists
    if exist(MATfilePath,'file')
        delete(MATfilePath);
    end
    % Save selected variables to MAT file in v7.0 MAT-format (less overhead than v7.3, but no support for files >2GB)
    save(MATfilePath, MATfileVars{:}, '-v7');

end
