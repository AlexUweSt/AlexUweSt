% Create MAT files for opening/closing curves and I-V sweeps of selected measurement
%==========================================================================================
% MAT FILE: _R_Opening_[FILES]_YYYY-MM-DD_HHMM.mat
% MAT FILE: _R_Closing_[FILES]_YYYY-MM-DD_HHMM.mat
% MAT FILE: _IV_Sweeps_[FILES]_YYYY-MM-DD_HHMM.mat      230
%=========================================================================================C:\Users\Filip\Desktop\2016-10-11_14.59.54 - [MuMBox]_fluid(GToluol)=


clearvars -except MeasurementPath; clc;
Fmat = MCBJ_FUNCTIONS_MATFILE_with_motor_pos;

% von mir rausgenommen am 12.03.2015 --> Messpfad angeben
%if ~exist('MeasurementPath','var'), load('_config.mat'); end % get local configuration data
%if ~exist('MeasurementPath','var') || ~ischar(MeasurementPath) || ~exist(MeasurementPath,'dir'), MeasurementPath = pwd; end
%MeasurementPath = uigetdir(MeasurementPath, 'Measurement Directory');
%if ~MeasurementPath, error('Selection cancelled!'); end

%MeasurementPath = 'C:\Users\FilipK\Desktop\2018-11-14_09.53.50 - [MuMBox]_fluid(Mesytulene)_molecules(Coranulene)'
%MeasurementPath = 'C:\Users\FilipK\Desktop\neu\2018-09-24_14.42.56 - [MuMBox]_fluid(Toluol)_molecules(Salen_Mn)';
%MeasurementPath = 'D:\001_C60_measurements_03.2019\new from 23.03.2019\2019-03-19_12.51.32 - unknownEnv'
MeasurementPath = 'D:\Data work\PhD measurements\cleanup from desktop-random data\salen i dobri fajlovi\Kompilacija podataka za Nature\IRON\MCBJ_salen_Fe_04_04\2016-04-04_16.20.33 - [MuMBox]_fluid(Toluol)';
MeasurementPath = uigetdir(MeasurementPath,'Measurement Directory');

ts = now;

fprintf('Creating MAT file for opening curves...\n'); tic;
Fmat.MATFILE_Creator('Opening', MeasurementPath, ts);
fprintf('   FINISHED! (after %.3f seconds)\n', toc);

fprintf('Creating MAT file for closing curves...\n'); tic;
Fmat.MATFILE_Creator('Closing', MeasurementPath, ts);
fprintf('   FINISHED! (after %.3f seconds)\n', toc);

fprintf('Creating MAT file for I-V sweeps...\n'); tic;
Fmat.MATFILE_Creator('IVsweeps', MeasurementPath, ts);
fprintf('   FINISHED! (after %.3f seconds)\n', toc);

%% Motor position(modded) 
fprintf('Creating MAT file for motor position...\n'); tic;
Fmat.MATFILE_Creator('MotorPos', MeasurementPath, ts);
fprintf('   FINISHED! (after %.3f seconds)\n', toc);

clearvars;