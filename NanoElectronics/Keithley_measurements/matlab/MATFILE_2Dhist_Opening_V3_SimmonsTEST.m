% THIS SCRIPT DEPENDS ON DATA FROM MAT-FILES THAT WERE CREATED BY "MATFILE_OpeningClosing.m" !!!



R0 = 2000; %R0=10KOhm %old 2000

nbinsX =100; % bins: time (distance)100 def200
nbinsY =100; % bins: conductance100
GG0trigger = 0.9; % trigger value, _last_ occurrence defines position of dt=0
dTrigger = 0.6;  % delta of trigger value
dt_min = -80; % minimum time difference in seconds (x-axis limit) ::::default dt_min=-1
%---
dt_max = 40*140*10^(-1); % !!! maximum time difference in seconds (x-axis limit)
dt_max = 140;
dx_o = (2.2421e-11/2)*2; % !!! slope during opening (motor speed = 50, motor = 161:1 !!!) --> calculated from I_03 first opening mat file; P_36_Au: 1.5449e-12,     2.2e-11 
%dx_o = (4e-11)
%---
a = -6; % 10^a = minimum G/G0 (y-axis limit)-8
b = 1;  % 10^b = maximum G/G0 (y-axis limit)
cmax =40; % upper colormap limit in counts/bin, set to 0 for autolimit cmax =55
DistanceBarSize = 1; % size of distance bar in nm (position: bottom right)
histnorm = 1; % (0/1) 1 = normalize histogram! [using (max(max(N)) as 100%]
black_grid = 0; % (0/1) disable/enable black borders around bins
debug_plots = 0; % (0/1) disable/enable additional plots
DispText = 0; % (0/1) disable/enable plotting of info text
HistVideo = 0; % (0/1) disable/enable video creation (needs FFDShow 32bit! --> VFW-Config: Encoder: select FFV1!)


%load('D:\Files\0000_OWN_STUFF\STUDIUM\Promotion\Bruchmechanik\Auswertung\MCBJ_DATA\MCBJ-PC-01\MCBJ_U-Au-N0108_fluid(Toluol)\2013-04-23_16.21.01\_Opening-FILES_2013-10-29_1128.mat')
%for i=[1:80 100:198], R{1,i} = []; end;
%---
%load('D:\Files\0000_OWN_STUFF\STUDIUM\Promotion\Bruchmechanik\Auswertung\MCBJ_DATA\MCBJ-PC-01\MCBJ_U-Au-N0013_Toluol\2011-05-27_09.34.35\_Opening-FILES_2014-02-25_2325.mat')
%for i=[1:1 50:FileCount], R{1,i} = []; end;
%for i=[1:30 40:FileCount], R{1,i} = []; end;

%% plot absolute TS and TS_K6430
% WARNING: K6430 sometimes sends the same timestamp value as before! Therefore we use the absolute unix timestamp in this script!
if debug_plots == 1
    for i=1:FileCount
        first_TS_of_File(i) = TS{1,i}(1);
        first_TS_K6430_of_File(i) = TS_K6430{1,i}(1);
    end;
    figure();
        plot(first_TS_K6430_of_File,'r.');
        title('first TS\_K6430 value of each file');
        xlabel('i=1:FileCount');
        ylabel('absolute TS\_K6430');
    figure();
        plot(first_TS_of_File,'r.');
        title('first TS value of each file');
        xlabel('i=1:FileCount');
        ylabel('absolute TS');
end


%% calc Rsample and G/G0
% e = 1.602 176 565(35) x 10^-19 C      [2013-02]
% h = 6.626 069 57(29) x 10^-34 Js      [2013-02]
G0 = 2 * (1.602176565e-19)^2 / 6.62606957e-34;
Rsample = cellfun(@(x) x - R0, R, 'UniformOutput',false); %subtract R0 from R
GG0 = cellfun(@(x) 1 ./ (x * G0), Rsample, 'UniformOutput',false); % calculate G/G0

%% find last occurrence of GG0trigger (+- dTrigger) and filter data (dt_min ... dt_max)
TS_filtered = cell(1, FileCount);
GG0_filtered = cell(1, FileCount);
for i=1:FileCount
    GG0filter = GG0{1,i} >= GG0trigger-dTrigger & GG0{1,i} <= GG0trigger+dTrigger;
    if sum(GG0filter) > 0
        TS_dt0_idx = find(GG0filter,1,'last'); % find last occurrence
        TStemp = TS{1,i}(:) - TS{1,i}(TS_dt0_idx);
        GG0_filtered{i} = GG0{1,i}( TStemp >= dt_min & TStemp <= dt_max );
        TS_filtered{i} = TStemp( TStemp >= dt_min & TStemp <= dt_max );
    end
end;
clear GG0filter TS_dt0_idx TStemp;


%% Convert time to electrode distance
x_filtered = cell(1, FileCount);
for i=1:FileCount
    x_filtered{i} = TS_filtered{i} * dx_o;
end;

%% Prepare data for histograms
TSall = vertcat(TS_filtered{1,:});
xall = vertcat(x_filtered{1,:});
GG0all = vertcat(GG0_filtered{1,:});
GG0all( TSall(:) > dt_max ) = [];
xall( TSall(:) > dt_max ) = [];
TSall( TSall(:) > dt_max ) = [];


%% PLOT 2D Opening-Histogram
if debug_plots == 1
    %Create plot
    figure();
    %hist3( [TSall,GG0all], [nbins nbins] );
    hist3( [xall,GG0all], [nbinsX nbinsY] );
    view(3);

    %Create 2D histogram
    figure();
    %[N,C] = hist3([TSall,GG0all],[nbins nbins]);
    [N,C] = hist3([xall,GG0all],[nbinsX nbinsY]);
    imagesc(C{1,1}(:),C{1,2}(:),N');
    set(gca,'YDir','normal'); % revert y-axis back to normal [see http://www.mathworks.com/support/solutions/en/data/1-16MFFQ/ ]
    set(findobj(gca,'Type','axes'),'TickDir','out','XMinorTick','on','YMinorTick','on','PlotBoxAspectRatio',[1 1 1]); % change plot settings
    colormap jet;
    caxis([0 10]);
    colorbar;
end


%% Create semi-log 2D histogram
figure();
set(gcf, 'Position', get(0, 'Screensize')); %setfullscreen
    %bin_edges = { linspace(dt_min,dt_max,nbins+1), logspace(a,b,nbins+1) };
    %bin_edges = { linspace(min(xall),max(xall),nbins+1), logspace(a,b,nbins+1) };
    bin_edges = { linspace(min(xall),dt_max*dx_o,nbinsX+1), logspace(a,b,nbinsY+1) };
    %[N,C] = hist3([TSall,GG0all],'Edges',bin_edges); % bin-center (C) calculation of hist3 is wrong!
    %N = hist3( [TSall,GG0all], 'Edges', bin_edges );
    N = hist3( [xall,GG0all], 'Edges', bin_edges );
    if histnorm == 1
        N = N / (max(max(N)) / 100);
    end;
    %imagesc(C{1,1}(:),C{1,2}(:),N');
    %pcolor(C{1,1}(:),C{1,2}(:),N');
    pcolor( bin_edges{1}, bin_edges{2}, N' );
    %set(gca,'YDir','normal'); % revert y-axis back to normal [see http://www.mathworks.com/support/solutions/en/data/1-16MFFQ/ ]
    %set(findobj(gca,'Type','axes'),'TickDir','out','XMinorTick','off','YMinorTick','on','PlotBoxAspectRatio',[1 1 1]); % change plot settings
    set(findobj(gca,'Type','axes'),'TickDir','out','XMinorTick','off','YMinorTick','on'); % change plot settings
   % set(gca,'XTick',0); % display only 0 on x-axis
    set(gca,'box','on'); % remove the ticks on top and on right side %off
    set(gca,'XTickLabel',get(gca,'XTick')); % suppress exponential formatting in x-axis figure ticks
    set(gca,'XTickLabelMode','auto');
    if black_grid == 0, shading flat; else shading faceted; end;
   
 shading interp;  
    colormap jet;
    
    
% load('MyColormaps','farbe')
% set(gcf,'Colormap',farbe)
    
    if cmax > 0, caxis([0 cmax]); end;
    h = colorbar;
    %caxis( [caxis( [1 100] );] ); %log scaling (2do!)
    %set(h,'YScale','log'); %log scaling (2do!)
    ylim( [10^a 10^b] );
    %xlim( [dt_min, dt_max] );
    set(gca,'yscale','log');
    %xlabel('\fontsize{14}{\itt} / s');
    %xlabel('\fontsize{14}distance'); ylabel('\fontsize{14}{\itG} / G_0');
    xlabel('\fontsize{50}distance'); ylabel('\fontsize{50}{\itG} (G_0)'); % [increased font size for paper!]
    %set(gca,'FontSize',12); % increase tick label font size (x-,y- and cbar-axis)
    set(gca,'FontSize',40); % increase tick label font size (x-,y- and cbar-axis) [increased font size for paper!]16
    set(h,'FontSize',40); % decrease font size of colorbar16
    %cmax = 50;
    caxis([0 cmax]);
    %shading interp;
    set(gca,'LineWidth',2);
        if histnorm == 1
        ylbl = get(h,'YTickLabel');
        set( h,'YTickLabel', [ylbl repmat(' %',size(ylbl,1),1)] );
        clear ylbl;
    end;
%ylbl([0 10^(-1) 10^(-2)]);
    set(gca,'YTick',[1E-6 1E-5 1E-4 1E-3 1E-2 1E-1 1E-0]); %setcorrect ticks
   % pbaspect([1 1 1]) %aspect ratio
    title('\fontsize{60}177s - 100 Curves')
    set(gcf, 'Position', [265 70 1140 910]);
    %--- annotation:
   % axPos = get(gca,'Position'); % axPos = [xMin,yMin,xExtent,yExtent] = [left bottom width height]
   % xMinMax = xlim();
   % annotationW = axPos(3) * ( DistanceBarSize*1e-9 / (xMinMax(2) - xMinMax(1)) ) ;
   % annotation('textbox', [axPos(1)+axPos(3)-annotationW 0.01 annotationW 0.055], 'String',{[num2str(DistanceBarSize) ' nm']},'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','top','FitBoxToText','off','BackgroundColor',[0 0 0],'Color',[1 1 1]);
    clear axPos xMinMax annotationW;


%% GET SAMPLE NAME AND ENVIRONMENT
ctemp = strsplit(MATINFO_Opening_FILES.Path, filesep);
SamplePath = cell(1,2*length(temp)-3);
SamplePath(1:2:end) = temp(1:end-1);
SamplePath(2:2:end-1) = {filesep};
SampleName = temp{end-1};
DateAndEnv = temp{end}; % Measurement starting time and sample environment
clear temp;
if DispText == 1, title([SampleName ' \\ ' DateAndEnv],'Color', [0 0 0.56],'FontSize',9,'Interpreter','none'); end;


%% GET DEFAULT EXPORT PATH: Documents \ MATLAB
ExportPath = [ winqueryreg('HKEY_CURRENT_USER','Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders','Personal') '\MATLAB\'];


%% SAVE FIGURE TO PNG FILE
% PNGfilePath = [SamplePath{:} '\' DateAndEnv ' - 2Dhist_Opening_V3.png'];
% % Delete old PNG file if it exists
% if exist(PNGfilePath,'file')
%     recycle('off');
%     delete(PNGfilePath);
% end
% FUNC_Save2PNG(PNGfilePath);


%% HISTOGRAM VIDEO
if HistVideo == 1
    ScrSize = get(0,'ScreenSize');
    clear MovieFrames; MovieFrames(1:sum(~cellfun(@isempty,x_filtered(1,:)))) = struct('cdata',[],'colormap',[]); %Preallocation for improved speed
    FrameIdx = 1;
    FrameSize = [ScrSize(3)/4 0.1*ScrSize(4) ScrSize(3)/2 0.8*ScrSize(4)];
    % Generate histogram plots and save them as movie frames
    figure('name','[avi] 2D-x Histogram VIDEO','numbertitle','off','Position',FrameSize)
    for FileNr = 1:FileCount
        if ~isempty( x_filtered{1,FileNr} )
            xall = x_filtered{1,FileNr};
            GG0all = GG0_filtered{1,FileNr};
            bin_edges = { linspace(min(xall),dt_max*dx_o,nbinsX+1), logspace(a,b,nbinsY+1) };
            N = hist3( [xall,GG0all], 'Edges', bin_edges );
            if histnorm == 1
                N = N / (max(max(N)) / 100);
            end;
            pcolor( bin_edges{1}, bin_edges{2}, N' );
            set(findobj(gca,'Type','axes'),'TickDir','out','XMinorTick','off','YMinorTick','on'); % change plot settings
            set(gca,'XTick',0); % display only 0 on x-axis
            set(gca,'box','off'); % remove the ticks on top and on right side
            set(gca,'XTickLabel',get(gca,'XTick')); % suppress exponential formatting in x-axis figure ticks
            if black_grid == 0, shading flat; else shading faceted; end;
            colormap jet;
            if cmax > 0, caxis([0 cmax]); end;
            h = colorbar;
            if histnorm == 1
                ylbl = get(h,'YTickLabel');
                set( h,'YTickLabel', [ylbl repmat(' %',size(ylbl,1),1)] );
                clear ylbl;
            end;
            %caxis( [caxis( [1 100] );] ); %log scaling (2do!)
            %set(h,'YScale','log'); %log scaling (2do!)
            ylim( [10^a 10^b] );
            %xlim( [dt_min, dt_max] );
            set(gca,'yscale','log');
            %xlabel('\fontsize{14}{\itt} / s');
            xlabel('\fontsize{14}distance');
            ylabel('\fontsize{14}{\itG} / G_0');
            %set(gca,'FontSize',12); % increase tick label font size
            %--- annotation:
            axPos = get(gca,'Position'); % axPos = [xMin,yMin,xExtent,yExtent] = [left bottom width height]
            xMinMax = xlim();
            annotationW = axPos(3) * ( DistanceBarSize*1e-9 / (xMinMax(2) - xMinMax(1)) ) ;
            annotation('textbox', [axPos(1)+axPos(3)-annotationW 0.01 annotationW 0.055], 'String',{[num2str(DistanceBarSize) ' nm']},'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','top','FitBoxToText','off','BackgroundColor',[0 0 0],'Color',[1 1 1]);
            clear axPos xMinMax annotationW;
            %infotext = { ['\fontsize{12}\color{black}  CurveNr. ',num2str(FileNr)] ['  ' FileInfo(FileNr).date] };
            %TextObj = text(min(V{FileNr}), max(I{FileNr}), infotext, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
            set(gcf,'Position',FrameSize); % make sure the frame has still the same size
            MovieFrames(FrameIdx) = getframe(gcf);
            FrameIdx = FrameIdx + 1;
        end
    end
    close(gcf); %close figure window
    % Save movie to file
%%
    filename = [ datestr(now,'yyyy-mm-dd_HHMMSSFFF ') SampleName ' - ' DateAndEnv '.avi' ];
    %movie2avi(MovieFrames, [filedir '\' filename], 'compression', 'FFDS', 'quality', 100, 'fps', 10);
    % Encoding Settings:  ffds = FFDShow Encoder (has to be configured with e.g. VirtualDub before use!!!), i420 = Intel Indeo
    %MovieFrames( cellfun(@isempty,MovieFrames) ) = []; % remove empty cells (Encoder doesn't like 0x0 frames ;-))
    movie2avi(MovieFrames, [ExportPath '\' filename], 'compression', 'FFDS', 'fps', 10);
    system(['explorer.exe ' ExportPath]); %open destination directory in windows explorer
end
