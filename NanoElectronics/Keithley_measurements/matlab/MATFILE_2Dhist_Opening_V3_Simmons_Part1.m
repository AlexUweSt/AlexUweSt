% THIS SCRIPT DEPENDS ON DATA FROM MAT-FILES THAT WERE CREATED BY "MATFILE_OpeningClosing.m" !!!



R0 = 2000;
nbins = 200;
GG0trigger = 0.9; % trigger value, _last_ occurence defines position of dt=0
dTrigger = 0.5;  % delta of trigger value
dt_min = -20; % minimum time difference in seconds (x-axis limit)
dt_max = 1000; % maximum time difference in seconds (x-axis limit)
cmax = 200; % upper colormap limit in counts/bin, set to 0 for autolimit
a = -5; % 10^a = minimum G/G0 (y-axis limit)
b = 0;  % 10^b = maximum G/G0 (y-axis limit)
black_grid = 0; % (0/1) disable/enable black borders around bins
debug_plots = 0; % (0/1) disable/enable additional plots
DispText = 0; % (0/1) disable/enable plotting of info text



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
WF_Au = 5.1;
Rsample = cellfun(@(x) x - R0, R, 'UniformOutput',false); %subtract R0 from R
GG0 = cellfun(@(x) 1 ./ (x * G0), Rsample, 'UniformOutput',false); % calculate G/G0
x_Distance = cellfun(@(x) real( SimmonsModel( 'Distance', x, WF_Au ) ), GG0, 'UniformOutput',false); % calculate Distance using Simmons Model
%2do: remove GG0<0 instead of using real part!


%% find last occurence of GG0trigger (+- dTrigger) and filter data (dt_min ... dt_max)
TS_t0_idx = cell(1, FileCount);
TS_filtered = cell(1, FileCount);
GG0_filtered = cell(1, FileCount);
GG0_TS_x_filtered = cell(3, FileCount);
for i=1:FileCount
    GG0filter = GG0{1,i} >= GG0trigger-dTrigger & GG0{1,i} <= GG0trigger+dTrigger;
    if sum(GG0filter) > 0
        TS_t0_idx{i} = find(GG0filter,1,'last'); % find last occurence
        TStemp = TS{1,i}(:) - TS{1,i}(TS_t0_idx{i});
        GG0_filtered{i} = GG0{1,i}( TStemp >= dt_min & TStemp <= dt_max );
        TS_filtered{i} = TStemp( TStemp >= dt_min & TStemp <= dt_max );
        GG0_TS_x_filtered{1,i} = GG0_filtered{i};
        GG0_TS_x_filtered{2,i} = TS_filtered{i};
        GG0_TS_x_filtered{3,i} = x_Distance{1,i}( TStemp >= dt_min & TStemp <= dt_max );
    end
end;
clear GG0filter TStemp;


%% PLOT GG0 and x vs TS
if debug_plots == 1
    figure(); % GG0 vs TS
        set(0,'DefaultAxesLineStyleOrder',{'.-'}); %change from '-' to '.-' or '.' to see datapoints!
        %set(0,'DefaultAxesColorOrder',jet(FileCount));
        plot( GG0_TS_x_filtered{[2 1],:} );
        %plot curve number at last data point (to identify curves with strange behaviour)
        for i=1:FileCount
            if ~isempty(GG0_TS_x_filtered{1,i})
                text(GG0_TS_x_filtered{2,i}(end),GG0_TS_x_filtered{1,i}(end),num2str(i));
            end
        end
        
        a = -2; %y-lim
        b = 0; 
        ylim( [10^a 10^b] );
        xlabel('\fontsize{14}{\itt} / s'); ylabel('\fontsize{14}{\itG} / G_0');
    figure(); % x vs TS
        set(0,'DefaultAxesLineStyleOrder',{'.'}); %change from '-'  to '.-' or '.' to see datapoints!
        plot( GG0_TS_x_filtered{2:3,:} );
        xlabel('\fontsize{14}{\itt} / s'); ylabel('\fontsize{14}Electrode Distance / m');
end


%% [2mod!] Remove Curves where x < 0
Slope_GG0_TS_x_filtered = GG0_TS_x_filtered;
if debug_plots == 1
    figure(); % GG0 vs TS
        set(0,'DefaultAxesLineStyleOrder',{'.-'}); %change from '-' to '.-' or '.' to see datapoints!
        %set(0,'DefaultAxesColorOrder',jet(FileCount));
        plot( Slope_GG0_TS_x_filtered{[2 1],:} );
        xlabel('\fontsize{14}{\itt} / s'); ylabel('\fontsize{14}{\itG} / G_0');
    figure(); % x vs TS
        set(0,'DefaultAxesLineStyleOrder',{'.'}); %change from '-'  to '.-' or '.' to see datapoints!
        plot( Slope_GG0_TS_x_filtered{2:3,:} );
        xlabel('\fontsize{14}{\itt} / s'); ylabel('\fontsize{14}Electrode Distance / m');
end


%% Estimate the slope of x vs TS
t1 = 10; %=-1
t2 = 1000; % old: 40 4
t2_min_x = 3e-10; % NEW! minimum distance to be reached at t2 (else curve gets filtered)0
t2_max_x = 6e-10; % NEW! maximum distance to be reached at t2 (else curve gets filtered)9e-8
Slope_GG0_TS_x_filtered = GG0_TS_x_filtered;
SlopeFilter = cell(1, FileCount);
SlopeFilter_Idx_Xintercept = cell(1, FileCount);
Slope_x_shift = cell(1, FileCount);
for i=1:FileCount %25
    if ~isempty(Slope_GG0_TS_x_filtered{1,i}) && ~isempty(find(Slope_GG0_TS_x_filtered{2,i}<=t2,1,'last')) && ~isempty(find(Slope_GG0_TS_x_filtered{2,i}<=t2,1,'last'))
        if Slope_GG0_TS_x_filtered{3,i}( find(Slope_GG0_TS_x_filtered{2,i}<=t2,1,'last') ) >= t2_min_x && Slope_GG0_TS_x_filtered{3,i}( find(Slope_GG0_TS_x_filtered{2,i}<=t2,1,'last') ) <= t2_max_x
            SlopeFilter{i} = Slope_GG0_TS_x_filtered{2,i} < t1 | Slope_GG0_TS_x_filtered{2,i} > t2;
            Slope_GG0_TS_x_filtered{1,i}(SlopeFilter{i}) = [];
            Slope_GG0_TS_x_filtered{2,i}(SlopeFilter{i}) = [];
            Slope_GG0_TS_x_filtered{3,i}(SlopeFilter{i}) = [];
        else
            Slope_GG0_TS_x_filtered{1,i} = [];
            Slope_GG0_TS_x_filtered{2,i} = [];
            Slope_GG0_TS_x_filtered{3,i} = [];
        end
        if ~isempty(Slope_GG0_TS_x_filtered{1,i})
            % move all curves to same origin of x:
            Slope_x_shift{i} = Slope_GG0_TS_x_filtered{3,i}(1);
            Slope_GG0_TS_x_filtered{3,i} = Slope_GG0_TS_x_filtered{3,i} - Slope_x_shift{i};
        end
    end
end
for i=1:FileCount
    SlopeFilterBool(i) = ~isempty(SlopeFilter{1})+1;
end
if sum(SlopeFilterBool)==0, error('NO REMAINING CURVES FOR SLOPE ESTIMATION!'); end;
Slope_TSall = vertcat( Slope_GG0_TS_x_filtered{2,:} );
Slope_x_all = vertcat( Slope_GG0_TS_x_filtered{3,:} );
%SlopeFit = polyfit(Slope_TSall,Slope_x_all,1);
SlopeFit = fit(Slope_TSall,Slope_x_all,'poly1','Robust','Bisquare');
SlopeFitCoeff = coeffvalues(SlopeFit);
figure(); % x vs TS
    set(0,'DefaultAxesLineStyleOrder',{'.'}); %change from '-'  to '.-' or '.' to see datapoints!
    plot( Slope_GG0_TS_x_filtered{2:3,:} );
    %hold on; h = plot((t1:t2),polyval(SlopeFit,(t1:t2)),'k--','LineWidth',2); legend(h,['\delta x = ' num2str(SlopeFit(1)) '*t + ' num2str(SlopeFit(2))]);
    hold on; h = plot(SlopeFit); set(h, 'Color','k','LineStyle','--','LineWidth',2); legend(h,['\delta x = ' num2str(SlopeFitCoeff(1)) '*t + ' num2str(SlopeFitCoeff(2))]);
    xlabel('\fontsize{14}{\itt} / s'); ylabel('\fontsize{14}{\delta} Electrode Distance / m');




%% STOP HERE (for now...)
break;



%% Prepare data for x vs TS histograms
TSall = vertcat( GG0_TS_x_filtered{2,:} );
x_all = vertcat( GG0_TS_x_filtered{3,:} );
x_all( TSall(:) > dt_max ) = [];
TSall( TSall(:) > dt_max ) = [];


%% PLOT 2D Histogram: x vs TS
%if debug_plots == 1
    %Create plot
    figure();
        hist3( [TSall,x_all], [nbins nbins] );
        view(3);

    %Create 2D histogram
    figure();
        [N,C] = hist3([TSall,x_all],[nbins nbins]);
        imagesc(C{1,1}(:),C{1,2}(:),N');
        set(gca,'YDir','normal'); % revert y-axis back to normal [see http://www.mathworks.com/support/solutions/en/data/1-16MFFQ/ ]
        set(findobj(gca,'Type','axes'),'TickDir','out','XMinorTick','on','YMinorTick','on','PlotBoxAspectRatio',[1 1 1]); % change plot settings
        xlabel('\fontsize{14}{\itt} / s'); ylabel('\fontsize{14}Electrode Distance / m');
        colormap jet;
        caxis([0 25]);
        colorbar;
%end



%%% --- UNMODIFIED CODE AFTER THIS LINE ---

%% Prepare data for histograms
TSall = vertcat(TS_filtered{1,:});
GG0all = vertcat(GG0_filtered{1,:});
GG0all( TSall(:) > dt_max ) = [];
TSall( TSall(:) > dt_max ) = [];


%% PLOT 2D Opening-Histogram
if debug_plots == 1
    %Create plot
    figure();
    hist3( [TSall,GG0all], [nbins nbins] );
    view(3);

    %Create 2D histogram
    figure();
    [N,C] = hist3([TSall,GG0all],[nbins nbins]);
    imagesc(C{1,1}(:),C{1,2}(:),N');
    set(gca,'YDir','normal'); % revert y-axis back to normal [see http://www.mathworks.com/support/solutions/en/data/1-16MFFQ/ ]
    set(findobj(gca,'Type','axes'),'TickDir','out','XMinorTick','on','YMinorTick','on','PlotBoxAspectRatio',[1 1 1]); % change plot settings
    colormap jet;
    caxis([0 10]);
    colorbar;
end


%% Create semi-log 2D histogram
figure();
    bin_edges = { linspace(dt_min,dt_max,nbins+1), logspace(a,b,nbins+1) };
    %[N,C] = hist3([TSall,GG0all],'Edges',bin_edges); % bin-center (C) calculation of hist3 is wrong!
    N = hist3( [TSall,GG0all], 'Edges', bin_edges ); 
    %imagesc(C{1,1}(:),C{1,2}(:),N');
    %pcolor(C{1,1}(:),C{1,2}(:),N');
    pcolor( bin_edges{1}, bin_edges{2}, N' );
    %set(gca,'YDir','normal'); % revert y-axis back to normal [see http://www.mathworks.com/support/solutions/en/data/1-16MFFQ/ ]
    set(findobj(gca,'Type','axes'),'TickDir','out','XMinorTick','on','YMinorTick','on','PlotBoxAspectRatio',[1 1 1]); % change plot settings
    if black_grid == 0, shading flat; else shading faceted; end;
    colormap jet;
    if cmax > 0, caxis([0 cmax]); end;
    h=colorbar;
    %caxis( [caxis( [1 100] );] ); %log scaling (2do!)
    %set(h,'YScale','log'); %log scaling (2do!)
    ylim( [10^a 10^b] );
    xlim( [dt_min, dt_max] );
    set(gca,'yscale','log');
    xlabel('\fontsize{14}{\itt} / s');
    ylabel('\fontsize{14}{\itG} / G_0');
    %set(gca,'FontSize',12); % increase tick label font size


%% GET SAMPLE NAME AND ENVIRONMENT
temp = strsplit(MATINFO_Opening_FILES.Path, filesep);
SamplePath = cell(1,2*length(temp)-3);
SamplePath(1:2:end) = temp(1:end-1);
SamplePath(2:2:end-1) = {filesep};
SampleName = temp{end-1};
DateAndEnv = temp{end}; % Measurement starting time and sample environment
clear temp;
if DispText == 1, title([SampleName ' \\ ' DateAndEnv],'Color', [0 0 0.56],'FontSize',9,'Interpreter','none'); end;


%% SAVE FIGURE TO PNG FILE
% PNGfilePath = [SamplePath{:} '\' DateAndEnv ' - 2Dhist_Opening_V3.png'];
% % Delete old PNG file if it exists
% if exist(PNGfilePath,'file')
%     recycle('off');
%     delete(PNGfilePath);
% end
% FUNC_Save2PNG(PNGfilePath);
