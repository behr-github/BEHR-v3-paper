classdef misc_wrf_update_plots
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static = true)
        function plot_lonweight()
            UsGrid = GlobeGrid(0.05, 'domain', 'us');
            optimal_hour = 13.5 - UsGrid.GridLon/15;
            figure; pcolor(UsGrid.GridLon, UsGrid.GridLat, optimal_hour);
            colorbar;
            state_outlines('k');
        end
        
        function plot_daily_times()
            start_date = ask_date('Enter the start date');
            end_date = ask_date('Enter the end date');
            
            dvec = datenum(start_date):datenum(end_date);
            
            for d=1:numel(dvec)
                fprintf('Now on %s\n', datestr(dvec(d)));
                Data = load_behr_file(dvec(d), 'daily', 'us');
                
                if d==1
                    avg_hour = nan(size(Data(1).Grid.GridLon));
                    min_hour = nan(size(Data(1).Grid.GridLon));
                    max_hour = nan(size(Data(1).Grid.GridLon));
                    avg_count = zeros(size(Data(1).Grid.GridLon));
                end
                
                for a=1:numel(Data)
                    fprintf('  Swath %d of %d\n', a, numel(Data));
                    loncorns = Data(a).FoV75CornerLongitude;
                    latcorns = Data(a).FoV75CornerLatitude;
                    
                    % Get the time of the WRF file that was used. The last
                    % 8 characters of the file name should be HH-MM-SS or
                    % HH:MM:SS in UTC, so convert that to a numerical hour
                    % appropriately
                    wrf_file = Data(a).BEHRWRFFile;
                    wrf_hour = str2double(wrf_file(end-7:end-6));
                    wrf_minute = str2double(wrf_file(end-4:end-3));
                    wrf_utc = wrf_hour + wrf_minute / 60 * ones(size(Data(a).Longitude));
                    
                    % Grid the UTC hours and add to the running average
                    gridded_hours = cvm_generic_wrapper(loncorns, latcorns, wrf_utc, Data(1).Grid);
                    
                    notnans = ~isnan(gridded_hours);
                    avg_hour = nansum2(cat(3, avg_hour, gridded_hours),3);
                    avg_count = avg_count + notnans;
                    min_hour = min(min_hour, gridded_hours);
                    max_hour = max(max_hour, gridded_hours);
                end
            end
            
            % Finish up the running average
            avg_hour = avg_hour ./ avg_count;
            
            figure;
            pcolor(Data(1).Grid.GridLon, Data(1).Grid.GridLat, avg_hour);
            shading flat
            cb=colorbar;
            cb.Label.String = 'Average UTC hour';
            state_outlines('k');
            
            figure;
            pcolor(Data(1).Grid.GridLon, Data(1).Grid.GridLat, min_hour);
            shading flat
            cb=colorbar;
            cb.Label.String = 'Minimum UTC hour';
            state_outlines('k');
            
            figure;
            pcolor(Data(1).Grid.GridLon, Data(1).Grid.GridLat, max_hour);
            shading flat
            cb=colorbar;
            cb.Label.String = 'Max UTC hour';
            state_outlines('k');
        end
    end
    
end

