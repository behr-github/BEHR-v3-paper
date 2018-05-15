classdef misc_alb_plots
    %MISC_ALB_PLOTS Miscellaneous plots for checking the implementation of
    %the BRDF albedo in BEHR
    
    properties(Constant)
        mydir = fileparts(mfilename('fullpath'));
        black_sky_dir = '/Volumes/share-sat/SAT/BEHR/AlbedoTestBRDF/BlackSky';
        brdf_dir = '/Volumes/share-sat/SAT/BEHR/AlbedoTestBRDF/BRDF';
        brdf_flip_dir = '/Volumes/share-sat/SAT/BEHR/AlbedoTestBRDF/BRDF_tomradRAA';
        workspace_dir = fullfile(behr_paths.behr_core, 'Workspaces', 'BRDF Albedo');
        
        mcd43c1_dir = '/Volumes/share-sat/SAT/MODIS/MCD43C1';
        mcd43c3_dir = '/Volumes/share-sat/SAT/MODIS/MCD43C3/v006';
        scia_modis_dir = '/Users/Josh/Documents/Fortran/Sciatran/Utils/LUTGen/Tests/modis';
    end
    
    methods(Static = true)
        function make_black_sky_behr
            % Make the OMI_SP and OMI_BEHR files using NASA SPv3 but the
            % old MCD43C3 albedo.
            data_dir = misc_alb_plots.black_sky_dir;
            G = GitChecker;
            G.Strict = true;
            G.addCommitRange(behr_repo_dir, 'fb5f363', 'fb5f363'); % commit where I'd updated to SPv3 but not changed the albedo, then cherry picked the corner fixes
            G.checkState();
            
            read_omno2_v_aug2012('start', '2005-04-01', 'end', '2005-09-30', 'sp_mat_dir', data_dir);
            BEHR_main('start', '2005-04-01', 'end', '2005-09-30', 'sp_mat_dir', data_dir, 'behr_mat_dir', data_dir);
        end
        
        function make_brdf_behr
            % Make the OMI_SP and OMI_BEHR files using NASA SPv3 and the
            % new MCD43C1 BRDF.
            data_dir = misc_alb_plots.brdf_dir;
            G = GitChecker;
            G.Strict = true;
            G.addCommitRange(behr_repo_dir, 'f88f57c1', 'f88f57c1'); % commit where I'd finished merging the BRDF branch
            G.checkState();
            
            read_omno2_v_aug2012('start', '2005-04-01', 'end', '2005-09-30', 'sp_mat_dir', data_dir);
            BEHR_main('start', '2005-04-01', 'end', '2005-09-30', 'sp_mat_dir', data_dir, 'behr_mat_dir', data_dir);
        end
        
        function make_soas_behr_bs
            % Make the OMI_SP and OMI_BEHR files using NASA SPv3 and the
            % new MCD43C1 BRDF for the 
            G = GitChecker;
            G.Strict = true;
            G.addCommitRange(behr_repo_dir, 'fb5f363', 'fb5f363'); % the commit where I'd fixed some of the corner calculations, but not implemented the BRDF albedo
            G.checkState();
            
            data_dir = misc_alb_plots.black_sky_dir;
            read_omno2_v_aug2012('start', '2013-05-31', 'end', '2013-07-10', 'sp_mat_dir', data_dir, 'modis_mcd43_dir', '/Volumes/share-sat/SAT/MODIS/MCD43C3');
            BEHR_main('start', '2013-05-031', 'end', '2013-07-10', 'sp_mat_dir', data_dir, 'behr_mat_dir', data_dir);
        end
        
        function make_soas_behr_brdf
            % Make the OMI_SP and OMI_BEHR files using NASA SPv3 and the
            % new MCD43C1 BRDF for the 
            G = GitChecker;
            G.Strict = true;
            G.addReqCommits(behr_repo_dir, '86b60e9'); % commit where I'd made the BRDF calculation use the correct RAA definition (0 = backscattering)
            G.checkState();
            
            data_dir = misc_alb_plots.brdf_dir;
            read_omno2_v_aug2012('start', '2013-05-31', 'end', '2013-07-10', 'sp_mat_dir', data_dir, 'modis_mcd43_dir', '/Volumes/share-sat/SAT/MODIS/MCD43C1');
            BEHR_main('start', '2013-05-031', 'end', '2013-07-10', 'sp_mat_dir', data_dir, 'behr_mat_dir', data_dir);
        end
        
        function make_behr_avgs
            save_dir = misc_alb_plots.workspace_dir;
            sdate = '2005-04-01';
            edate = '2005-09-30';
            lonlim = [-125 -65];
            latlim = [25 50];
            loaddir = misc_alb_plots.black_sky_dir;
            [~,NO2_GRID, LONGRID, LATGRID] = no2_column_map_2014(sdate, edate, lonlim, latlim, 'behrdir', loaddir, 'rows', [0 29], 'makefig', false);
            save(fullfile(save_dir, 'BlackSkyLeft.mat'), 'NO2_GRID', 'LONGRID', 'LATGRID');
            [~,NO2_GRID, LONGRID, LATGRID] = no2_column_map_2014(sdate, edate, lonlim, latlim, 'behrdir', loaddir, 'rows', [30 59], 'makefig', false);
            save(fullfile(save_dir, 'BlackSkyRight.mat'), 'NO2_GRID', 'LONGRID', 'LATGRID');
            
            loaddir = misc_alb_plots.brdf_dir;
            [~,NO2_GRID, LONGRID, LATGRID] = no2_column_map_2014(sdate, edate, lonlim, latlim, 'behrdir', loaddir, 'rows', [0 29], 'makefig', false);
            save(fullfile(save_dir, 'BRDFLeft.mat'), 'NO2_GRID', 'LONGRID', 'LATGRID');
            [~,NO2_GRID, LONGRID, LATGRID] = no2_column_map_2014(sdate, edate, lonlim, latlim, 'behrdir', loaddir, 'rows', [30 59], 'makefig', false);
            save(fullfile(save_dir, 'BRDFRight.mat'), 'NO2_GRID', 'LONGRID', 'LATGRID');
        end
        
        function test_raa_definition
            % Typically, the relative azimuth angle is defined as the
            % absolute value of the difference between the solar and
            % viewing azimuth angles; it is usually constrained to be
            % between 0 and 180 deg. This means that backscattering should
            % be defined as RAA = 0 and forward scattering as RAA = 180,
            % especially since Roujean et al. 1992 (which sets up the
            % kernels used in the MODIS BRDF) says on p. 20457 (last
            % paragraph first column) "...the backscattering
            % direction...phi = 0" (phi is represents the RAA).
            %
            % However, in the TOMRAD lookup table, the definition of RAA is
            % flipped so that RAA = 180 is backscattering (I believe), so
            % in BEHR our usual definition of RAA matches this. Therefore,
            % I want to test that flipping the RAA back is the correct
            % definition for the BRDF. I'm hypothesizing that if we only
            % average elements in the backscattering direction (RAA_BEHR >
            % 90) that the BRDF-derived albedo should be greater when the
            % definition passed to the BRDF kernels is s.t. RAA = 0 is
            % backscattering.
            if ask_yn('Generate the albedo averages?')
                [back_eq_0, lon, lat] = manual_average(misc_alb_plots.brdf_dir, 'MODISAlbedo', @(raa) raa > 90);
                back_eq_180 = manual_average(misc_alb_plots.brdf_flip_dir, 'MODISAlbedo', @(raa) raa > 90);
                save('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/BRDF Albedo/BRDF Implementation/raa_definition.mat', 'back_eq_0', 'back_eq_180', 'lon', 'lat');
            else
                D = load('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/BRDF Albedo/BRDF Implementation/raa_definition.mat');
                back_eq_0 = D.back_eq_0;
                back_eq_180 = D.back_eq_180;
                lon = D.lon;
                lat = D.lat;
            end
            
            figure; 
            pcolor(lon,lat,back_eq_0);
            shading flat
            colorbar
            state_outlines('k','not','ak','hi');
            title('Backscatter phi = 0')
            
            figure; 
            pcolor(lon,lat,back_eq_180);
            shading flat
            colorbar
            state_outlines('k','not','ak','hi');
            title('Backscatter phi = 180')
            
            figure; 
            pcolor(lon,lat,back_eq_0-back_eq_180);
            shading flat
            colorbar
            state_outlines('k','not','ak','hi');
            title('Backscatter phi = 0 - phi = 180')
        end
        
        function soas_scatterplots
            % Makes scatter plots of BEHR with BRDF and black sky albedo
            % vs. aircraft profiles from SOAS
            
            G = GitChecker;
            G.Strict = true;
            G.addReqCommits(no2_prof_repo, 'f52e3f12');
            G.checkState();
            
            common_opts = {'campaign', 'soas',...
                           'profnums', 'all',...
                           'profile_input', 'ranges',...
                           'merge_dir', '',...
                           'behr_prefix', 'OMI_BEHR_v2-1C_',...
                           'startdate', '',...
                           'enddate', '',...
                           'starttime', '12:00',...
                           'endtime', '15:00',...
                           'timezone', 'auto',...
                           'no2field', 'cl',...
                           'no2conv', 1e-9,...
                           'altfield', '',...
                           'radarfield', '',...
                           'presfield', '',...
                           'tempfield', '',...
                           'minheight', 0,...
                           'numBLpts', 20,...
                           'minagl', 0.5,...
                           'useground', 0,...
                           'cloudtype', 'omi',...
                           'cloudfrac', 0.2,...
                           'rowanomaly', 'XTrackFlags',...
                           'behrfield', 'BEHRColumnAmountNO2Trop',...
                           'debug', 0,...
                           'clean', 1};
            
            [prof_lon, prof_lat, ~, brdf_behr_no2, brdf_air_no2, brdf_db] = Run_Spiral_Verification('behr_dir', misc_alb_plots.brdf_dir, common_opts{:});
            [~, ~, ~, bs_behr_no2, bs_air_no2, bs_db] = Run_Spiral_Verification('behr_dir', misc_alb_plots.black_sky_dir, common_opts{:});
            
            % Scatter plot
            [brdf_fit_x, brdf_fit_y, brdf_legstr] = calc_fit_line(brdf_air_no2, brdf_behr_no2, 'regression', 'rma', 'xcoord', 0:2e15:2e16);
            [bs_fit_x, bs_fit_y, bs_legstr] = calc_fit_line(bs_air_no2, bs_behr_no2, 'regression', 'rma', 'xcoord', 0:2e15:2e16);
            l = gobjects(4,1);
            figure;
            l(1) = plot(bs_air_no2, bs_behr_no2, 'bo', 'linewidth', 2);
            hold on
            l(2) = plot(bs_fit_x, bs_fit_y, 'b--', 'linewidth', 2);
            
            l(3) = plot(brdf_air_no2, brdf_behr_no2, 'ro', 'linewidth', 2);
            l(4) = plot(brdf_fit_x, brdf_fit_y, 'r--', 'linewidth', 2);
            
            xylims;
            xlabel('Aircraft NO_2');
            ylabel('BEHR NO_2');
            set(gca,'fontsize',16);
            legend(l, {'Black sky', bs_legstr, 'BRDF', brdf_legstr});
            
            % Scatter plot with lines through origin
            [brdf_fit_x, brdf_fit_y, brdf_legstr] = calc_fit_line(brdf_air_no2, brdf_behr_no2, 'regression', 'orth-origin', 'xcoord', 0:2e15:2e16);
            [bs_fit_x, bs_fit_y, bs_legstr] = calc_fit_line(bs_air_no2, bs_behr_no2, 'regression', 'orth-origin', 'xcoord', 0:2e15:2e16);
            l = gobjects(4,1);
            figure;
            l(1) = plot(bs_air_no2, bs_behr_no2, 'bo', 'linewidth', 2);
            hold on
            l(2) = plot(bs_fit_x, bs_fit_y, 'b--', 'linewidth', 2);
            
            l(3) = plot(brdf_air_no2, brdf_behr_no2, 'ro', 'linewidth', 2);
            l(4) = plot(brdf_fit_x, brdf_fit_y, 'r--', 'linewidth', 2);
            
            xylims;
            xlabel('Aircraft NO_2');
            ylabel('BEHR NO_2');
            set(gca,'fontsize',16);
            legend(l, {'Black sky', bs_legstr, 'BRDF', brdf_legstr});
            
            % Map plots
            xl = [min(prof_lon) - 2, max(prof_lon) + 2];
            yl = [min(prof_lat) - 2, max(prof_lat) + 2];
            
            figure; scatter(prof_lon, prof_lat, [], brdf_air_no2 - brdf_behr_no2, 'filled');
            cb = colorbar; cb.Label.String = 'BRDF VCD (Aircraft - BEHR)';
            colormap(jet); 
            state_outlines('k'); xlim(xl); ylim(yl); 
            
            figure; scatter(prof_lon, prof_lat, [], brdf_air_no2, 'filled');
            cb = colorbar; cb.Label.String = 'BRDF VCD (Aircraft)';
            colormap(jet); 
            state_outlines('k'); xlim(xl); ylim(yl);
            
            figure; scatter(prof_lon, prof_lat, [], brdf_behr_no2, 'filled');
            cb = colorbar; cb.Label.String = 'BRDF VCD (BEHR)';
            colormap(jet); 
            state_outlines('k'); xlim(xl); ylim(yl);
            
            figure; scatter(prof_lon, prof_lat, [], bs_behr_no2, 'filled');
            cb = colorbar; cb.Label.String = 'BRDF VCD (BEHR)';
            colormap(jet); 
            state_outlines('k'); xlim(xl); ylim(yl);
        end
        
        function varargout = ler_vs_simple_brdf(colorby, aerosol_type, DEBUG_LEVEL)
            allowed_colorbys = {'Longitude', 'Latitude', 'Site', 'SiteType', 'Month', 'SZA', 'VZA', 'RAA', 'None'};
            allowed_aer_types = {'continental', 'maritime'};
            
            if ~exist('colorby','var')
                colorby = ask_multichoice('Color the scatter plot by a variable?', allowed_colorbys, 'list', true);
            elseif ~ismember(colorby, allowed_colorbys)
                error('misc_alb_plots:bad_input', 'COLORBY must be one of %s', strjoin(allowed_colorbys, ', '));
            end
            if ~exist('aerosol_type', 'var')
                aerosol_type = ask_multichoice('Which type of aerosols for the bottom layer?', allowed_aer_types, 'list', true, 'default', 'continental');
            elseif ~ismember(aerosol_type, allowed_aer_types)
                error('misc_alb_plots:bad_input', 'AEROSOL_TYPE must be one of %s', strjoin(allowed_aer_types, ', '));
            end
            if ~exist('DEBUG_LEVEL', 'var')
                DEBUG_LEVEL = 1;
            end
            
            % Get the site output files
            aerosol_subdir = sprintf('aerosols_%s_mcd43d_more_szas', aerosol_type);
            locs_info = ncinfo(fullfile(misc_alb_plots.scia_modis_dir, 'loc_coeffs_D_new.nc'));
            loc_names = cellstr(ncread(locs_info.Filename, 'ShortName'));
            loc_coeffs = ncread(locs_info.Filename, 'ModisCoefficients');
            loc_months = double(ncread(locs_info.Filename, 'Months'));
                        
            for a=1:numel(loc_names)
                if DEBUG_LEVEL > 0
                    fprintf('Working on %s (%d of %d)\n', loc_names{a}, a, numel(loc_names));
                end
                loc_file = ncinfo(fullfile(misc_alb_plots.scia_modis_dir, aerosol_subdir, sprintf('%s-LER.nc', loc_names{a})));
                loc_lers = ncread(loc_file.Filename, 'LER');
                loc_lers = permute(loc_lers, [4,3,2,1]);
                
                if a==1
                    szas = ncread(loc_file.Filename, 'sza');
                    vzas = ncread(loc_file.Filename, 'vza');
                    raas = 180 - ncread(loc_file.Filename, 'raa'); % SCIATRAN RAA definition is backwards from MODIS
                    geometries = combvec(szas', vzas', raas');
                    modis_albs = nan(numel(loc_names), numel(loc_months), size(geometries,2));
                    modis_lers = nan(size(modis_albs));
                    color_vals = nan(size(modis_albs));
                end
                
                
                for b=1:numel(loc_months)
                    this_coeffs = loc_coeffs(:,3,b,a); % always use band 3
                    modis_albs(a,b,:) = modis_brdf_alb(this_coeffs(1), this_coeffs(2), this_coeffs(3), geometries(1,:), geometries(2,:), geometries(3,:));
                    modis_lers(a,b,:) = reshape(loc_lers(:,:,:,b),[],1);
                    switch lower(colorby)
                        case 'none'
                            color_vals(a,b,:) = 0;
                        case 'site'
                            color_vals(a,b,:) = a;
                        case 'sitetype'
                            site_tmp = find(strcmp(ncreadatt(loc_file.Filename, '/', 'sitetype'), {'Cities','PowerPlants','RuralAreas'}));
                            if isempty(site_tmp)
                                error('misc_alb_plots:site_type', 'Cannot identify site type for %s', loc_names{a})
                            end
                            color_vals(a,b,:) = site_tmp;
                        case 'month'
                            color_vals(a,b,:) = loc_months(b);
                        case 'sza'
                            color_vals(a,b,:) = geometries(1,:);
                        case 'vza'
                            color_vals(a,b,:) = geometries(2,:);
                        case 'raa'
                            color_vals(a,b,:) = geometries(3,:);
                        case {'longitude', 'latitude'}
                            color_vals(a,b,:) = ncreadatt(loc_file.Filename, '/', lower(colorby));
                        otherwise
                            error('misc_alb_plots:not_implemented', 'Color by %s not implemented', colorby);
                    end
                end
            end
            
            % Print the slopes and their standard deviations by bin
            switch lower(colorby)
                case 'none'
                    bin_edges = [-1 1];
                case 'site'
                    bin_edges = 1:(numel(loc_names)+1);
                case 'sitetype'
                    bin_edges = 1:4;
                case 'month'
                    bin_edges = 1:13;
                case 'sza'
                    bin_edges = 0:10:100;
                case 'vza'
                    bin_edges = 0:10:100;
                case 'raa'
                    bin_edges = 0:10:180;
                case 'longitude'
                    bin_edges = -125:5:-65;
                case 'latitude'
                    bin_edges = 25:5:50;
                otherwise
                    error('misc_alb_plots:not_implemented', 'Binning %s not implemented', colorby);
            end
            
            [alb_bins, bin_centers] = bin_data(color_vals, modis_albs(:), bin_edges);
            ler_bins = bin_data(color_vals, modis_lers(:), bin_edges);
            
            for i_bin = 1:numel(bin_centers)
                [~,~,~,fit] = calc_fit_line(alb_bins{i_bin}, ler_bins{i_bin}, 'regression', 'rma');
                fprintf('Slope for %s bin centered on %f = %f +/- %f\n', colorby, bin_centers(i_bin), fit.P(1), fit.StdDevM);
            end
            
            f1 = figure;
            scatter(modis_albs(:), modis_lers(:), [], color_vals(:));
            if ~strcmpi(colorby, 'None')
                cb = colorbar;
                cb.Label.String = colorby;
            end
            xlabel('BRDF');
            ylabel('LER');
            set(gca,'fontsize',16);
            xylims();
            plot_fit_line(gca,[],'regression','rma');
            f1.Children(1).Location = 'northwest';
            
            f4 = figure;
            scatter(modis_albs(:)/pi, modis_lers(:), [], color_vals(:));
            if ~strcmpi(colorby, 'None')
                cb = colorbar;
                cb.Label.String = colorby;
            end
            xlabel('BRF');
            ylabel('LER');
            set(gca,'fontsize',16);
            xylims();
            plot_fit_line(gca,[],'regression','rma');
            f1.Children(1).Location = 'northwest';
            
            perdel = reldiff(modis_lers(:), modis_albs(:))*100;
            f2 = figure; boxplot(perdel);
            ylabel('%\Delta: LER - BRDF albedo');
            set(gca,'xtick',[],'fontsize',16, 'ygrid','on');
            
            del = modis_lers(:) - modis_albs(:);
            f3 = figure; boxplot(del);
            ylabel('LER - BRDF albedo');
            set(gca,'xtick',[],'fontsize',16, 'ygrid','on');
            
            if nargout > 0
                varargout = {f1, f2, f3, f4};
            end
        end
        
        function test_black_sky_integration(DEBUG_LEVEL)
            % I'm curious how the BRDF can be greater than the black sky
            % albedo, it seems like since the black sky albedo is the BRDF
            % integrated over all viewing angles, then the albedo should be
            % greater. On the other hand, in the definition in Lucht et al.
            % 2000 (IEEE Trans. Geosci. Rem. Sense., p 977), the albedo is
            % the integral over viewing directions divided by pi, so
            % perhaps that factor is what makes the difference.

            if ~exist('DEBUG_LEVEL', 'var')
                DEBUG_LEVEL = 1;
            end
            
            % Get the site output files
            aerosol_subdir = 'aerosols_continental_mcd43d_more_szas';
            locs_info = ncinfo(fullfile(misc_alb_plots.scia_modis_dir, 'loc_coeffs_D_new.nc'));
            loc_names = cellstr(ncread(locs_info.Filename, 'ShortName'));
            loc_coeffs = ncread(locs_info.Filename, 'ModisCoefficients');
            loc_coeffs = squeeze(loc_coeffs(:,3,:,:)); % always use band 3
            loc_months = double(ncread(locs_info.Filename, 'Months'));
            
            i_month = 12; % month, later should vary this
            i_site = 2; % start with atlanta, vary this later
            
            this_coeffs = loc_coeffs(:,i_month,i_site);
            loc_file = ncinfo(fullfile(misc_alb_plots.scia_modis_dir, aerosol_subdir, sprintf('%s-LER.nc', loc_names{i_site})));
            
            sza_vec = ncread(loc_file.Filename, 'sza');
            vza_vec = ncread(loc_file.Filename, 'vza');
            raa_vec = 180 - ncread(loc_file.Filename, 'raa'); % SCIATRAN RAA definition is backwards from MODIS
            black_sky = nan(size(sza_vec));
            black_sky_2 = nan(size(sza_vec));
            for a=1:numel(sza_vec)
                [vzas, raas] = meshgrid(vza_vec, raa_vec);
                modis_albs = modis_brdf_alb(this_coeffs(1), this_coeffs(2), this_coeffs(3), sza_vec(a), vzas, raas);
                black_sky(a) = integrate_black_sky(this_coeffs, sza_vec(a));
                black_sky_2(a) = integrate_black_sky_roman(this_coeffs, sza_vec(a));
                
                figure; 
                contour(vzas, raas, modis_albs,'fill','on');
                hold on
                [c_matrix, c_handle] = contour(vzas, raas, modis_albs, [black_sky(a), black_sky(a)], 'color','k', 'linewidth', 2, 'linestyle', '--');
                clabel(c_matrix, c_handle);
                xlabel('Viewing zenith angle');
                ylabel('Relative azimuth angle');
                cb = colorbar;
                cb.Label.String = 'BRDF reflectance';
                title(sprintf('SZA = %g', sza_vec(a)));
            end
            
            figure;
            plot(sza_vec, black_sky, 'ko');
            hold on
            plot(sza_vec, black_sky_2, 'bx');
            legend('Black sky alb', 'BSA, alt calc.');
            xlabel('Solar zenith angle');
            ylabel('Albedo');
        end
        
        function compare_c1integral_c3(DEBUG_LEVEL)
            % Compute the integrate of the MCD43C1 BRDF product against the
            % MCD43C3 black sky albedo product. In theory, these should be
            % similar/the same, since according to
            % https://www.umb.edu/spectralmass/terra_aqua_modis/v006/mcd43c3,
            % the C3 product is calcuated from the C1 product. The reason I
            % care is that I found that definitions of black-sky albedo may
            % or may not include a factor of 1/pi in the integral over all
            % viewing directions. Since that factor is what makes the BSA <
            % BRDF in summer, it matters whether it is in the MODIS BSA
            % product.
            
            
            
            test_date = ask_date('Enter the date to test (2012-06-01 to 2012-06-08)');
            % Load the MODIS C1 BRDF coefficients and black sky albedo for
            % one day over the US
            year_str = datestr(test_date, 'yyyy');
            modis_day = modis_date_to_day(test_date);
            c1_file = dirff(fullfile(misc_alb_plots.mcd43c1_dir, year_str, sprintf('MCD43C1.A%s%03d*.hdf', year_str, modis_day)));
            c1_info = hdfinfo(c1_file(1).name);
            
            c3_file = dirff(fullfile(misc_alb_plots.mcd43c3_dir, year_str, sprintf('MCD43C3.A%s%03d*.hdf', year_str, modis_day)));
            c3_info = hdfinfo(c3_file(1).name);
            
            
            [modis_lon, modis_lat, modis_xx, modis_yy] = modis_cmg_latlon(0.05, [-125 -65], [25 50], 'grid');
            
            for a=1:3
                c1_coeffs(:,:,a) = hdfreadmodis(c1_info.Filename, hdfdsetname(c1_info, 1, 1, sprintf('BRDF_Albedo_Parameter%d_Band3',a)));
            end
            c1_coeffs = c1_coeffs(modis_yy, modis_xx, :);
            c1_coeffs = permute(c1_coeffs, [3 1 2]);
            
            c3_bsa = hdfreadmodis(c3_info.Filename, hdfdsetname(c3_info, 1, 1, 'Albedo_BSA_Band3'));
            c3_bsa = c3_bsa(modis_yy, modis_xx);
            
            % Calculate the integrals from the C1 BRDF
            c1_integrated_bsa = nan(size(c3_bsa));
            for a=1:numel(c3_bsa)
                if mod(a,1000) == 1
                    fprintf('Integral %d of %d\n', a, numel(c3_bsa));
                end
                
                if all(~isnan(c1_coeffs(:,a)))
                    % p. 29 of the MOD43 theoretical basis document
                    % indicated that the MCD43C3 BSA is for SZA = 45 deg
                    sza = solar_zenith_from_time(test_date, modis_lon(a), modis_lat(a), 'noon');
                    c1_integrated_bsa(a) = integrate_black_sky(c1_coeffs(:,a), sza);
                end
            end
            
            figure;
            pcolor(modis_lon, modis_lat, c3_bsa);
            shading flat
            colorbar;
            state_outlines('k');
            title(sprintf('MCD43C3 BSA %s', datestr(test_date)));
            
            
            figure;
            pcolor(modis_lon, modis_lat, c1_integrated_bsa);
            shading flat
            colorbar;
            state_outlines('k');
            title(sprintf('MCD43C1 integrated %s', datestr(test_date)));
            
            bsa_del = c3_bsa - c1_integrated_bsa;
            figure;
            pcolor(modis_lon, modis_lat, bsa_del);
            shading flat
            colorbar;
            colormap(blue_red_cmap);
            caxis(calc_plot_limits(bsa_del(:), 'difference'));
            state_outlines('k');
            title(sprintf('MCD43C3 BSA - MCD43C1 integrated %s', datestr(test_date)));
        end
    end
    
end

function [avg_value, lon, lat] = manual_average(filepath, avg_field, raa_rule)
p = inputParser;
p.addParameter('avg_field', 'BEHRColumnAmountNO2Trop');
p.addParameter('raa_rule', @(raa) true(size(raa)) );

F = dir(fullfile(filepath,'OMI_BEHR_*.mat'));
do_init = true;
for a=1:numel(F)
    fprintf('Averaging %s\n', F(a).name);
    D = load(fullfile(filepath, F(a).name), 'OMI');
    if do_init
        sum_value = nan(size(D.OMI(1).Longitude));
        sum_weight = nan(size(D.OMI(1).Latitude));
        lon = D.OMI(1).Longitude;
        lat = D.OMI(1).Latitude;
        do_init = false;
    end
    for b=1:numel(D.OMI)
        omi = omi_pixel_reject(D.OMI(b),'omi',0.2,'XTrackFlags');
        
        rr = raa_rule(omi.RelativeAzimuthAngle);
        omi.Areaweight(~rr) = 0;
        sum_value = nansum2(cat(3, sum_value, omi.(avg_field) .* omi.Areaweight),3);
        sum_weight = nansum2(cat(3, sum_weight, omi.Areaweight),3);
    end
end

avg_value = sum_value ./ sum_weight;

end

function bs_alb = integrate_black_sky(modis_coeffs, sza)
v_vec = linspace(0,pi/2,90);
r_vec = linspace(0,2*pi,360);
[vzas, raas] = meshgrid(v_vec, r_vec);

K = cell(1,3);
K{1} = ones(size(vzas));
[K{2}, K{3}] = modis_brdf_kernels(sza, vzas*180/pi, raas*180/pi);

h_k = nan(size(K));
for a=1:numel(K)
    vza_integral = trapz(vzas(1,:), K{a}.*sin(vzas).*cos(vzas), 2);
    h_k(a) = trapz(raas(:,1), vza_integral) / pi;
end

bs_alb = sum(modis_coeffs .* h_k');
end

function bs_alb = integrate_black_sky_alt(modis_coeffs, sza)
% This was a check that whether I added up the three parameters first and
% then integrated or integrated each parameter separately I got the same
% answer. I did, at least for Atlanta in June.
v_vec = linspace(0,pi/2,90);
r_vec = linspace(0,2*pi,360);
[vzas, raas] = meshgrid(v_vec, r_vec);

K = modis_brdf_alb(modis_coeffs(1), modis_coeffs(2), modis_coeffs(3), repmat(sza,size(vzas)), vzas*180/pi, raas*180/pi);
vza_integral = trapz(vzas(1,:), K .* sin(vzas) .* cos(vzas), 2);
bs_alb = 1/pi .* trapz(raas(:,1), vza_integral);
end

function bs_alb = integrate_black_sky_roman(modis_coeffs, sza)
% Roman et al. (Rem. Sense. Environ. 114 p. 738) also defines black sky
% albedo in Eq. 54a. However, there's some differences at first glance:
%   1) Their R is the BRF, not the BRDF. But, they define the BRDF in Eq.
%   41 the same way as the MOD43 TBD does the BRDF, so I think there's a
%   terminology issue - I think the MODIS BRDF may actually be a BRF.
%
%   2) They appear to be missing a factor a sin(vza), but that is taken up
%   in the change of variables from integrating over vza to over cos(vza).
v_vec = linspace(0,pi/2,90);
r_vec = linspace(0,2*pi,360);
[vzas, raas] = meshgrid(v_vec, r_vec);

K = modis_brdf_alb(modis_coeffs(1), modis_coeffs(2), modis_coeffs(3), repmat(sza,size(vzas)), vzas*180/pi, raas*180/pi);
% The change of variable from vza ( = \theta_v) to cos(vza) = \mu changes
% the differential:
%   d\mu = (\partial \mu / \partial \theta_v) * d\theta_v
%        = (\partial cos(\theta_v) / \partial \theta_v) * d\theta_v
%        = -sin(theta_v) * d\theta_v
% The integration limits are also flipped (they integrate from cos(pi/2) =
% 0 to cos(0) = 1 instead of 0 to pi/2) which is where the negative sign
% goes.
vza_integral = trapz(cos(fliplr(vzas(1,:))), fliplr(K .* cos(vzas)), 2);
bs_alb = 1/pi .* trapz(raas(:,1), vza_integral);
end
