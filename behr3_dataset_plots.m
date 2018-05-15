function [  ] = behr3_dataset_plots(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

HOMEDIR = getenv('HOME');
root_dir = fullfile(HOMEDIR, 'Dropbox', 'Berkeley', 'Research', 'My Papers', 'BEHR3');
save_root = fullfile(root_dir, 'Images');
main_tex_file = fullfile(root_dir, 'BEHR-ESSD.tex');

albedo_jja = 'MODISAlbedo-JJA-Monthly.mat';
albedo_djf = 'MODISAlbedo-DJF-Monthly.mat';

no2vcd_jja = 'BEHRColumnAmountNO2Trop-JJA-Monthly.mat';
no2vcd_d_jja = 'BEHRColumnAmountNO2Trop-JJA-Daily.mat';
no2vcd_djf = 'BEHRColumnAmountNO2Trop-DJF-Monthly.mat';
no2vcd_d_djf = 'BEHRColumnAmountNO2Trop-DJF-Daily.mat';

no2visvcd_jja = 'BEHRColumnAmountNO2TropVisOnly-JJA-Monthly.mat';
no2visvcd_djf = 'BEHRColumnAmountNO2TropVisOnly-DJF-Monthly.mat';
no2visvcd_d_jja = 'BEHRColumnAmountNO2TropVisOnly-JJA-Daily.mat';
no2visvcd_d_djf = 'BEHRColumnAmountNO2TropVisOnly-DJF-Daily.mat';

do_save = true;
do_close = true;

make_long_plots = ask_yn('Some figures take a long time to produce - remake them now?');
if make_long_plots
    make_very_long_plots = ask_yn('Even the ones that take >12 hr?');
else
    make_very_long_plots = false;
end


insert_change_stats_table('average', main_tex_file, 'INCR-TABLE');

temperature_error_effect();
incr_no2_changes(no2vcd_jja, no2visvcd_jja, 'MonthlyIncrDiffs', 'map');
incr_no2_changes(no2visvcd_jja, no2visvcd_jja, fullfile('Supplement','MonthlyVisIncrDiffs'), 'map');
incr_no2_changes(no2vcd_d_jja, no2visvcd_d_jja, 'DailyIncrDiffs', 'map', 4:9);

incr_no2_changes(no2vcd_djf, no2visvcd_djf, fullfile('Supplement', 'DJF-MonthlyIncrDiffs'), 'map');
incr_no2_changes(no2vcd_d_djf, no2visvcd_d_djf, fullfile('Supplement', 'DJF-DailyIncrDiffs'), 'map', 4:9);
incr_no2_changes(no2vcd_jja, no2visvcd_jja, fullfile('Supplement', 'MonthlyIncrDiffsHists'), 'hist', [], {'remove_outliers', true});
incr_no2_changes(no2vcd_djf, no2visvcd_djf, fullfile('Supplement', 'DJF-MonthlyIncrDiffsHists'), 'hist', [], {'remove_outliers', true});
incr_no2_changes(no2vcd_d_jja, no2visvcd_d_jja, fullfile('Supplement', 'DailyIncrDiffsHists'), 'hist', 4:9, {'remove_outliers', true});
incr_no2_changes(no2vcd_d_djf, no2visvcd_d_djf, fullfile('Supplement', 'DJF-DailyIncrDiffsHists'), 'hist', 4:9, {'remove_outliers', true});

overall_trop_changes('OverallDiffs', 'map');
overall_vis_changes(fullfile('Supplement', 'OverallVisDiffs'), 'map');
alb_changes();
alb_supplement_changes(make_long_plots);
monthly_v_daily_final();
land_cover();


% The following figures take a long time so ask if we should spend the time
% doing them
if make_long_plots
    amf_reproduction_skill();
    daily_vs_monthly_ut();
    vis_cause_plots();
else
    warning('Some figures (that take a long time) not reproduced');
end

if make_very_long_plots
    insert_change_stats_table('individual', main_tex_file, 'INCR-INDIV-TABLE');
end

    function incr_no2_changes(no2_variable, vis_no2_variable, save_name, plot_type, increment_indices, extra_args)
        if ~exist('increment_indices', 'var')
            increment_indices = [];
        end
        
        if ~exist('plot_type', 'var')
            plot_type = 'map';
        end
        
        if ~exist('extra_args', 'var')
            extra_args = {};
        end
        
        figs_vector = gobjects(0);
        
        if do_plot_increment(1, increment_indices)
            figs_vector(1) = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v2_1C_dir, 'new_dir', misc_behr_update_plots.behr_nasa_only_dir,...
                'base_comparison_file', no2_variable, 'new_comparison_file', no2_variable, 'diff_type', 'rel', 'plot_type', plot_type, 'only_ocean', false, extra_args{:});
            title('');
        end
        
        if do_plot_increment(2, increment_indices)
            figs_vector(2) = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_nasa_only_dir, 'new_dir', misc_behr_update_plots.behr_nasa_brdfD_dir,...
                'base_comparison_file', no2_variable, 'new_comparison_file', no2_variable, 'diff_type', 'rel', 'plot_type', plot_type, 'only_ocean', false, extra_args{:});
            title('');
            if strcmpi(plot_type, 'map')
                caxis([-20 20]);
            end
        end
        
        if do_plot_increment(3, increment_indices)
            figs_vector(3) = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_nasa_brdfD_dir, 'new_dir', misc_behr_update_plots.behr_nasa_brdf_vis_dir,...
                'base_comparison_file', vis_no2_variable, 'new_comparison_file', vis_no2_variable, 'diff_type', 'rel', 'plot_type', plot_type, 'only_ocean', false, extra_args{:});
            title('');
        end
        
        if do_plot_increment(4, increment_indices)
            % This is a special case that if we want to plot the difference
            % the new daily profiles made, we need the base variable to be
            % the 'monthly' profiles average because there were no daily
            % profiles before this
            if regcmp(no2_variable, 'Daily')
                base_variable = strrep(no2_variable, 'Daily', 'Monthly');
            else
                base_variable = no2_variable;
            end
            
            figs_vector(4) = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_nasa_brdf_vis_dir, 'new_dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir,...
                'base_comparison_file', base_variable, 'new_comparison_file', no2_variable, 'diff_type', 'rel', 'plot_type', plot_type, 'only_ocean', false, extra_args{:});
            title('');
        end
        
        if do_plot_increment(5, increment_indices)
            figs_vector(5) = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_tempfix_dir, 'new_dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_wrftemp_dir,...
                'base_comparison_file', no2_variable, 'new_comparison_file', no2_variable, 'diff_type', 'rel', 'plot_type', plot_type, 'only_ocean', false, extra_args{:});
            title('');
            if strcmpi(plot_type, 'map')
                caxis([-10 10])
            end
        end
        
        if do_plot_increment(6, increment_indices)
            figs_vector(6) = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_wrftemp_dir, 'new_dir', misc_behr_update_plots.behr_final_dir,...
                'base_comparison_file', no2_variable, 'new_comparison_file', no2_variable, 'diff_type', 'rel', 'plot_type', plot_type, 'only_ocean', false, extra_args{:});
            title('');
        end
        
        if do_plot_increment(7, increment_indices)
            if regcmp(no2_variable, 'Monthly')
                % This arg will not affect the map plots, it is only used
                % by histograms
                only_ocean = true;
            else
                only_ocean = false;
            end
            
            figs_vector(7) = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_final_dir, 'new_dir', misc_behr_update_plots.behr_v3B_daily_fix_dir,...
                'base_comparison_file', no2_variable, 'new_comparison_file', no2_variable, 'diff_type', 'rel', 'plot_type', plot_type, 'only_ocean', only_ocean, extra_args{:});
            title('');
            if strcmpi(plot_type, 'map')
                % For the monthly profiles, only the ocean LUT changed
                % which is a small effect. For the daily profiles, the time
                % selection criteria also changed, which has a larger
                % effect.
                if regcmp(no2_variable, 'Monthly')
                    caxis([-1 1]);
                else
                    caxis([-10 10]);
                end
            end
        end
        
        if do_plot_increment(8, increment_indices)
            figs_vector(8) = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v3B_daily_fix_dir, 'new_dir', misc_behr_update_plots.behr_v3B_var_trop_dir,...
                'base_comparison_file', no2_variable, 'new_comparison_file', no2_variable, 'diff_type', 'rel', 'plot_type', plot_type, 'only_ocean', false, extra_args{:});
            title('');
            if strcmpi(plot_type, 'map')
                caxis([-20 20])
            end
        end
        
        if do_plot_increment(9, increment_indices)
            figs_vector(9) = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v3B_var_trop_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
                'base_comparison_file', no2_variable, 'new_comparison_file', no2_variable, 'diff_type', 'rel', 'plot_type', plot_type, 'only_ocean', false, extra_args{:});
            title('');
            if strcmpi(plot_type, 'map')
                caxis([-20 20])
            end
        end
        
        figs_vector(~ishandle(figs_vector)) = [];
        if numel(figs_vector) > 6
            dims = [0 3];
            do_center = false;
        else
            dims = [0 2];
            do_center = true;
        end
        fig_combo = combine_plots(figs_vector, 'dims', dims, 'scale', 1);
        %fig_combo.Position(3:4) = [2 2.5] .* fig_combo.Position(3:4);
        colormap(blue_red_cmap);
        close(figs_vector);
        
        if mod(numel(figs_vector),2) ~= 0 && do_center
            % if we had an odd number of figures, then center the bottom
            % one. Assume that the children are in last-in-first-out order.
            ch = findall(fig_combo.Children, 'type', 'axes');
            center_axes(ch(1), ch(3), ch(2));
        end
        
        label_subfigs(fig_combo, 'xshift', 0.2)
        
        if do_save
            save_all_formats(fig_combo, fullfile(save_root, save_name));
        end
        if do_close
            close(fig_combo);
        end
    end

    function overall_trop_changes(save_name, plot_type)
        if ~exist('plot_type', 'var')
            plot_type = 'map';
        end
        
        fig_jja_monthly = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v2_1C_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
            'base_comparison_file', no2vcd_jja, 'new_comparison_file', no2vcd_jja, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        fig_djf_monthly = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v2_1C_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
            'base_comparison_file', no2vcd_djf, 'new_comparison_file', no2vcd_djf, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        fig_jja_daily = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v2_1C_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
            'base_comparison_file', no2vcd_jja, 'new_comparison_file', no2vcd_d_jja, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        fig_djf_daily = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v2_1C_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
            'base_comparison_file', no2vcd_djf, 'new_comparison_file', no2vcd_d_djf, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        
        figs_vector = [fig_jja_monthly, fig_djf_monthly, fig_jja_daily, fig_djf_daily];
        fig_combo = combine_plots(figs_vector, 'dims', [2 2]);
        fig_combo.Position(3:4) = 2 .* fig_combo.Position(3:4);
        colormap(blue_red_cmap);
        close(figs_vector);
        
        label_subfigs(fig_combo, 'xshift', 0.15);
        if do_save
            save_all_formats(fig_combo, fullfile(save_root, save_name));
        end
        if do_close
            close(fig_combo);
        end
    end

    function overall_vis_changes(save_name, plot_type)
        if ~exist('plot_type', 'var')
            plot_type = 'map';
        end
        
        fig_jja_monthly_vis = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v2_1C_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
            'base_comparison_file', no2visvcd_jja, 'new_comparison_file', no2visvcd_jja, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        fig_djf_monthly_vis = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v2_1C_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
            'base_comparison_file', no2visvcd_djf, 'new_comparison_file', no2visvcd_djf, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        fig_jja_daily_vis = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v2_1C_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
            'base_comparison_file', no2visvcd_jja, 'new_comparison_file', no2visvcd_d_jja, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        fig_djf_daily_vis = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v2_1C_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
            'base_comparison_file', no2visvcd_djf, 'new_comparison_file', no2visvcd_d_djf, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        
        figs_vector = [fig_jja_monthly_vis, fig_djf_monthly_vis, fig_jja_daily_vis, fig_djf_daily_vis];
        fig_combo = combine_plots(figs_vector, 'dims', [2 2]);
        fig_combo.Position(3:4) = 2 .* fig_combo.Position(3:4);
        colormap(blue_red_cmap);
        close(figs_vector);
        
        label_subfigs(fig_combo, 'xshift', 0.15);
        if do_save
            save_all_formats(fig_combo, fullfile(save_root, save_name));
        end
        if do_close
            close(fig_combo);
        end
    end

    function alb_changes(plot_type)
        if ~exist('plot_type', 'var')
            plot_type = 'map';
        end
        
        fig_jja = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_nasa_only_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
            'base_comparison_file', albedo_jja, 'new_comparison_file', albedo_jja, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        fig_djf = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_nasa_only_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
            'base_comparison_file', albedo_djf, 'new_comparison_file', albedo_djf, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        
        fig_combo = combine_plots([fig_jja; fig_djf], 'dims', [1 2]);
        fig_combo.Position(3) = 2*fig_combo.Position(3);
        colormap(blue_red_cmap);
        close([fig_jja, fig_djf]);
        label_subfigs(fig_combo, 'xshift', 0.15);
        
        if do_save
            save_all_formats(fig_combo, fullfile(save_root, 'SurfReflDiffs'));
        end
        if do_close
            close(fig_combo);
        end
    end

    function alb_supplement_changes(regen_from_behr_files, plot_type)
        if ~exist('plot_type', 'var')
            plot_type = 'map';
        end
        
        fig_v5v6_jja = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_nasa_only_dir, 'new_dir', misc_behr_update_plots.behr_modisv6_dir,...
            'base_comparison_file', albedo_jja, 'new_comparison_file', albedo_jja, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        fig_v5v6_djf = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_nasa_only_dir, 'new_dir', misc_behr_update_plots.behr_modisv6_dir,...
            'base_comparison_file', albedo_djf, 'new_comparison_file', albedo_djf, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        fig_brdf_bsa_jja = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_modisv6_dir, 'new_dir', misc_behr_update_plots.behr_nasa_brdfD_dir,...
            'base_comparison_file', albedo_jja, 'new_comparison_file', albedo_jja, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        fig_brdf_bsa_djf = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_modisv6_dir, 'new_dir', misc_behr_update_plots.behr_nasa_brdfD_dir,...
            'base_comparison_file', albedo_djf, 'new_comparison_file', albedo_djf, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        [fig_brdf_bsa_boxplot, fig_tmp(1)] = misc_behr_update_plots.plot_bsa_to_brdf_indiv_pixel_changes(regen_from_behr_files);
        
        fig_vec = [fig_v5v6_jja, fig_v5v6_djf, fig_brdf_bsa_jja, fig_brdf_bsa_djf, fig_brdf_bsa_boxplot];
        close(fig_tmp);
        
        fig_combo = combine_plots(fig_vec, 'dims', [3 2]);
        close(fig_vec);
        fig_combo.Position(3:4) = [4/3, 2] .* fig_combo.Position(3:4);
        % The boxplot should be the first child (since Matlab's graphics
        % objects seem to follow a first-in/last-out order. Expand it so
        % that it fills both columns of the figure. Also move it up a bit
        % so the labels are inside the figure bounds. Then center it
        % between the colorbar of the right hand figure above it and the
        % left hand figure proper above it.
        fig_combo.Children(1).Position(2:3) = [1.2, 2] .* fig_combo.Children(1).Position(2:3);
        center_axes(fig_combo.Children(1), fig_combo.Children(5), fig_combo.Children(2));
        colormap(blue_red_cmap);
        label_subfigs(fig_combo, 'xshift', 0.25);
        % The bottom figure needs its text label shifted back to look nice
        fig_combo.Children(1).Children(1).Position(1) = -0.1;
        % Also its Y ticks aren't great by default
        fig_combo.Children(1).YTick = -80:40:80;
        
        if do_save
            save_all_formats(fig_combo, fullfile(save_root, 'Supplement', 'SuppSurfReflDiffs'));
        end
        if do_close
            close(fig_combo);
        end
    end

    function temperature_error_effect(plot_type)
        if ~exist('plot_type', 'var')
            plot_type = 'map';
        end
        
        fig_jja = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir, 'new_dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_tempfix_dir,...
            'base_comparison_file', no2vcd_jja, 'new_comparison_file', no2vcd_jja, 'diff_type', 'rel', 'plot_type', plot_type);
        caxis([-10 10]); title('');
        fig_djf = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_dir, 'new_dir', misc_behr_update_plots.behr_nasa_brdf_vis_profs_tempfix_dir,...
            'base_comparison_file', no2vcd_djf, 'new_comparison_file', no2vcd_djf, 'diff_type', 'rel', 'plot_type', plot_type);
        caxis([-10 10]); title('');
        
        fig_vector = [fig_jja, fig_djf];
        fig_combo = combine_plots(fig_vector, 'dims', [1 2], 'scale', 1);
        %fig_combo.Position(3) = 2*fig_combo.Position(3);
        colormap(blue_red_cmap);
        close(fig_vector);
        label_subfigs(fig_combo, 'xshift', 0.15);
        
        if do_save
            save_all_formats(fig_combo, fullfile(save_root, 'Supplement', 'TemperatureErrorIncrement'));
        end
        if do_close
            close(fig_combo);
        end
    end

    function monthly_v_daily_final(plot_type)
        if ~exist('plot_type', 'var')
            plot_type = 'map';
        end
        
        fig_jja = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v3B_surfpres_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
            'base_comparison_file', no2vcd_jja, 'new_comparison_file', no2vcd_d_jja, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        fig_djf = misc_behr_update_plots.plot_single_incr_diff('base_dir', misc_behr_update_plots.behr_v3B_surfpres_dir, 'new_dir', misc_behr_update_plots.behr_v3B_surfpres_dir,...
            'base_comparison_file', no2vcd_djf, 'new_comparison_file', no2vcd_d_djf, 'diff_type', 'rel', 'plot_type', plot_type);
        title('');
        
        all_plots = [fig_jja, fig_djf];
        fig_combo = combine_plots(all_plots, 'dims', [1 2]);
        fig_combo.Position(3) = 2*fig_combo.Position(3);
        colormap(blue_red_cmap);
        close(all_plots);
        label_subfigs(fig_combo, 'xshift', 0.15);
        
        if do_save
            save_all_formats(fig_combo, fullfile(save_root, 'DailyVsMonthlyVCDs'));
        end
        if do_close
            close(fig_combo);
        end
    end

    function daily_vs_monthly_ut()
        fprintf('SE US profiles...\n');
        [se_profs, se_prof_med, se_freq_dist, se_mean_bin, se_med_bin, fig_regions] = misc_behr_update_plots.plot_no2_vs_cloudfrac('2012-06-01','2012-08-31',false,'SE');
        drawnow nocallbacks
        
        fprintf('SE US shape factors...\n');
        [se_shape_facts, fig_tmp(1), fig_tmp(2), fig_tmp(3), fig_tmp(4), fig_tmp(5)] = misc_behr_update_plots.plot_no2_vs_cloudfrac('2012-06-01','2012-08-31',true,'SE');
        % Had some trouble with "AWT-EventQueue-0" exceptions, which some
        % Googling suggested to me it might have to do with too many
        % graphics draws trying to be deferred. So I'm going to try to
        % force the graphics draws after each plotting function and clean
        % up the unneeded figures right away.
        drawnow nocallbacks
        close(fig_tmp);
        clear fig_tmp
        
        fprintf('NW US profiles...\n');
        [nw_profs, nw_prof_med, nw_freq_dist, nw_mean_bin, nw_med_bin, fig_tmp(1)] = misc_behr_update_plots.plot_no2_vs_cloudfrac('2012-06-01', '2012-08-31',false,'NW');
        drawnow nocallbacks
        close(fig_tmp);
        clear fig_tmp
        
        fprintf('NW US shape factors...\n');
        [nw_shape_facts, fig_tmp(1), fig_tmp(2), fig_tmp(3), fig_tmp(4), fig_tmp(5)] = misc_behr_update_plots.plot_no2_vs_cloudfrac('2012-06-01', '2012-08-31',true,'NW');
        drawnow nocallbacks
        close(fig_tmp);
        
        prof_figs = [se_profs, se_prof_med, nw_profs, nw_prof_med];
        fig_combo = combine_plots(prof_figs, 'dims', [2 2], 'scale', 1);
        %fig_combo.Position(3:4) = [3 2] .* fig_combo.Position(3:4);
        label_subfigs(fig_combo, 'xshift', 0.15);
        close(prof_figs);
        
        freq_figs = [se_freq_dist, nw_freq_dist, se_shape_facts, nw_shape_facts];
        fig_freq_combo = combine_plots(freq_figs, 'dims', [2 2], 'scale', 1);
        %fig_freq_combo.Position(3:4) = 2*fig_freq_combo.Position(3:4);
        label_subfigs(fig_freq_combo, 'xshift', 0.15);
        close(freq_figs);
        
        bin_figs = [se_mean_bin, nw_mean_bin, se_med_bin, nw_med_bin];
        fig_bin_combo = combine_plots(bin_figs, 'dims', [2 2], 'scale', 1);
        %fig_bin_combo.Position(3:4) = 2*fig_bin_combo.Position(3:4);
        label_subfigs(fig_bin_combo, 'xshift', 0.15);
        close(bin_figs);
        
        if do_save
            save_all_formats(fig_combo, fullfile(save_root, 'ProfileMeansAndMedians'));
            save_all_formats(fig_freq_combo, fullfile(save_root, 'UTFreqDist'));
            save_all_formats(fig_bin_combo, fullfile(save_root, 'Supplement', 'UTBinnedByCF'));
            save_all_formats(fig_regions, fullfile(save_root, 'Supplement', 'Profile-Average-Regions'));
        end
        if do_close
            close([fig_combo, fig_freq_combo, fig_bin_combo, fig_regions]);
        end
    end

    function land_cover
        fig = misc_behr_update_plots.plot_ground_cover;
        if do_save
            save_all_formats(fig, fullfile(save_root, 'Supplement', 'MODIS-Land-Type'));
        end
        if do_close
            close(fig);
        end
    end

    function insert_change_stats_table(avg_or_indiv, tex_file, tex_marker)
        args = {'JJA', 'Monthly';...
            'DJF', 'Monthly';...
            'JJA', 'Daily';...
            'DJF', 'Daily'};
        
        is_daily = strcmpi('daily', args(:,2));
        first_monthly = find(~is_daily, 1);
        first_daily = find(is_daily, 1);
        
        
        for i=1:size(args,1)
            fprintf('  Tabulating %s %s\n', args{i,:});
            [all_diff_arrays{i}, this_rownames, this_colnames] = misc_behr_update_plots.tabulate_average_differences('variables', {'BEHRColumnAmountNO2Trop'},...
                'comparison_type', avg_or_indiv, 'months', args{i,1}, 'prof_type', lower(args{i,2}), 'use_vis_incr', true);
            rownames{i} = cat(1, repmat(args(i,2), 1, numel(this_rownames)), this_rownames);
            colnames{i} = cat(1, repmat(args(i,1), 1, size(this_colnames,2)), this_colnames(2:end,:));
        end
        
        n_cols = unique(cellfun(@(x) size(x,2), all_diff_arrays));
        if numel(n_cols) ~= 1
            error('Not all the output diff arrays had the same number of columns');
        end
        
        monthly_diff_array = cat(2, all_diff_arrays{~is_daily});
        monthly_rownames = rownames{first_monthly}';
        all_colnames = cat(2, colnames{~is_daily});
        
        daily_diff_array = cat(2, all_diff_arrays{is_daily});
        daily_rownames = rownames{first_daily}';
        
        diff_array = cat(1, monthly_diff_array, daily_diff_array);
        all_rownames = cat(1, monthly_rownames, daily_rownames);
        
        make_latex_table(diff_array, 'rownames', all_rownames, 'colnames', all_colnames, 'm2l', {'%.2g', 2, 'asym'},...
            'lines', {'\tophline', '\middlehline', '\bottomhline'}, 'extra_hlines', size(monthly_diff_array,1),...
            'file', tex_file, 'insert', true, 'marker', tex_marker,...
            'caption', ['Percent differences in \chem{NO_2} VCDs for each increment. Means are given with $1 \sigma$ uncertainties; ',...
                'medians are given with uncertainties as the distance to the upper and lower quartiles. Outliers were removed before calculating these statistics. ',...
                '*Statistics for visible-only \chem{NO_2} column. **Statistics only for ocean pixels.'],...
            'label', 'tab:incr-diffs');
    end

    function amf_reproduction_skill()
        fig = misc_behr_update_plots.plot_amf_reproduction_skill('diff_type','rel','plot_type','hist');
        label_subfigs(fig, 'xshift', 0.1);
        ch = get(fig, 'children');
        for i = 2:2:numel(ch)
            ch(i).Title.String = '';
        end
        
        if do_save
            save_all_formats(fig, fullfile(save_root, 'Supplement', 'AMF-Reproduction-Skill'));
        end
        if do_close
            close(fig);
        end
    end

    function vis_cause_plots()
        [fig_s3, fig_cf, tmpfigs(1), tmpfigs(2), tmpfigs(3)] = misc_behr_update_plots.plot_vis_scatter(false);
        close(tmpfigs);
        
        vis_figs = [fig_s3, fig_cf];
        fig_combo = combine_plots(vis_figs, 'dims', [1 2]);
        fig_combo.Position(3) = 2*fig_combo.Position(3);
        label_subfigs(fig_combo, 'xshift', 0.15);
        close(vis_figs);
        
        if do_save
            save_all_formats(fig_combo, fullfile(save_root, 'VisChangeCauses'));
        end
        if do_close
            close(fig_combo);
        end
    end
end

function b = do_plot_increment(increment, incr_indices)
b = isempty(incr_indices) || ismember(increment, incr_indices);
end

function save_all_formats(fig, save_name)
save_dir = fileparts(save_name);
if ~exist(save_dir, 'dir')
    fprintf('Making directory "%s"\n', save_dir);
    mkdir(save_dir)
end

set(fig,'paperpositionmode','auto');
savefig(fig, [save_name, '.fig'], 'compact');
print(fig, save_name, '-dpng','-r0');
print_eps(fig, save_name);
end
