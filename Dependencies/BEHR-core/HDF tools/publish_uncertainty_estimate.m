function publish_uncertainty_estimate(change_field, varargin)
E = JLLErrors;

p = inputParser;
p.addOptional('output_location', '.');
p.addParameter('region', 'us');

p.parse(varargin{:});
pout = p.Results;

region = pout.region;
output_file = pout.output_location;
if isdir(output_file)
    output_file = fullfile(output_file, sprintf('BEHR-%s-uncertainty.hdf', region));
end

allowed_change_fields = {'PercentChangeNO2', 'PercentChangeNO2Vis', 'PercentChangeAMF', 'PercentChangeAMFVis'};
if ~ismember(change_field, allowed_change_fields)
    E.badinput('CHANGE_FIELD must be one of: %s', strjoin(allowed_change_fields, ', '));
end

% First we need to get what uncertainty parameters are
% available and for which months. The months should be the same
% for all parameters. We'll load the uncertainty at the same
% time

F = dir(behr_paths.BEHRUncertSubdir(region));
uncert_params = {F([F.isdir]).name};
% Need to remove ., .., and BaseCase
uncert_params(regcmp(uncert_params, '\.+')) = [];
uncert_params(strcmp(uncert_params, 'BaseCase')) = [];

% Assume for now that there's four month of uncertainty
uncertainties = make_empty_struct_from_cell(uncert_params);

% Load the lat/lon grid from a BEHR file - it's not included in the
% uncertainty files, which is a bit of a nuisance.
GridFile = load(fullfile(behr_paths.BEHRMatSubdir(region,'monthly'), behr_filename('2012-01-01','monthly',region)),'OMI');
grid_lon = GridFile.OMI(1).Longitude;
grid_lat = GridFile.OMI(1).Latitude;
grid_size = size(grid_lon);
lightning_lon = -95;
lightning_lat = 37.5;

for i_param = 1:numel(uncert_params)
    fprintf('Loading %s files: ', uncert_params{i_param});
    
    uncert_files_tmp = dirff(fullfile(behr_paths.BEHRUncertSubdir(region), uncert_params{i_param}, 'BEHR*.mat'));
    
    % Order the files DJF, MAM, JJA, SON. This way the order in the 3D
    % array makes sense
    if numel(uncert_files_tmp) ~= 4
        E.notimplemented('Not 4 uncertainty files')
    end
    regexes = {'DJF', 'MAM', 'JJA', 'SON'};
    uncert_names = {uncert_files_tmp.name};
    perm_vec = nan(size(uncert_files_tmp));
    for i_re = 1:numel(regexes)
        perm_vec(i_re) = find(regcmp(uncert_names, regexes{i_re}));
    end
    uncert_files = uncert_files_tmp(perm_vec);
    
    file_dates = cellfun(@(x) regexp(x, '(DJF|MAM|JJA|SON)', 'match', 'once'), {uncert_files.name}, 'UniformOutput', false);
    if i_param == 1
        check_file_dates = file_dates;
    elseif ~isequal(file_dates, check_file_dates)
        E.notimplemented('Different parameters produced for different months');
    end
    
    substruct = make_empty_struct_from_cell({'date', 'percent_diff'}, cell(size(uncert_files)));
    
    for i_file = 1:numel(uncert_files)
        UAvg = load(uncert_files(i_file).name);
        ErrorAvg = UAvg.ErrorAvg;
        
        fprintf('%d ', i_file);
        
        substruct(i_file).date = file_dates{i_file};
        if numel(ErrorAvg) == 1
            perdiff = ErrorAvg.(change_field);
        elseif numel(ErrorAvg) == 2
            % If we calculated a percent difference by raising
            % and lowering the perturbed value, this will
            % effectively calculate a mean percent difference.
            % I do it this way rather than summing and dividing
            % by two or summing the absolute value and dividing
            % by two to keep the sign. This assumes implicitly
            % that the change w.r.t. the input parameter is
            % monotonic; if, e.g. increasing or decreasing the
            % input parameter both increased the NO2 VCDs, then
            % this method will calculate a reduced uncertainty
            % than if we summed or just took one side. However,
            % I argue that is both unlikely and, even if it
            % happens, correct because the range over which the
            % VCDs varies is smaller in that case than if it
            % was monotonically increasing or decreasing with
            % the input parameter.
            perdiff = (ErrorAvg(2).(change_field) - ErrorAvg(1).(change_field))./2;
        else
            E.notimplemented('numel(ErrorAvg) > 2')
        end
        % Sometime the perturbation produces an absurdly high
        % uncertainty. Restrict to 3 sigma.
        %perdiff(abs(perdiff) > 3*nanstd(perdiff)) = nan;
        substruct(i_file).percent_diff = perdiff;
    end
    uncertainties.(uncert_params{i_param}) = substruct;
    fprintf('\n');
end


n_months = numel(check_file_dates);
n_params = numel(uncert_params);


if n_months ~= 4
    E.notimplemented('# of months ~= 4');
end

for i_month = 1:n_months
    pAvg = RunningAverage();
    pAvgAdjusted = RunningAverage();
    
    for i_param = 1:n_params
        this_struct = uncertainties.(uncert_params{i_param})(i_month);
        if ~strcmp(this_struct.date, check_file_dates{i_month})
            E.callError('wrong_date', 'Dates in substruct out of order compared to what was expected')
        end
        pAvg.addData((this_struct.percent_diff).^2);
        if strcmpi(uncert_params{i_param}, 'ProfileTime')
            % ProfileTime chooses a random daily profile from the same
            % month as the retrieval, which is probably overestimating the
            % error due to wind direction because we see in Laughner et al.
            % 2018 AMT that WRF gets the plumes right ~70% of the time. As
            % a better estimate then, we should reduce this error by ~70%.
            pAvgAdjusted.addData((this_struct.percent_diff .* 0.3).^2);
        else
            pAvgAdjusted.addData((this_struct.percent_diff).^2);
        end
    end
    perdiff = sqrt(pAvg.values);
    uncertainties.Total(i_month).date = this_struct.date; % assume its okay to use the last substructure's date
    uncertainties.Total(i_month).percent_diff = perdiff;
    uncertainties.AdjustedTotal(i_month).date = this_struct.date;
    
    adjusted_perdiff = sqrt(pAvgAdjusted.values);
    lightning_xx = grid_lon > lightning_lon & grid_lat < lightning_lat;
    adjusted_perdiff(lightning_xx) = max(adjusted_perdiff(lightning_xx), 100);
    uncertainties.AdjustedTotal(i_month).percent_diff = adjusted_perdiff;
end

h5create(output_file, '/Longitude', grid_size, 'Datatype', 'single', 'FillValue', single(behr_fill_val));
h5write(output_file, '/Longitude', single(grid_lon))
h5create(output_file, '/Latitude', grid_size, 'Datatype', 'single', 'FillValue', single(behr_fill_val));
h5write(output_file, '/Latitude', single(grid_lat))

uncert_params{end+1} = 'Total';
uncert_params{end+1} = 'AdjustedTotal';
for i_param = 1:numel(uncert_params)
    if ~isequal(size(uncertainties.(uncert_params{i_param})(1).percent_diff), grid_size)
        E.callError('inconsistent_grids', 'The grid size in the uncertainty files is different from the standard US domain grid');
    end
    
    if strcmpi(uncert_params{i_param}, 'total')
        dset_name = '/Total_Uncertainty';
        attributes.description = 'Total uncertainty (in percent) of BEHRColumnAmountNO2Trop: quadrature sum of the individual contributors';
    elseif strcmpi(uncert_params{i_param}, 'AdjustedTotal')
        dset_name = '/Adjusted_Total_Uncertainty';
        attributes.description = sprintf('Total uncertainty (in percent) of BEHRColumnAmountNO2Trop: adjusted for overestimated plume direction contribution in the calculation and error in SE US lightning. Uncertainty is the greater of 100%% or the quadrature sum (ProfileTime reduced by 70%%) in the area east of %.1f W and south of %.1f N', lightning_lon, lightning_lat);
    else
        dset_name = sprintf('/Uncertainty_from_%s', uncert_params{i_param}); 
        attributes.description = sprintf('Uncertainty (in percent) of BEHRColumnAmountNO2Trop due to %s', uncert_params{i_param});
    end
    
    perdiff_array = cat(3, uncertainties.(uncert_params{i_param}).percent_diff);
    perdiff_array(isnan(perdiff_array) | isinf(perdiff_array)) = behr_fill_val;
    perdiff_size = size(perdiff_array);
    h5create(output_file, dset_name, perdiff_size, 'Datatype', 'single', 'Fillvalue', single(behr_fill_val));
    h5write(output_file, dset_name, single(perdiff_array));
    
    attributes.dim_description = sprintf('1: (length %d) Latitude, 2: (length %d) Longitude, 3: (length %d) Time', perdiff_size(1), perdiff_size(2), perdiff_size(3));
    
    param_dates = {uncertainties.(uncert_params{i_param}).date};
    date_indices = num2cell(1:numel(param_dates));
    attributes.time_periods = strjoin(sprintfmulti('%d: %s', date_indices, param_dates), ', ');
    
    att_names = fieldnames(attributes);
    for i_att = 1:numel(att_names)
        h5writeatt(output_file, dset_name, att_names{i_att}, attributes.(att_names{i_att}));
    end
    
end

end