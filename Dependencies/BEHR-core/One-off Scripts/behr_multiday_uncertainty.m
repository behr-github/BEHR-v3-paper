function [ ] = behr_multiday_uncertainty(start_date, end_date, prof_mode, region, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

start_date = validate_date(start_date);
end_date = validate_date(end_date);

dvec = start_date:end_date;
for d=1:numel(dvec)
    fprintf('Working on %s\n', datestr(dvec(d)));
    [Data,OMI] = load_behr_file(dvec(d), prof_mode, region);
    [Delta, DeltaGrid] = behr_uncertainty_estimation(Data,OMI,'MODISAlbedo',15,varargin{:});
    %{
    avg_change = BEHR_day_no2(DeltaGrid,'PercentChangeNO2');
    figure; 
    pcolor(DeltaGrid(2).Longitude, DeltaGrid(2).Latitude, avg_change);
    shading flat; 
    caxis([-10 10]);
    cb=colorbar; 
    cb.Label.String = '%\Delta VCD for +15% surf. refl.';
    state_outlines('k')
    title(datestr(dvec(d)));
    %}
    save(sprintf('BEHR_uncert_%s.mat', datestr(dvec(d),'yyyymmdd')), 'Delta', 'DeltaGrid');
end

end

