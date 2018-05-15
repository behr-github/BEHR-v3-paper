function [ fighandle ] = plot_error_envelope_x( y_in, x_lower, x_upper, varargin )
%plot_error_envelope_x(y_in, x_lower, x_upper, [opt] colorspec, [opt] figure number): Uses the fill function to plot an error envelope for error in x value.
%   This plots an error envelope on your plot for error in the x_values.
%   This is set up to allow the use of standard deviation/std. error or
%   quartiles, hence the lower and upper bounds are input separately.  The
%   default color of the envelope is light grey, but the color can be
%   specified as an optional argument.  Most likely this will need to be
%   called before you plot your actual data, as the fill will probably
%   cover the line if plotted afterwards.
%
%   Returns the figure number it used to plot.

p = inputParser;
p.addRequired('y',@isnumeric);
p.addRequired('x_lower',@isnumeric);
p.addRequired('x_upper',@isnumeric);
p.addParamValue('colorspec',[0.8 0.8 0.8]);
p.addParamValue('facealpha',1,@isscalar);
p.addParamValue('fignum',-1);

p.parse(y_in, x_lower, x_upper, varargin{:});
pout = p.Results;

y = pout.y;
xl = pout.x_lower;
xu = pout.x_upper;

colspec = pout.colorspec;
nfig = pout.fignum;
facealpha = pout.facealpha;

% If no figure handle is passed, create one. If one is, switch to that
% figure and turn hold on so that we don't delete existing data.
if nfig == -1
    fighandle = figure();
else
    fighandle = figure(nfig);
    hold on
end

X = [xl fliplr(xu)];
Y = [y fliplr(y)];

fill(X,Y,colspec,'edgecolor','none','FaceAlpha',facealpha);
hold off

end

