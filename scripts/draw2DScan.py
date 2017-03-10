#! /bin/env python

# Making various plots given the input limits
# This script is kindly provided by Olivier Bondu
# 
# USAGE: ./draw2DScan.py limits_grid.json lambda yt -o out2Ddraw
#
#
# NOTE: You may need to run a setup script to get missing python packages.
# The script should be in the current directory. Run like this:
# source setup_for_2Dlimits.sh

import os, sys, argparse
import copy
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
import matplotlib.path as mpath
import matplotlib.contour as mcontour

import HiggsAnalysis.bbggLimits.ParametersGrid as pg

# To Do: incorparate style file in here
# import cp3_llbb.CommonTools.CMSStyle as CMSStyle

#os.system('source setup_for_2Dlimits.sh')

class EmptyHandle(object):
    """
    A dummy class representing a empty Handle in the legend
    """
    pass

class EmptyHandleHandler(object):
    """
    A Legend handler for EmptyHandle
    """

    def legend_artist(self, legend, orig_handle,
                       fontsize, handlebox):
        """
        Set handlebox width at 0 and return a dummy artist

        Parameters
        ----------
        legend : :class:`matplotlib.legend.Legend` instance
            The legend for which these legend artists are being created.
        orig_handle : :class:`matplotlib.artist.Artist` or similar
            The object for which these legend artists are being created.
        fontsize : float or int
            The fontsize in pixels. The artists being created should
            be scaled according to the given fontsize.
        handlebox : :class:`matplotlib.offsetbox.OffsetBox` instance
            The box which has been created to hold this legend entry's
            artists. Artists created in the `legend_artist` method must
            be added to this handlebox inside this method.
        """

        handlebox.set_width(0)

        artist = mlines.Line2D([], [], visible=False)

        handlebox.add_artist(artist)

        return artist

## Matplotlib hack
# Override algorithm in matplotlib.contour to check if the contour label is not overlaid by one of our point

def locate_label(self, linecontour, labelwidth):
    """
    Find a good place to plot a label (relatively flat
    part of the contour).
    """

    nsize = len(linecontour)

    # SB: Always cut the contour line in 50 segments
    # if labelwidth > 1:
        # xsize = int(np.ceil(nsize / labelwidth))
    # else:
        # xsize = 1

    xsize = 50

    if xsize == 1:
        ysize = nsize
    else:
        ysize = int(labelwidth)

    XX = np.resize(linecontour[:, 0], (xsize, ysize))
    YY = np.resize(linecontour[:, 1], (xsize, ysize))
    # I might have fouled up the following:
    yfirst = YY[:, 0].reshape(xsize, 1)
    ylast = YY[:, -1].reshape(xsize, 1)
    xfirst = XX[:, 0].reshape(xsize, 1)
    xlast = XX[:, -1].reshape(xsize, 1)
    s = (yfirst - YY) * (xlast - xfirst) - (xfirst - XX) * (ylast - yfirst)
    L = np.sqrt((xlast - xfirst) ** 2 + (ylast - yfirst) ** 2).ravel()
    dist = np.add.reduce(([(abs(s)[i] / L[i]) for i in range(xsize)]), -1)
    x, y, ind = self.get_label_coords(dist, XX, YY, ysize, labelwidth)

    # There must be a more efficient way...
    lc = [tuple(l) for l in linecontour]
    dind = lc.index((x, y))

    return x, y, dind

def too_close(self, x, y, lw):
    "Return *True* if a label is already near this location."

    global ax

    # SB: check if the label is close to a point
    all_points = [
            (data['observed_allowed_x'], data['observed_allowed_y']),
            (data['observed_excluded_x'], data['observed_excluded_y']),
            (data['expected_allowed_x'], data['expected_allowed_y']),
            (data['expected_excluded_x'], data['expected_excluded_y'])
            ]

    t = ax.transData

    for points in all_points:
        for loc in zip(*points):
            loc = t.transform(loc)
            d = np.sqrt((x - loc[0]) ** 2 + (y - loc[1]) ** 2)
            if d < 1 * lw:
                return True

    return False

# Change the two matplotlib fuctions to ours
mcontour.ContourLabeler.locate_label = locate_label
mcontour.ContourLabeler.too_close = too_close

def getHalfDiamondPath(invert=False):
    t = mtransforms.Affine2D().translate(-0.5, -0.5).rotate_deg(45).scale(1.1, 1.2)

    # Half square
    vertices = [
            [0.0, 0.0], [1.0, 0.0],
            [1.0, 0.5], [0.0, 0.5],
            [0.0, 0.0]
            ]

    codes = [
            mpath.Path.MOVETO, mpath.Path.LINETO,
            mpath.Path.LINETO, mpath.Path.LINETO,
            mpath.Path.CLOSEPOLY,
            ]

    p = mpath.Path(vertices, codes)

    if invert:
        t.rotate_deg(180)

    return t.transform_path(p)

def getDiamondPath():
    t = mtransforms.Affine2D().translate(-0.5, -0.5).rotate_deg(45).scale(1.1, 1.2)

    # Full square
    vertices = [
            [0.0, 0.0], [1.0, 0.0],
            [1.0, 1.0], [0.0, 1.0],
            [0.0, 0.0]
            ]

    codes = [
            mpath.Path.MOVETO, mpath.Path.LINETO,
            mpath.Path.LINETO, mpath.Path.LINETO,
            mpath.Path.CLOSEPOLY,
            ]

    p = mpath.Path(vertices, codes)

    return t.transform_path(p)


# Options

parser = argparse.ArgumentParser(description='Draw comparison of expected limits for all clusters & MVA')
parser.add_argument('input', action='store', type=str, help='JSON file containing the limits for all the points')
parser.add_argument('-o', '--output', action='store', type=str, help='Output directory', required=True)
parser.add_argument('-u', '--unblind', action='store_true', dest='unblind', help='If set, draw also observed upper limits')
parser.add_argument('-r', '--rescale-to-hh-br', action='store_true', dest='rescale_to_hh', help='If set, limits are rescaled to the HH BR')
parser.add_argument('--no-latex', action='store_true', dest='no_latex', help='Do not create LaTeX table of limits')

parser.add_argument('x', action='store', type=str, help='X parameter')
parser.add_argument('y', action='store', type=str, help='Y parameter')

parser.add_argument('--lambda', action='store', type=float, dest='kl', default=1, help='default value for k_l')
parser.add_argument('--yt', action='store', type=float, dest='yt', default=1, help='default value for k_t')
parser.add_argument('--c2', action='store', type=float, dest='c2', default=0, help='default value for c2')
parser.add_argument('--cg', action='store', type=float, dest='cg', default=0, help='default value for cg')
parser.add_argument('--c2g', action='store', type=float, dest='c2g', default=0, help='default value for c2g')

parser.add_argument('--pdf', action='store_true', dest='pdf', default=False, help='Make PDF images as well')

options = parser.parse_args()

limits = None
with open(options.input) as f:
    limits = json.load(f)

parameters = ['lambda', 'yt', 'c2', 'cg', 'c2g']

default_values = {
        "lambda": 1,
        "yt": 1,
        "c2": 0,
        "cg": 0,
        "c2g": 0
    }

parameter_values = {
        'lambda': options.kl,
        'yt': options.yt,
        'c2': options.c2,
        'cg': options.cg,
        'c2g': options.c2g
        }

parameter_position = {
        'lambda': 0,
        'yt': 1,
        'c2': 2,
        'cg': 3,
        'c2g': 4
        }

parameter_legend = {
        'lambda': r'\kappa_\lambda',
        'yt': r'\kappa_t',
        'c2': r'c_2',
        'cg': r'c_g',
        'c2g': r'c_{2g}'
        }

if not options.x in parameters:
    raise Exception("Invalid X parameter")

if not options.y in parameters:
    raise Exception("Invalid Y parameter")

fixed_parameters = parameters[:]
fixed_parameters.remove(options.x)
fixed_parameters.remove(options.y)

def filterPoints(x):
    global fixed_parameters, parameter_values

    for p in fixed_parameters:
        if x[p] != parameter_values[p]:
            return False

    return True

points = pg.getPoints(filterPoints)

print points

br = 0.26 / 100.

data = {}
has_SM = False
sm_observed_excluded = False
sm_expected_excluded = False

for point in points:

    # print("Working on point %d" % point)

    is_SM = False
    if point == 324:
        has_SM = True
        is_SM = True

    all_x = data.setdefault('all_x', [])
    all_y = data.setdefault('all_y', [])

    observed_allowed_x = data.setdefault('observed_allowed_x', [])
    observed_allowed_y = data.setdefault('observed_allowed_y', [])
    expected_allowed_x = data.setdefault('expected_allowed_x', [])
    expected_allowed_y = data.setdefault('expected_allowed_y', [])

    observed_excluded_x = data.setdefault('observed_excluded_x', [])
    observed_excluded_y = data.setdefault('observed_excluded_y', [])
    expected_excluded_x = data.setdefault('expected_excluded_x', [])
    expected_excluded_y = data.setdefault('expected_excluded_y', [])

    expected = data.setdefault('expected', [])
    # Index 0 is DOWN error, index 1 is UP error
    one_sigma = data.setdefault('one_sigma', [[], []])
    two_sigma = data.setdefault('two_sigma', [[], []])
    observed = data.setdefault('observed', [])

    if not str(point) in limits:
        continue

    params = pg.getParametersFromPoint(point)

    point = str(point)

    expected.append(limits[point]['expected'])
    observed.append(limits[point].get('observed', 0))

    # UP error
    one_sigma[1].append(limits[point]['one_sigma'][1])
    # DOWN error
    one_sigma[0].append(limits[point]['one_sigma'][0])

    # UP error
    two_sigma[1].append(limits[point]['two_sigma'][1])
    two_sigma[0].append(limits[point]['two_sigma'][0])

    x_value = params[parameter_position[options.x]]
    y_value = params[parameter_position[options.y]]

    xs, down, up = pg.getCrossSectionForParameters(*params)
    if is_SM:
        if observed[-1] < (xs * br) :
            sm_observed_excluded = True

        if expected[-1] < (xs * br) :
            sm_expected_excluded = True

    else:
        if observed[-1] < (xs * br) :
            observed_excluded_x.append(x_value)
            observed_excluded_y.append(y_value)
        else:
            observed_allowed_x.append(x_value)
            observed_allowed_y.append(y_value)

        if expected[-1] < (xs * br) :
            expected_excluded_x.append(x_value)
            expected_excluded_y.append(y_value)
        else:
            expected_allowed_x.append(x_value)
            expected_allowed_y.append(y_value)

    all_x.append(x_value)
    all_y.append(y_value)

data['expected'] = np.asarray(data['expected'])
data['observed'] = np.asarray(data['observed'])
data['one_sigma'] = np.asarray(data['one_sigma'])
data['two_sigma'] = np.asarray(data['two_sigma'])

data['observed_allowed_x'] = np.asarray(data['observed_allowed_x'])
data['observed_allowed_y'] = np.asarray(data['observed_allowed_y'])
data['expected_allowed_x'] = np.asarray(data['expected_allowed_x'])
data['expected_allowed_y'] = np.asarray(data['expected_allowed_y'])
data['observed_excluded_x'] = np.asarray(data['observed_excluded_x'])
data['observed_excluded_y'] = np.asarray(data['observed_excluded_y'])
data['expected_excluded_x'] = np.asarray(data['expected_excluded_x'])
data['expected_excluded_y'] = np.asarray(data['expected_excluded_y'])

if options.rescale_to_hh:
    data['expected'] /= br
    data['observed'] /= br
    data['one_sigma'] /= br
    data['two_sigma'] /= br

# CMSStyle.changeFont()

# Create a figure instance
fig = plt.figure(1, figsize=(7, 7), dpi=80)

# Create an axes instance
ax = fig.add_subplot(111)

ax.set_ylabel('$%s$' % parameter_legend[options.y], fontsize='x-large', y=0.85)
ax.set_xlabel('$%s$' % parameter_legend[options.x], fontsize='x-large', x=0.85, labelpad=-2)
ax.margins(0.1, 0.1)

fig.tight_layout()

# ax.set_autoscale_on(False)

# Create isoline of cross-sections

xmin, xmax = min(data['all_x']), max(data['all_x'])
# padding = (xmax - xmin) * 0.1
# xmin -= padding
# xmax += padding

ymin, ymax = min(data['all_y']), max(data['all_y'])
# padding = (ymax - ymin) * 0.1
# ymin -= padding
# ymax += padding

xaxis = np.linspace(xmin, xmax, 100)
yaxis = np.linspace(ymin, ymax, 100)
x, y = np.meshgrid(xaxis, yaxis)

params = [parameter_values['lambda'], parameter_values['yt'], parameter_values['c2'], parameter_values['cg'], parameter_values['c2g']]

params[parameter_position[options.x]] = x
params[parameter_position[options.y]] = y

z = pg.getCrossSectionForParameters(*params)

isolines = np.array([1, 5, 10, 25, 50, 75, 100, 150, 200, 500])

if options.x == "lambda" and options.y == "yt" and options.c2 == -3:
    isolines = np.array([10, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400])

if options.rescale_to_hh:
    isolines *= 100

iso = ax.contour(x, y, z[0] * br if not options.rescale_to_hh else z[0], levels=isolines, colors='gray', linestyles='dashdot')

# ax.grid()

angle = 45
theta1 = angle
theta2 = angle + 180

expected_markers = mpath.Path.wedge(theta1, theta2)
observed_markers = mpath.Path.wedge(theta2, theta1)

expected_sm_marker = getHalfDiamondPath(True)
observed_sm_marker = getHalfDiamondPath()

excluded_color = '#8E2800'
allowed_color = '#468966'

markersize = 10

if mpl.__version__ == '1.2.1':
    # Older version of matplotlib have a bug when using custom Path as marker. The size of the symbol is twice larger
    markersize /= 2

# Allowed
ax.plot(data['observed_allowed_x'], data['observed_allowed_y'], ls='', lw=0, marker=observed_markers, 
        ms=markersize, mew=1, mec=allowed_color, mfc='none', zorder=100, label='Allowed')
if has_SM and not sm_observed_excluded:
    ax.plot(default_values[options.x], default_values[options.y], ls='', lw=0, marker=observed_sm_marker, 
            ms=markersize*1.3, mew=2, mec=allowed_color, mfc='none', zorder=100, label='Allowed')
ax.plot(data['expected_allowed_x'], data['expected_allowed_y'], ls='', lw=0, marker=expected_markers,
        ms=markersize, mew=1, mec=allowed_color, mfc='none', zorder=100, label='Allowed')
if has_SM and not sm_expected_excluded:
    ax.plot(default_values[options.x], default_values[options.y], ls='', lw=0, marker=expected_sm_marker,
            ms=markersize*1.15, mew=2, mec=allowed_color, mfc='none', zorder=100, label='Allowed')

# Excluded
ax.plot(data['observed_excluded_x'], data['observed_excluded_y'], ls='', ms=markersize, lw=0, 
        marker=observed_markers, color=excluded_color, mec=excluded_color, mew=1, zorder=100, label='Excluded (95% CL)')
if has_SM and sm_observed_excluded:
    ax.plot(default_values[options.x], default_values[options.y], ls='', ms=markersize*1.15, lw=0, 
            marker=observed_sm_marker, color=excluded_color, mew=2, mec=excluded_color, zorder=100)
ax.plot(data['expected_excluded_x'], data['expected_excluded_y'], ls='', ms=markersize, lw=0, 
        marker=expected_markers, color=excluded_color, mec=excluded_color, mew=1, zorder=100, label='Excluded (95% CL)')
if has_SM and sm_expected_excluded:
    ax.plot(default_values[options.x], default_values[options.y], ls='', ms=markersize*1.15, lw=0, 
            marker=expected_sm_marker, color=excluded_color, mew=2, mec=excluded_color, zorder=100)

# SM

ax.clabel(iso, inline=1, fontsize='medium', fmt='%.0f fb')

ax.minorticks_on()

parameters_formatted_text = []
for p in fixed_parameters:
    parameter_value = parameter_values[p]
    if parameter_value == default_values[p]:
        parameter_value = parameter_legend[p] + '^{SM}'
    parameters_formatted_text.append("${} = {}$".format(parameter_legend[p], parameter_value))

theory_text = ', '.join(parameters_formatted_text)

#ax.text(0.080, 0.955, r"$pp \rightarrow hh \rightarrow b\bar{b}VV \rightarrow b\bar{b}l\nu l\bar{\nu}$", transform=ax.transAxes, ha='left', va='baseline')
ax.text(0.00, 1.03, r"CMS", transform=ax.transAxes, ha='left', va='baseline', size='xx-large')
ax.text(0.15, 1.03, r"$pp \rightarrow HH \rightarrow b\bar{b}\gamma\gamma$    "+theory_text, 
        transform=ax.transAxes, ha='left', va='center')
#ax.text(0.70, 1.03, theory_text, transform=ax.transAxes, ha='right', va='baseline')
ax.text(1.00, 1.03, '$36.5\,fb^{-1}$ (13 TeV)', transform=ax.transAxes, ha='right', va='baseline', size='large')

# Legend

fig.tight_layout()

fig.subplots_adjust(bottom=0.14, top=0.93)

allowed_expected_legend_marker = mlines.Line2D([], [], mec=allowed_color, marker=expected_markers, lw=0, mew=1, ms=markersize, mfc='none')
excluded_expected_legend_marker = mlines.Line2D([], [], color=excluded_color, mec=excluded_color, marker=expected_markers, lw=0, mew=1, ms=markersize)

allowed_observed_legend_marker = mlines.Line2D([], [], mec=allowed_color, marker=observed_markers, lw=0, mew=1, ms=markersize, mfc='none')
excluded_observed_legend_marker = mlines.Line2D([], [], mfc=excluded_color, mec=excluded_color, marker=observed_markers, lw=0, mew=1, ms=markersize)

handles = [EmptyHandle(), EmptyHandle()]
labels  = ["Expected:", "Observed:"]

handles += [allowed_expected_legend_marker, allowed_observed_legend_marker, excluded_expected_legend_marker, excluded_observed_legend_marker]
labels += ["allowed", "allowed", "excluded (95% CL)", "excluded (95% CL)"]

if has_SM:
    sm_legend_marker = mlines.Line2D([], [], mec='black', marker=getDiamondPath(), lw=0, mew=2, ms=markersize, mfc='none')
    handles += [sm_legend_marker]
    labels += ["SM"]
# else:
    # handles += [EmptyHandle()]
    # labels += [""]

# handles += [EmptyHandle()]
# labels += ["(95% CL)"]

handles += [mlines.Line2D([], [], ls='dashdot', marker=None, lw=1, color='gray')]
labels += ["Theory"]

text_x_pos = 0.045

lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(-0.04, -0.06),
        frameon=False, numpoints=1, fontsize='medium', ncol=4, columnspacing=2,
        handler_map={EmptyHandle: EmptyHandleHandler()})

# ax.text(text_x_pos, -0.072, theory_text, transform=ax.transAxes, ha='left')

# CMSStyle.applyStyle(fig, ax, 2300)

plot_name = os.path.dirname(options.input).replace('Cards/', '').replace('/', '_') + '_2d_%s_vs_%s_scan' % (options.x, options.y)

names = []
for p in fixed_parameters:
    names.append('%s=%.1f' % (p, parameter_values[p]))

plot_name += '_'+'_'.join([x.replace('.', 'p') for x in names])

if options.rescale_to_hh:
    plot_name += '_rescaled_to_hh'


fig.savefig(os.path.join(options.output, plot_name + '.png'))
print("Plot saved as %r") % os.path.join(options.output, plot_name + '.png')
if options.pdf:
  fig.savefig(os.path.join(options.output, plot_name + '.pdf'))
  print("Plot saved as %r") % os.path.join(options.output, plot_name + '.pdf')

