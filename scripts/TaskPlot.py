################################################################################
# This file is part of CMacIonize
# Copyright (C) 2018, 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# CMacIonize is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CMacIonize is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
################################################################################

##
# @file plot_tasks.py
#
# @brief Script to plot the task plot for a given file with task output.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# import modules
import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as pl

name = "ScalingTest_08_tasks.txt"

# change the default matplotlib settings to get nicer plots
pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (9.0, 5.0)
pl.rcParams["font.size"] = 14

# name labels and colours for the various task types
task_types = np.loadtxt(
    "ScalingTest_08_types.txt",
    delimiter="\t",
    dtype={"names": ("type", "label"), "formats": ("i4", "S100")},
    ndmin=1,
)
task_types.sort(0)
task_names = task_types["label"]
# GCC tends to prepend typeinfo strings with the number of characters in the
# string. We remove this bit from the task name.
for i in range(len(task_names)):
    if len(task_names[i]) > 10:
        task_names[i] = task_names[i][2:]
    else:
        task_names[i] = task_names[i][1:]
task_colors = {
    "NBasedResources": ["C0", 1.0],
    "ParticleGeometryResource": ["C0", 0.9],
    "GaussBasedResources": ["C0", 0.8],
    "InteractionTask": ["C0", 0.7],
    "TMatrixM0Task": ["C0", 0.6],
    "TMatrixMAllTask": ["C0", 0.5],
    "WignerDResources": ["C0", 0.4],
    "AlignmentAverageTask": ["C1", 1.0],
    "ResetTMatrixResourceTask": ["C2", 1.0],
    "ResetInteractionResourceTask": ["C2", 0.9],
    "AbsorptionCoefficientGrid": ["C5", 1.0],
    "AbsorptionCoefficientTask": ["C6", 1.0],
    "AbsorptionSpecialWignerDResources": ["C6", 0.9],
    "AbsorptionShapeAveragingTask": ["C7", 1.0],
}
for color in task_colors.values():
    c = matplotlib.colors.to_rgba(color[0])
    color[0] = (c[0], c[1], c[2], color[1])

# load the data
data = np.loadtxt(
    name,
    delimiter="\t",
    dtype={
        "names": ("thread", "start", "end", "type"),
        "formats": ("i4", "u8", "u8", "i4"),
    },
)

# get information about the system
nthread = data["thread"].max() + 1

## make the plot

fig, ax = pl.subplots(1, 1, sharex=True)

ttot_min = data["start"].min()
ttot_max = data["end"].max()
ttot_int = ttot_max - ttot_min

ax.set_xlim(
    ttot_min - 0.05 * (ttot_max - ttot_min),
    ttot_max + 0.05 * (ttot_max - ttot_min),
)
ax.axvline(x=ttot_min, linestyle="--", color="k", linewidth=0.8)
ax.axvline(x=ttot_max, linestyle="--", color="k", linewidth=0.8)

# loop over the threads
for i in range(nthread):
    # filter out the data for this thread
    thread = data[data["thread"] == i]

    # create the task plot
    bar = [((task["start"]), (task["end"] - task["start"])) for task in thread]
    colors = [
        task_colors[task_names[task["type"]].decode("utf-8")][0]
        for task in thread
    ]
    ax.broken_barh(bar, (i - 0.4, 0.8), facecolors=colors, edgecolor="none")

# add empty blocks for the legend
step_names = {
    "Single T-matrix": "C0",
    "Ensemble averaging": "C1",
    "Resource recyling": "C2",
    #  "Scattering coefficient": "C5",
    "Absorption coefficient": "C6",
    #  "Shape average": "C7",
}
for step_name in step_names:
    ax.plot([], [], color=step_names[step_name], label=step_name)

# add the legend and clean up the axes
ax.legend(loc="upper center", ncol=min(2, len(task_colors)))
ax.set_ylim(-1.0, nthread * 1.2)
yticks = np.arange(8)
yticklabels = ["Thread {0}".format(ytick + 1) for ytick in yticks]
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels)

xticks = np.linspace(ttot_min, ttot_max, 5)
xticklabels = np.linspace(0.0, 1.0, 5)
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)

# get the total idle time and add information to the title
# alltime /= nthread * nproc
# ax.set_title("Total empty fraction: {0:.2f} \%".format((1.0 - alltime) * 100.0))
ax.set_xlabel("Simulation time fraction")

ax.xaxis.set_minor_locator(AutoMinorLocator())

# finalize and save the plot
pl.tight_layout()
pl.savefig("taskplot.png", dpi=300, bbox_inches="tight")
pl.savefig("taskplot.pdf", dpi=300, bbox_inches="tight")
