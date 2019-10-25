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
import argparse

# parse the command line arguments
argparser = argparse.ArgumentParser(
    description="Plot task plot based on a given task output file."
)

argparser.add_argument("-n", "--name", action="store", required=True)
argparser.add_argument("-t", "--type-file", action="store", required=True)
argparser.add_argument("-o", "--output", action="store", required=True)
argparser.add_argument("-l", "--labels", action="store_true")

args = argparser.parse_args()

name = args.name

# change the default matplotlib settings to get nicer plots
pl.rcParams["text.usetex"] = True
pl.rcParams["figure.figsize"] = (12, 10)
pl.rcParams["font.size"] = 14

# name labels and colours for the various task types
task_types = np.loadtxt(
    args.type_file,
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
task_colors = pl.cm.ScalarMappable(cmap="tab20").to_rgba(
    np.linspace(0.0, 1.0, len(task_names))
)

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
    colors = [task_colors[task["type"]] for task in thread]
    ax.broken_barh(bar, (i - 0.4, 0.8), facecolors=colors, edgecolor="none")
    # optionally add labels
    if args.labels:
        # status text
        label = ""
        if nthread > 1:
            label += "thread {0} - ".format(i)
        tottime = sum([task["end"] - task["start"] for task in thread])
        tottime /= ttot_max - ttot_min
        label += "{0:.2f} \% load".format(tottime * 100.0)
        ax.text(
            0.5 * (ttot_min + ttot_max),
            i + 0.2,
            label,
            ha="center",
            bbox=dict(facecolor="white", alpha=0.9),
        )

# add empty blocks for the legend
for i in range(len(task_colors)):
    ax.plot([], [], color=task_colors[i], label=task_names[i].decode("utf-8"))

# add the legend and clean up the axes
ax.legend(loc="upper center", ncol=min(3, len(task_colors)))
ax.set_ylim(-1.0, nthread * 1.1)
ax.set_yticks([])

# get the total idle time and add information to the title
# alltime /= nthread * nproc
# ax.set_title("Total empty fraction: {0:.2f} \%".format((1.0 - alltime) * 100.0))
ax.set_xlabel("Simulation time (ticks)")

ax.xaxis.set_minor_locator(AutoMinorLocator())

# finalize and save the plot
pl.tight_layout()
pl.savefig(args.output, dpi=300, bbox_inches="tight")
