# To create a PNG of the graph, use
#   dot -Tpng OUTPUT > NAME.png

import argparse

argparser = argparse.ArgumentParser(
    "Create a Grapviz graph representation of a task-based calculation."
)
argparser.add_argument(
    "--file",
    "-f",
    action="store",
    help="Name of the task-based log output file that holds task information.",
    required=True,
)
argparser.add_argument(
    "--output",
    "-o",
    action="store",
    help="Name of the output file.",
    required=True,
)
args = argparser.parse_args()


def trim_label(label):
    if len(label) > 10:
        return label[2:]
    else:
        return label[1:]


dfile = open(args.file, "r")

resources = {}
tasks = {}
readlinks = []
writelinks = []
tasklinks = []
children = []
for line in dfile.readlines():
    data = line.split()
    if data[0] == "resource":
        resources[int(data[1])] = trim_label(data[2])
    elif data[0] == "task":
        tasks[int(data[1])] = trim_label(data[2])
    elif data[0] == "readlink":
        readlinks.append([int(data[1]), int(data[2])])
    elif data[0] == "writelink":
        writelinks.append([int(data[1]), int(data[2])])
    elif data[0] == "tasklink":
        tasklinks.append([int(data[1]), int(data[2])])
    elif data[0] == "child":
        children.append([int(data[1]), int(data[2])])

lines = []
nodes = []
for link in tasklinks:
    if link[0] in tasks and link[1] in tasks:
        line = tasks[link[0]] + "->" + tasks[link[1]] + ";"
        if not line in lines:
            lines.append(line)
for link in readlinks:
    if link[0] in tasks and link[1] in resources:
        res = "res" + resources[link[1]]
        if not res in nodes:
            nodes.append(res)
        line = tasks[link[0]] + "->" + res + '[color="blue"];'
        if not line in lines:
            lines.append(line)
for link in writelinks:
    if link[0] in tasks and link[1] in resources:
        res = "res" + resources[link[1]]
        if not res in nodes:
            nodes.append(res)
        line = tasks[link[0]] + "->" + res + '[color="red"];'
        if not line in lines:
            lines.append(line)
for child in children:
    if child[0] in resources and child[1] in resources:
        res1 = "res" + resources[child[0]]
        res2 = "res" + resources[child[1]]
        line = res1 + "->" + res2 + '[color="green"];'
        if not line in lines:
            lines.append(line)

ofile = open(args.output, "w")
ofile.write("digraph task_flow{\n")
for node in nodes:
    ofile.write(node + '[label="' + node[3:] + '" color="blue"];\n')
for line in lines:
    ofile.write(line + "\n")
ofile.write("}")
