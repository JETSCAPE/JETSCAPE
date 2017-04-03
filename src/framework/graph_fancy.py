#!/usr/local/bin/python

from graph_tool.all import *
from pylab import *

print 'Graph-tool test program'

g = load_graph("~/JetScape/framework.cmake/my_test.graphml")

print 'Graph is directed = %s '%(g.is_directed())
print 'my_test.graphml loaded!'

print 'Graph Veretx/Node and Edge properties:'
g.list_properties()

# --------------------------
# Setup for graphViz output usage from graph-tool (best for directed graphs "dot")

pT=g.edge_properties["pT"]
E=g.edge_properties["E"]

t=g.vertex_properties["t"]
label=g.vertex_properties["label"]

# easy here to change the color scales ... (how to do directly in dot format not clear, certainly not easy via a script)
E.a=1/(E.a/10)

gdict={'rankdir' : 'LR'} #or top botton 'TB'
edict={'dir' : 'forward', 'arrowhead' : 'normal','arrowsize' : '1.0','label' : pT}
vdict={'shape' : "plain", 'label' : label}

graphviz_draw(g,layout="dot",gprops=gdict, vprops=vdict, eprops=edict, ecmap=matplotlib.cm.autumn, penwidth=E, ecolor=E, vcolor="#ffffff", output="my_test_fancy.pdf")

# vcolor=t,vcmap=matplotlib.cm.gist_gray
# vnorm=True,

# --------------------------

# ##! python

#for e in g.edges():
#    print(e)

#edges = g.get_edges()

#for e in g.edges():
#    print '%s and pT= %s'%(e,g.edge_properties.pT[e])

#deg = g.degree_property_map("in")

# --------------------------

#graph_draw(g,edge_pen_width=pT,output="test-draw.pdf")

#pos=planar_layout(g)

#graph_draw(g, pos=pos,output="test-draw.pdf")

#pos = radial_tree_layout(g, g.vertex(0))
#graph_draw(g, pos=pos, output="graph-draw-radial.pdf")

#g.save("test.dot")

#tree = min_spanning_tree(g)
#graph_draw(g, edge_color=tree, output="min_tree.pdf")
