/* This software is distributed under the GNU Lesser General Public License */
//==========================================================================
//
//   bellman_ford_test.cpp
//
//==========================================================================
// $Id: bellman_ford_test.cpp,v 1.2 2002/11/07 13:38:37 raitner Exp $

#include <iostream>

#include <GTL/graph.h>
#include <GTL/bellman_ford.h>
#include <GTL/edge_map.h>

#ifdef __GTL_MSVCC
#   ifdef _DEBUG
#	define new DEBUG_NEW
#	undef THIS_FILE
	static char THIS_FILE[] = __FILE__;
#   endif   // _DEBUG
#endif	// __GTL_MSVCC

int main (int argc, char* argv[])
{
    graph G;
    G.make_directed();
    node n1 = G.new_node();
    node n2 = G.new_node();
    node n3 = G.new_node();
    node n4 = G.new_node();
    node n5 = G.new_node();
    node n6 = G.new_node();

    edge e1_2 = G.new_edge(n1, n2);
    edge e2_3 = G.new_edge(n2, n3);
    edge e1_6 = G.new_edge(n1, n6);
    edge e3_6 = G.new_edge(n3, n6);
    edge e4_3 = G.new_edge(n4, n3);
    edge e6_4 = G.new_edge(n6, n4);
    edge e6_5 = G.new_edge(n6, n5);
    edge e5_4 = G.new_edge(n5, n4);
    edge e5_1 = G.new_edge(n5, n1);

    edge_map<double> w(G);
    w[e1_2] = 1;
    w[e2_3] = 2;
    w[e1_6] = 8;
    w[e3_6] = 3;
    w[e4_3] = 2;
    w[e6_4] = 1;
    w[e6_5] = 3;
    w[e5_4] = 10;
    w[e5_1] = 2;

    node_map<double> result(G);
    result[n1] = 0;
    result[n2] = 1;
    result[n3] = 3;
    result[n4] = 7;
    result[n5] = 9;
    result[n6] = 6;

    node_map<node> preds(G);
    preds[n1] = node();
    preds[n2] = n1;
    preds[n3] = n2;
    preds[n4] = n6;
    preds[n5] = n6;
    preds[n6] = n3;

    bellman_ford B;
    B.store_preds(true);
    B.source(n1);
    B.weights(w);

    G.save("test.gml");
    
    cout << "Bellman Ford with positive edge weights" << endl;
    
    if (B.check(G) == algorithm::GTL_OK) 
    {
	cout << "check OK" << endl; 
    } 
    else 
    {
	cout << "check FAILED" << endl;
	exit(1);
    }

    if (B.run(G) == algorithm::GTL_OK) 
    {
	cout << "run OK" << endl; 
    } 
    else 
    {
	cout << "run FAILED" << endl;
	exit(1);
    }    

    graph::node_iterator it, end;
    
    for (it = G.nodes_begin(), end = G.nodes_end(); it != end; ++it)
    {
	if(result[*it] != B.distance(*it))
	{
	    cout << "Distance for node " << *it << " FAILED" << endl;
	    exit(1);
	}
    }

    cout << "Distances OK" << endl;

    for (it = G.nodes_begin(), end = G.nodes_end(); it != end; ++it)
    {
	if(preds[*it] != B.predecessor_node(*it))
	{
	    cout << "Predecessor for node " << *it << " FAILED" << endl;
	    exit(1);
	}
    }

    cout << "Predecessors OK" << endl;

    if (B.negative_cycle()) 
    {
	cout << "Negative Cycle FAILED" << endl;
    }
    else 
    {
	cout << "Negative Cycle OK" << endl;	
    }

    //----------------------------------------------------------------------
    //   negative edge weights
    //----------------------------------------------------------------------

    B.reset();

    cout << "Bellman Ford with some negative edge weights" << endl;

    w[e1_2] = 1;
    w[e2_3] = 2;
    w[e1_6] = -3;
    w[e3_6] = 3;
    w[e4_3] = 2;
    w[e6_4] = 1;
    w[e6_5] = 3;
    w[e5_4] = 10;
    w[e5_1] = 2;

    B.weights(w);

    result[n1] = 0;
    result[n2] = 1;
    result[n3] = 0;
    result[n4] = -2;
    result[n5] = 0;
    result[n6] = -3;

    preds[n1] = node();
    preds[n2] = n1;
    preds[n3] = n4;
    preds[n4] = n6;
    preds[n5] = n6;
    preds[n6] = n1;

     G.save("test2.gml");
    

    if (B.check(G) == algorithm::GTL_OK) 
    {
	cout << "check OK" << endl; 
    } 
    else 
    {
	cout << "check FAILED" << endl;
	exit(1);
    }

    if (B.run(G) == algorithm::GTL_OK) 
    {
	cout << "run OK" << endl; 
    } 
    else 
    {
	cout << "run FAILED" << endl;
	exit(1);
    }
    
    for (it = G.nodes_begin(), end = G.nodes_end(); it != end; ++it)
    {
	if(result[*it] != B.distance(*it))
	{
	    cout << "Distance for node " << *it << " FAILED" << endl;
	    exit(1);
	}
    }

    cout << "Distances OK" << endl;

    for (it = G.nodes_begin(), end = G.nodes_end(); it != end; ++it)
    {
	if(preds[*it] != B.predecessor_node(*it))
	{
	    cout << "Predecessor for node " << *it << " FAILED" << endl;
	    exit(1);
	}
    }

    cout << "Predecessors OK" << endl;

    if (B.negative_cycle()) 
    {
	cout << "Negative Cycle FAILED" << endl;
    }
    else 
    {
	cout << "Negative Cycle OK" << endl;	
    }

    //----------------------------------------------------------------------
    //   negative cycle
    //----------------------------------------------------------------------

    B.reset();

    cout << "Bellman Ford with negative cycle" << endl;

    w[e1_2] = 1;
    w[e2_3] = 2;
    w[e1_6] = -8;
    w[e3_6] = 3;
    w[e4_3] = 2;
    w[e6_4] = 1;
    w[e6_5] = 3;
    w[e5_4] = 10;
    w[e5_1] = 2;

    B.weights(w);

    if (B.check(G) == algorithm::GTL_OK) 
    {
	cout << "check OK" << endl; 
    } 
    else 
    {
	cout << "check FAILED" << endl;
	exit(1);
    }

    if (B.run(G) == algorithm::GTL_OK) 
    {
	cout << "run OK" << endl; 
    } 
    else 
    {
	cout << "run FAILED" << endl;
	exit(1);
    }
    
    if (B.negative_cycle()) 
    {
	cout << "Negative Cycle OK" << endl;	
    }
    else 
    {
	cout << "Negative Cycle FAILED" << endl;
    }
}
