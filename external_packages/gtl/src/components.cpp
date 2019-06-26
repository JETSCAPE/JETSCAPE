/* This software is distributed under the GNU Lesser General Public License */
//==========================================================================
//
//   components.cpp
//
//==========================================================================
// $Id: components.cpp,v 1.5 2001/11/07 13:58:09 pick Exp $

#include <GTL/components.h>

#ifdef __GTL_MSVCC
#   ifdef _DEBUG
#	ifndef SEARCH_MEMORY_LEAKS_ENABLED
#	error SEARCH NOT ENABLED
#	endif
#	define new DEBUG_NEW
#	undef THIS_FILE
	static char THIS_FILE[] = __FILE__;
#   endif   // _DEBUG
#endif	// __GTL_MSVCC

__GTL_BEGIN_NAMESPACE

components::components () : dfs ()
{
    scan_whole_graph (true);
    num_of_components = 0;
}

void components::reset () 
{ 
    dfs::reset ();
    comp.erase (comp.begin(), comp.end());
    num_of_components = 0;
}

int components::check (graph& G) 
{
    return G.is_undirected() && whole_graph && 
	dfs::check (G) == GTL_OK ? GTL_OK : GTL_ERROR;
}
    

//--------------------------------------------------------------------------
//   Handler
//--------------------------------------------------------------------------


void components::new_start_handler (graph& G, node& st) 
{
    li = comp.insert (comp.end(), 
	pair<list<node>,list<edge> > (list<node> (), list<edge> ()));
    (*li).first.push_back (st);
    ++num_of_components;
}

void components::before_recursive_call_handler (graph& G, edge& e, node& n)
{
    (*li).first.push_back (n);
    // (*li).second.push_back (e);    
}


void components::old_adj_node_handler (graph& G, edge& e, node& n) 
{
    node curr = n.opposite (e);

    //
    // Store backedges at lower endpoint
    //

    if (dfs_num (curr) > dfs_num (n)) { 
	(*li).second.push_back (e);    
    }
}


__GTL_END_NAMESPACE

//--------------------------------------------------------------------------
//   end of file
//--------------------------------------------------------------------------
