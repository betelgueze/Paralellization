/**
 * Serial and parallel alogrithms for detecting
 * edge connectivity in undirected graphs.
 *
 * Author: Jiri Hon, xhonji01@stud.fit.vutbr.cz
 * Date: 2014-12-09
 */

#include <ogdf/basic/Graph.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/graphalg/MinimumCut.h>
#include <ogdf/basic/BinaryHeap.h>
#include <list>
#include <omp.h>
#include <unistd.h>

using namespace ogdf;
using namespace std;


typedef struct flow
{// Structure for saving flow through edge in both directions
	int s_t;
	int t_s;
} flow_t;


/**
 * Get virtual memory usage of current process.
 * NOTE: This works only for Linux and does not give precise results.
 * @author Don Wakefield @see http://stackoverflow.com/a/671389
 * @param vm_usage
 * @param resident_set
 */
void process_mem_usage(double &vm_usage, double &resident_set)
{
	vm_usage     = 0.0;
	resident_set = 0.0;

	// 'file' stat seems to give the most reliable results
	ifstream stat_stream("/proc/self/stat",ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
					>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
					>> utime >> stime >> cutime >> cstime >> priority >> nice
					>> O >> itrealvalue >> starttime >> vsize >> rss;
					// don't care about the rest

	stat_stream.close();

	// in case x86-64 is configured to use 2MB pages
	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
	vm_usage     = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}

/**
 * Get the other node connected by given edge
 * @param u First node
 * @param e Edge
 * @return Second node connected by edge
 */
static inline node other_node(node u, edge e)
{
	if (e->source() == u)
		return e->target();
	else
		return e->source();
}

/**
 * @brief Print graph adjacency list
 * @param G Graph
 */
static inline void print_adj_list(Graph &G)
{
	node u; edge e;
	forall_nodes(u, G)
	{
		cout << u->index() << ": ";
		forall_adj_edges(e, u)
			cout << other_node(u, e) << "[" << e->index() << "] ";
		cout << endl;
	}
}

/**
 * Print graph summary info
 * @param G
 */
static inline void print_graph_stats(Graph &G)
{
	cout << "|V| = " << G.numberOfNodes()
		  << " |E| = " << G.numberOfEdges() << endl;
}

/**
 * Get flow between u,v
 * @param F Flow for every edge
 * @param u First node
 * @param e Edge leading to v
 * @return Flow in given direction
 */
static inline int get_f(EdgeArray<flow_t> &F, node u, edge e)
{
	if (u == e->source())
		return F[e].s_t;
	else
		return F[e].t_s;
}

/**
 * Add flow between u,v
 * @param F Flow for every edge
 * @param u First node
 * @param e Edge leading to v
 * @param f Flow to be added
 */
static inline void add_f(EdgeArray<flow_t> &F, node u, edge e, int f)
{
	if (u == e->source())
		F[e].s_t += f;
	else
		F[e].t_s += f;
}

/**
 * Breadth First Search as a part of Edmonds-Karp algorithm.
 * @param G Graph
 * @param s Source node
 * @param t Target node
 * @param P Reverse path from target to source
 * @param F Flow for every edge
 * @return Flow increase
 */
int EdmondsKarpBFS(
	Graph &G, node s, node t, NodeArray<edge> &P, EdgeArray<flow_t> &F)
{
	node u, v;
	edge e;
	int f;

	forall_nodes(u, G)
	// Initialize backward path array
		P[u] = NULL;

	list<node> Q;
	Q.push_back(s);

	while (!Q.empty())
	{
		u = Q.front(); Q.pop_front();

		forall_adj_edges(e, u)
		{
			v = other_node(u, e);
			if (v == s)
				continue;

			f = get_f(F, u, e);

			if (1 - f > 0 && P[v] == NULL)
			{
				P[v] = e;
				if (v != t)
					Q.push_back(v);
				else
					return 1;
			}
		}
	}
	return 0;
}

/**
 * Edmonds-Karp algorithm implementation for simple graphs
 * @param G Graph
 * @param s Source node
 * @param t Target node
 * @return Maximal flow from s to t
 */
int EdmondsKarp(Graph &G, node s, node t)
{
	node u, v;
	edge e;
	int m, f = 0;

	// Array indicating flow through each edge
	flow_t def_flow = {0, 0};
	EdgeArray<flow_t> F(G, def_flow);
	NodeArray<edge> P(G); // Backward path from BFS

	while (true)
	{
		m = EdmondsKarpBFS(G, s, t, P, F);
		if (m == 0)
			break;
		f = f + m;

		v = t;
		while (v != s)
		{
			e = P[v];
			u = other_node(v, e);
			add_f(F, u, e, m);
			add_f(F, v, e, -m);
			v = u;
		}
	}
	return f;
}

/**
 * Serial implementation of Edge connectivity algorithm using
 * N calls to Edmonds-Karp algorithm.
 * @param G Graph
 * @return Edge connectivity
 */
int MinMaxFlowSerial(Graph &G)
{
	node s = G.firstNode();
	int min = G.numberOfEdges();
	node t;
	int m;

	forall_nodes(t, G)
	{
		if (t == s)
			continue;

		m = EdmondsKarp(G, s, t);
		if (m == 0)
			break;
		if (m < min)
			min = m;
	}
	return min;
}

/**
 * Parallel implementation of Edge connectivity algorithm using
 * N calls to Edmonds-Karp algorithm.
 * @param G Graph
 * @return Edge connectivity
 */
int MinMaxFlowParallel(Graph &G)
{
	node s = G.firstNode();
	int min = G.numberOfEdges();

	// Preprocess nodes to make parallel-for possible
	node nodes[G.numberOfNodes()];
	node u;
	int i = 0;
	forall_nodes(u, G)
		nodes[i++] = u;

	#pragma omp parallel default(shared) firstprivate(i,s) reduction(min:min)
	{
		node t;
		int m;

		#pragma omp for schedule(dynamic, 10)
		for(i = 0; i < G.numberOfNodes(); ++i)
		{
			t = nodes[i];
			if (t == s)
				continue;

			m = EdmondsKarp(G, s, t);
			if (m < min)
				min = m;
		}
	}
	return min;
}

/**
 * Print edge lists
 * @param E Edge lists
 */
static inline void print_edge_lists(Array<list<edge> > &E)
{
	for (int i = 0; i < E.size(); i++)
	{
		cout << i << ": ";
		for (list<edge>::const_iterator it = E[i].begin(), end = E[i].end(); it != end; ++it)
			cout << (*it)->index() << " ";
		cout << endl;
	}
}

/**
 * Reinitiali all edge lists
 * @param E Edge lists
 */
static inline void clear_edge_lists(Array<list<edge> > &E)
{
	for (int i = 0; i < E.size(); i++)
		E[i].clear();
}

/**
 * Contract nodes connected with given edge and also
 * remove possible self loops.
 * @param G Graph
 * @param e Contracted edge
 */
static inline void contract_nodes(Graph &G, edge e)
{
	node s = e->source();
	G.contract(e);

	forall_adj_edges(e, s)
	{// Delete self loops
		if (e->isSelfLoop())
			G.hideEdge(e);
	}
}

/**
 * Get minimal degree of graph
 * @param G Graph
 * @return Minimal degree
 */
static inline int min_degree(Graph &G)
{
	node u;
	int min = G.numberOfEdges();

	forall_nodes(u, G)
	{
		if (u->degree() < min)
			min = u->degree();
	}
	return min;
}

/**
 * Implementation of FOREST algorithm of Nagamochi & Ibaraki
 * @param G Graph
 * @return Last scanned edge
 */
edge NagamochiIbarakiForest(Graph &G)
{
	node u, v; edge e, le;
	EdgeArray<bool> es(G, false);

	BinaryHeap<node> Q(G.numberOfNodes());
	NodeArray<const BinaryHeap<node>::Element *> R(G, NULL);

	forall_nodes(u, G)
		R[u] = &(Q.insert(u, 0.0));

	while (!Q.empty())
	{// Unscanned node exists
		u = Q.extractMin();

		forall_adj_edges(e, u)
		{
			if (!es[e])
			{// Unscanned edge exists
				v = other_node(u, e);
				Q.decPriority(*(R[v]), R[v]->getPriority() - 1);
				es[e] = true;
				le = e;
			}
		}
	}
	return le;
}

/**
 * Serial implementation of edge connectivity algorithm
 * by Nagamochi & Ibaraki
 * @param G Graph
 * @return Edge connectivity
 */
int NagamochiIbarakiSerial(Graph &G)
{
	Graph Gc = G;

	int k = Gc.numberOfEdges();
	edge le;

	while (Gc.numberOfNodes() > 2)
	{
		le = NagamochiIbarakiForest(Gc);
		k = min(k, min_degree(Gc));
		contract_nodes(Gc, le);
	}
	// Dirty hack to bypass a possible bug with hiding edges that
	// was observed when there are only two nodes left in Gc
	Gc.restoreAllEdges();
	return min(k, min_degree(Gc));
}

/**
 * Print brief help and usage suggestions
 */
void help()
{
	cout << "Serial and parallel alogrithms for detecting edge connectivity" << endl
		  << "in undirected graphs." << endl
		  << "Authors: Jiri Hon and Martin Risa" << endl
		  << "USAGE: ./econn -i input.gml [-a algorithm] [-t threads]" << endl
		  << "  -i Path to input GML file." << endl
		  << "  -a Used algorithm. Possible values: minmaxflow, nagamochiibaraki." << endl
		  << "     Default is nagamochiibaraki." << endl
		  << "  -t Number of threads to be used, but only minmaxflow algorithm" << endl
		  << "     supports parallel computation." << endl
		  << "  -h Print this help message." << endl
		  << "OUTPUT FORMAT:" << endl
		  << "[#threads] [#nodes] [duration] [exp.duration] [edge conn.] [memory used]" << endl;
}


/**
 * Application start
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv)
{
	Graph G;
	int opt, econn, refeconn, n, m;
	int num_thread = 1;
	string alg = "nagamochiibaraki";
	string path;
	double start, duration, exp_duration, rss, vms;

	/**
	 * Parse command-line options
	 */
	while ((opt = getopt(argc, argv, "a:i:t:h")) != -1)
	{
		switch (opt)
		{
			case 'a':
				alg = optarg;
				break;
			case 't':
				num_thread = stoi(optarg);
				omp_set_num_threads(num_thread);
				break;
			case 'i':
				path = optarg;
				break;
			case 'h':
			default:
				help();
				return 0;
				break;
		}
	}
	if (path.empty())
	{
		help();
		return 1;
	}

	if (!GraphIO::readGML(G, path))
		return 2;
	n = G.numberOfNodes();
	m = G.numberOfEdges();

	/**
	 * Run chosen algorithm
	 */
	start = omp_get_wtime();

	if (alg == "minmaxflow")
	{
		exp_duration = ((double) n*n * m*m) / 0.2e14;
		if (num_thread > 1)
			econn = MinMaxFlowParallel(G);
		else
			econn = MinMaxFlowSerial(G);
	}
	else if (alg == "nagamochiibaraki")
	{
		exp_duration = ((double) n*m) / 1e7;
		econn = NagamochiIbarakiSerial(G);
	}
	duration = omp_get_wtime() - start;

	if (num_thread > 1)
		exp_duration /= num_thread;

	process_mem_usage(vms, rss);

	/**
	 * Compute reference edge connectivity for check
	 */
	start = omp_get_wtime();
	EdgeArray<double> w(G, 1);
	MinCut mc = MinCut(G, w);
	refeconn = mc.minimumCut();
	if (refeconn < 0)
		refeconn = 0;

	/**
	 * Generate output
	 */
	if (econn != refeconn)
	// Invalid edge connectivity
		cout << "FAIL " << econn << " != " << refeconn << endl;
	else
	// Computed edge connectivity match reference
		cout << num_thread << " " << n << " " << duration << " "
			  << exp_duration << " " << econn << " " << rss << endl;
	return 0;
}
