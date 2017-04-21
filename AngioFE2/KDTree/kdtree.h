#pragma once

#include <tuple>
#include <algorithm>
#include <type_traits>
#include <functional>
#include <vector>
#include <list>
#include <queue>
#include <stack>
#include <cassert>
#include <iterator>
#include <string>
#include <iostream>
#include <limits>

//helper methods that provide sensible defaults to the kd-tree
//computes the distance squared rather than the absolute distance
inline double ndim_distance(std::vector<double> & first, std::vector<double> & second)
{
    double rv =0;
    assert(first.size() == second.size());
    for(size_t i =0; i < first.size(); i++)
    {
        double diff = first[i] - second[i];
        rv += diff*diff;
    }
    return rv;
}

inline std::vector<double> null_accessor(std::vector<double> d)
{
    return d;
}

inline std::string vector_printer(std::vector<double> d)
{
    std::string res ="";
    for(size_t i=0; i < d.size(); i++)
    {
        res += std::to_string(d[i]);
        res += " ";
    }
    return res;
}

//this distance needs to relate to the distance computed by ndim_distance
//so in this case this also used the distance squared
inline double ndim_distance_to_plane(std::vector<double> & point,std::vector<double> & point_in_plane,std::vector<double> & normal)
{
    assert(point.size() ==  point_in_plane.size());
    assert(point.size() == normal.size());
    //consider checking that normal is a unit vector
    //computes dot((point - point_in_plane), normal)
    double rv = 0.0;
    for(size_t i=0; i < point.size(); i++)
    {
        rv += (point[i] - point_in_plane[i]) *normal[i];
    }
    return rv*rv;
}



//implements a K-d tree which automatically rebuilds itself as items are inserted into the structure
//the members of DIMR must have the less than operator defined on them, the default constructor must be defined and initialize each value to it's zero,
template <typename DIM, typename DIMR>
class KDTree
{
private:
	class KDNode
	{
	public:
		KDNode(int dimc, DIM d){ parent = nullptr; childLess = nullptr, childGreater = nullptr; dimchoice = dimc; depth = 1; dimensions = d; }
		int lessChildDepth()
		{
			if (childLess)
			{
				return childLess->depth;
			}
			return 0;
		}
		int greaterChildDepth()
		{
			if (childGreater)
			{
				return childGreater->depth;
			}
			return 0;
		}
		void ReComputeDepth()
		{
			depth = std::max(lessChildDepth(), greaterChildDepth()) + 1;
		}
		DIM dimensions;
		KDNode * parent;
		KDNode * childLess;
		KDNode * childGreater;
		int dimchoice;

		int depth; //the max of each subtree's depth + 1
	};
public:
	//distancef the distance function between two points
	//distancetoplane is: point, point in plane, normal
	//units a DIMR filled with whatever the unit value is in each dimension
	KDTree(std::function<DIMR(DIM)> accessor, std::function<double(DIMR &, DIMR &)> distancef, std::function<double(DIMR &, DIMR &, DIMR &)> distancetoplane, DIMR units);//ideally I want to remove dimensions later and have it be enforced in the constructor
	~KDTree();
	void insert(DIM item);//insert a single item
	void insert(std::list<DIM> items);//can have better performace  than individual insertions
	DIM nearest(DIM item);//returns the DIM closest to item, returns item if the tree is empty/ a best node cannot otherwise be found
	DIM nearestCondition(DIM item, std::function<bool(DIM)> condition);//returns the closest DIM to item meeting the condition returns item if there are no valid DIM's in the tree
	std::vector<DIM> within(DIM item, double dist, bool sorted = false);//returns all of the nodes that are closer to item than dist
	// an inplace rebuild of the tree only call this when the underlying nodes have moved eg mechanical step(s)
	void rebuild();

	//deletes the tree and allows the structure to be reused
	void clear();
	
#if !defined(NDEBUG) && defined(KDDEBUG)
        void verifyTree();
#else
	void verifyTree(){ ; }
#endif
        
        void PrintTree(std::function<std::string(DIMR)> printer);
private:
	KDNode * KDify(typename  std::vector< KDNode *>::iterator start, typename std::vector<KDNode*>::iterator end, int dimstart);
	void DFT(std::function<void(KDNode *)> func, KDNode * node);//depth first travelsal less, greater, node
        void DFTPre(std::function<void(KDNode *)> func, KDNode * node);//depth first travelsal less, greater, node
        void DFTPrePost(std::function<void(KDNode *)> prefunc,std::function<void(KDNode *)> postfunc, KDNode * node);
        void BFS(std::function<bool(KDNode *)> func, KDNode * node);
        KDNode * nearest_leaf_node(DIMR target,  KDNode * start);
        
	KDNode * root;
	std::function<DIMR(DIM)> _accessor;
	std::function<double(DIMR &, DIMR &)> _distancef;
	std::function<double(DIMR &, DIMR &, DIMR &)> _distancetoplane;
	size_t ndim;
	DIMR _units;

};

template <typename DIM, typename DIMR>
KDTree<DIM, DIMR>::KDTree(std::function<DIMR(DIM)> accessor, std::function<double(DIMR &, DIMR &)> distancef, std::function<double(DIMR &, DIMR &, DIMR &)> distancetoplane, DIMR units)
{
	//make sure this line is okay
	ndim = units.size();
	assert(ndim > 0);
	root = nullptr;
	_accessor = accessor;
	_distancef = distancef;
	_distancetoplane = distancetoplane;
	_units = units;
}

template <typename DIM, typename DIMR>
KDTree<DIM, DIMR>::~KDTree()
{
	//recursively delete nodes
	DFT([](KDNode * node)
	{
		delete node;
	}, root);
}
template <typename DIM, typename DIMR>
void KDTree<DIM, DIMR>::clear()
{
	//recursively delete nodes
	DFT([](KDNode * node)
	{
		delete node;
	}, root);

	root = nullptr;
}

template <typename DIM, typename DIMR>
void KDTree<DIM, DIMR>::DFT(std::function<void(KDNode *)> func, KDNode * node)//depth first travelsal less, greater, node
{
	if (node)
	{
		KDNode * child_less = node->childLess;
		KDNode * child_greater = node->childGreater;
		if (child_less)
		{
			DFT(func, child_less);
		}
		if (child_greater)
		{
			DFT(func, child_greater);
		}
		func(node);
	}
}
template <typename DIM, typename DIMR>
void KDTree<DIM, DIMR>::DFTPre(std::function<void(KDNode *)> func, KDNode * node)//depth first travelsal less, greater, node
{
	if (node)
	{
                func(node);
		KDNode * child_less = node->childLess;
		KDNode * child_greater = node->childGreater;
		if (child_less)
		{
			DFTPre(func, child_less);
		}
		if (child_greater)
		{
			DFTPre(func, child_greater);
		}
	}
}

template <typename DIM, typename DIMR>
void KDTree<DIM, DIMR>::DFTPrePost(std::function<void(KDNode *)> prefunc,std::function<void(KDNode *)> postfunc, KDNode * node)
{
    if (node)
	{
                prefunc(node);
		KDNode * child_less = node->childLess;
		KDNode * child_greater = node->childGreater;
		if (child_less)
		{
			DFTPrePost(prefunc, postfunc,child_less);
		}
		if (child_greater)
		{
			DFTPrePost(prefunc, postfunc, child_greater);
		}
		postfunc(node);
	}
}

template <typename DIM, typename DIMR>
void KDTree<DIM, DIMR>::BFS(std::function<bool(KDNode *)> func, KDNode * node)
{
    if(node)
    {
        std::queue<KDNode *> process;
        process.emplace(node);
        while(!process.empty())
        {
            if(func(process.front()))
            {
                //add the children of the node
                if(process.front()->childGreater)
                {
                    process.emplace(process.front()->childGreater);
                }
                if(process.front()->childLess)
                {
                    process.emplace(process.front()->childLess);
                }
                
            }
            process.pop();
        }
    }
}

template <typename DIM, typename DIMR>
void KDTree<DIM, DIMR>::insert(DIM item)
{
	//do the insertion then update the depths, then rebuild any trees that need it
	//the algorithm keeps track of the nodes that were traversed so their depths can be updated
	std::list<KDNode *> touched;
        
	DIMR goal = _accessor(item);
	//assert(goal.tuple_size() == ndim);
        assert(goal.size() == ndim);
	int cdim = 0;
	KDNode * current = root;

	if (!root)
	{
		KDNode *node = new KDNode(0, item);
		root = node;
                verifyTree();
                return;
	}
        bool found = false;
        while (!found)
        {
                DIMR ccut = _accessor(current->dimensions);//current cut
                touched.emplace_front(current);
                if(cdim != current->dimchoice)
                {
                    PrintTree(vector_printer);
                }
                assert(cdim == current->dimchoice);
                //check the nodes exist before comparing in that direction
                if (goal[cdim] < ccut[cdim] )
                {
                        if (current->childLess)
                        {
                                current = current->childLess;
                        }
                        else
                        {
                                //add the node
                                KDNode *node = new KDNode((cdim +1 ) % ndim, item);
                                node->parent = current;
                                current->childLess = node;
                                found = true;
                                assert(node->parent->dimchoice == cdim);
                        }
                }
                else
                {
                        if (current->childGreater)
                        {
                                current = current->childGreater;
                        }
                        else
                        {
                                //add the node
                                KDNode *node = new KDNode((cdim +1 ) % ndim, item);
                                node->parent = current;
                                current->childGreater = node;
                                found = true;
                                assert(node->parent->dimchoice == cdim);
                        }
                        
                }
                cdim = (cdim + 1) % ndim;
        }

        //now using the touched nodes to update the depth information

        //note this must be done in a bottom up manner
        auto iterd = touched.begin();
        while (iterd != touched.end())
        {
                (*iterd)->ReComputeDepth();
                ++iterd;
        }
        
        
        assert(current->depth >= 2);
        verifyTree();
        //where touched is actually needed rebuild any portion of the tree where the depths are too different
        //this must be done in a top down manner to limit the number of rebuilds done
        auto iterr = touched.rbegin();
        while (iterr != touched.rend())
        {
                //if the depths have a difference greater than or equal to 3 
                if (abs((*iterr)->lessChildDepth() - (*iterr)->greaterChildDepth()) >= 3)
                {
                        //rebuild and  return;
                        KDNode * parent = (*iterr)->parent;
                        KDNode ** pside;
                        //need to take into account if parent is nullptr
                        if(parent)
                        {
                            if (parent->childGreater == *iterr)
                            {
                                    pside = &parent->childGreater;
                            }
                            else
                            {
                                    pside = &parent->childLess;
                            }
                        }
                        else
                        {
                            pside = &root;
                        }
                        
                        cdim = (*iterr)->dimchoice;
                        std::vector<KDNode *> nodes;
                        DFT([&nodes](KDNode * node)
                        {
                                nodes.emplace_back(node);
                        }, *iterr);

                        //fix the pointers that weren't fixed by KDify
                        KDNode * subroot = KDify(nodes.begin(), nodes.end(), cdim);
                        subroot->parent = parent;
                        *pside = subroot;
                        DFT([](KDNode * node)
                        {
                                node->ReComputeDepth();
                        }, subroot);
                        //then recompute the upper depths
                        KDNode * upper = parent;
                        
                        if(subroot->childLess)
                        {
                            subroot->childLess->ReComputeDepth();
                        }
                        if(subroot->childGreater)
                        {
                            subroot->childGreater->ReComputeDepth();
                        }
                        subroot->ReComputeDepth();
                        
                        while (upper)
                        {
                                upper->ReComputeDepth();
                                upper = upper->parent;
                        }
                        verifyTree();
                        return;
                }
                ++iterr;
        }
        verifyTree();
}
//this will do at most the saame number of rebuilds as individual insertions
//it is likely to do fewer rebuilds 
template <typename DIM, typename DIMR>
void KDTree<DIM, DIMR>::insert(std::list<DIM> items)
{
    //do the insertion then update the depths, then rebuild any trees that need it
    //the algorithm keeps track of the nodes that were traversed so their depths can be updated
    std::vector<KDNode *> touched;
    auto iter = items.begin();
    while (iter !=  items.end())
    {
        DIMR goal = _accessor(*iter);
        //assert(goal.tuple_size() == ndim);
        assert(goal.size() == ndim);
        int cdim = 0;
        KDNode * current = root;

        if (!root)
        {
            KDNode *node = new KDNode(0, *iter);
            root = node;
            touched.emplace_back(root);
            verifyTree();
        }
        else
        {
            bool found = false;
            while (!found)
            {
                DIMR ccut = _accessor(current->dimensions); //current cut
                touched.emplace_back(current);
                if(cdim != current->dimchoice)
                {
                    PrintTree(vector_printer);
                }
                assert(cdim == current->dimchoice);
                //check the nodes exist before comparing in that direction
                if (goal[cdim] < ccut[cdim] )
                {
                    if (current->childLess)
                    {
                        current = current->childLess;
                    }
                    else
                    {
                        //add the node
                        KDNode *node = new KDNode((cdim +1 ) % ndim, *iter);
                        node->parent = current;
                        current->childLess = node;
                        found = true;
                        assert(node->parent->dimchoice == cdim);
                    }
                }
                else
                {
                    if (current->childGreater)
                    {
                        current = current->childGreater;
                    }
                    else
                    {
                        //add the node
                        KDNode *node = new KDNode((cdim +1 ) % ndim, *iter);
                        node->parent = current;
                        current->childGreater = node;
                        found = true;
                        assert(node->parent->dimchoice == cdim);
                    }

                }
                cdim = (cdim + 1) % ndim;
            }
        }
        ++iter;
    }

    // sort the nodes so the ones with lower depths are updated first
    std::stable_sort(touched.begin(),  touched.end(),  [](KDNode * first,  KDNode * second){
        return first->depth < second->depth;
    });
    
    //now using the touched nodes to update the depth information

    //note this must be done in a bottom up manner
    auto iterd = touched.begin();
    while (iterd != touched.end())
    {
        (*iterd)->ReComputeDepth();
        ++iterd;
    }

    verifyTree();
    //where touched is actually needed rebuild any portion of the tree where the depths are too different
    //this must be done in a top down manner to limit the number of rebuilds done
    auto iterr = touched.rbegin();
    while (iterr != touched.rend())
    {
        //if the depths have a difference greater than or equal to 3 
        if (abs((*iterr)->lessChildDepth() - (*iterr)->greaterChildDepth()) >= 3)
        {
            //rebuild and  return;
            KDNode * parent = (*iterr)->parent;
            KDNode ** pside;
            //need to take into account if parent is nullptr
            if(parent)
            {
                if (parent->childGreater == *iterr)
                {
                    pside = &parent->childGreater;
                }
                else
                {
                    pside = &parent->childLess;
                }
            }
            else
            {
                pside = &root;
            }

            int cdim = (*iterr)->dimchoice;
            std::vector<KDNode *> nodes;
            DFT([&nodes](KDNode * node)
                {
                    nodes.emplace_back(node);
                }, *iterr);

            //fix the pointers that weren't fixed by KDify
            KDNode * subroot = KDify(nodes.begin(), nodes.end(), cdim);
            subroot->parent = parent;
            *pside = subroot;
            DFT([](KDNode * node)
                {
                    node->ReComputeDepth();
                }, subroot);
            //then recompute the upper depths
            KDNode * upper = parent;

            if(subroot->childLess)
            {
                subroot->childLess->ReComputeDepth();
            }
            if(subroot->childGreater)
            {
                subroot->childGreater->ReComputeDepth();
            }
            subroot->ReComputeDepth();

            while (upper)
            {
                upper->ReComputeDepth();
                upper = upper->parent;
            }
            verifyTree();
            return;
        }
        ++iterr;
    }
    verifyTree();
}

template <typename DIM, typename DIMR>
void KDTree<DIM, DIMR>::rebuild()
{
    // avoids lots of deallocation and allocation
    std::vector<KDNode *> nodes;
    DFT([&nodes](KDNode * node)
        {
            nodes.emplace_back(node);
        }, root);
    KDNode * new_root = KDify(nodes.begin(),  nodes.end(),  0);
    root = new_root;
	if (new_root)
	{
		root->parent = nullptr;
		DFT([&nodes](KDNode * node)
		{
			node->ReComputeDepth();
		}, root);
	}
}


//returns a kd tree with starting dimension dimstart and given a start and end pointers 
template <typename DIM, typename DIMR>
typename KDTree<DIM, DIMR>::KDNode * KDTree<DIM, DIMR>::KDify(typename  std::vector<KDNode *>::iterator start, 
	typename std::vector< KDNode*>::iterator end, int dimstart)
{
	//make sure end is set correctly to not lose elements
	//check that start iterator does not move
	int size = std::distance(start, end);
	if (!size)
	{
		return nullptr;
	}
	int choice = size / 2;//check this is correct and does not need offset
	std::function<bool(KDNode *, KDNode * )> comparitor = [this, dimstart]
	(KDNode * first, KDNode * second)
        {
            return _accessor(first->dimensions)[dimstart] < _accessor(second->dimensions)[dimstart];
        };
	std::nth_element(start, start + choice, end, comparitor);
	auto median = start + choice;
	KDNode * medianNode = *median;
        //set the cdim this with recursion will set the cdim properly for all nodes
        medianNode->dimchoice = dimstart;
	//just process the halves that remain
	KDNode * lessChild = KDify(start, median, (dimstart +1) %ndim);
	KDNode * greaterChild = KDify(median + 1, end, (dimstart + 1) % ndim);
	if (lessChild)
	{
		lessChild->parent = medianNode;
		medianNode->childLess = lessChild;
	}
	else
	{
		medianNode->childLess = nullptr;
	}
	if (greaterChild)
	{
		greaterChild->parent = medianNode;
		medianNode->childGreater = greaterChild;
	}
	else
	{
		medianNode->childGreater = nullptr;
	}
	
	return medianNode;
}

template <typename DIM, typename DIMR>
typename KDTree<DIM, DIMR>::KDNode * KDTree<DIM, DIMR>::nearest_leaf_node(DIMR target,  KDNode * start)
{
    DIMR goal = target;
    assert(goal.size() == ndim);
    int cdim = start->dimchoice;
    KDNode * current = start;
    if (!root)
    {
        return nullptr;
    }
    while (true)
    {
        DIMR ccut = _accessor(current->dimensions); //current cut
        assert(cdim == current->dimchoice);
        //check the nodes exist before comparing in that direction
        if (goal[cdim] < ccut[cdim] )
        {
            if (current->childLess)
            {
                current = current->childLess;
            }
            else if(current->childGreater)
            {
                current = current->childGreater;
            }
        }
        else
        {
            if (current->childGreater)
            {
                current = current->childGreater;
            }
            else if(current->childLess)
            {
                current = current->childLess;
            }

        }
        cdim = (cdim + 1) % ndim;
        if((current->childLess ==  nullptr) && (current->childGreater ==  nullptr))
            return current;
    }
    return nullptr;
}

template <typename DIM, typename DIMR>
DIM KDTree<DIM, DIMR>::nearest(DIM item)
//returns a list of the closest DIM to item
{
    DIMR goal = _accessor(item);
    KDNode * startNode = nearest_leaf_node(goal,  root);
    KDNode * best_node = nullptr;
    if(startNode)
    {
        auto temp = _accessor(startNode->dimensions);
        double cbest = std::numeric_limits<double>::infinity();
        std::stack<KDNode *> n2p;
        std::stack<KDNode *> end_nodes;
        n2p.emplace(startNode);
        //end_nodes.push(root);
        // see wikipedia on nn in KDTree
        while(!n2p.empty())
        {
            KDNode * current = n2p.top();
            n2p.pop();
            temp = _accessor( current->dimensions);
            double d = _distancef(goal,  temp);
            if(d < cbest)
            {
                cbest = d;
                best_node = current;
            }
            
            if((!end_nodes.empty()) && (current == end_nodes.top()))
            {
                end_nodes.pop();
                continue;
            }
            
            //if points on the other side of spliting plane that could be closer
            // goto the leaf of that tree and continue
            DIMR zeros;
            zeros.resize(ndim,  0.0);
            int cdim = current->dimchoice;
            zeros[cdim] = _units[cdim];
            double d_pt_pl = _distancetoplane(goal,  temp,  zeros);
            auto bn = _accessor(best_node->dimensions);
            double d_g_b = _distancef(bn, goal);
            if(current->parent)
                n2p.emplace(current->parent);
            if(d_pt_pl <= d_g_b)
            {
                // there may be points on the other side of the plane
                if (goal[cdim] < temp[cdim] )
                {
                    if (current->childGreater)
                    {
                        KDNode * next_subtree = nearest_leaf_node(goal, current->childGreater);
                        if(next_subtree && (next_subtree !=  current))
                        {
                            n2p.emplace(next_subtree);
                            end_nodes.emplace(current);
                            continue;
                        }
                    } 
                }
                else
                {
                    if (current->childLess)
                    {
                        KDNode * next_subtree = nearest_leaf_node(goal, current->childLess);
                        if(next_subtree && (next_subtree != current))
                        {
                            n2p.emplace(next_subtree);
                            end_nodes.emplace(current);
                            continue;
                        }
                    }
                }
            }
        }
        
    }
    
    if(best_node)
    {
        return (best_node->dimensions);
    }
	assert(false);
    return item;
}

template <typename DIM, typename DIMR>
DIM KDTree<DIM, DIMR>::nearestCondition(DIM item, std::function<bool(DIM)> condition)
{
	DIMR goal = _accessor(item);
	KDNode * startNode = nearest_leaf_node(goal, root);
	KDNode * best_node = nullptr;
	if (startNode)
	{
		auto temp = _accessor(startNode->dimensions);
		double cbest = std::numeric_limits<double>::infinity();
		std::stack<KDNode *> n2p;
		std::stack<KDNode *> end_nodes;
		n2p.emplace(startNode);
		//end_nodes.push(root);
		// see wikipedia on nn in KDTree
		while (!n2p.empty())
		{
			KDNode * current = n2p.top();
			n2p.pop();
			temp = _accessor(current->dimensions);
			double d = _distancef(goal, temp);
			if ((d < cbest) && condition(current->dimensions))
			{
				cbest = d;
				best_node = current;
			}

			if ((!end_nodes.empty()) && (current == end_nodes.top()))
			{
				end_nodes.pop();
				continue;
			}

			//if points on the other side of spliting plane that could be closer
			// goto the leaf of that tree and continue
			DIMR zeros;
			zeros.resize(ndim, 0.0);
			int cdim = current->dimchoice;
			zeros[cdim] = _units[cdim];
			if (current->parent)
				n2p.emplace(current->parent);

			double d_pt_pl = _distancetoplane(goal, temp, zeros);
			double d_g_b;
			
			if (best_node)
			{
				DIMR bn = _accessor(best_node->dimensions);
				d_g_b = _distancef(bn, goal);
			}
			else
			{
				d_g_b = std::numeric_limits<double>::infinity();
			}

			if (d_pt_pl <= d_g_b)
			{
				// there may be points on the other side of the plane
				if (goal[cdim] < temp[cdim])
				{
					if (current->childGreater)
					{
						KDNode * next_subtree = nearest_leaf_node(goal, current->childGreater);
						if (next_subtree && (next_subtree != current))
						{
							n2p.emplace(next_subtree);
							end_nodes.emplace(current);
							continue;
						}
					}
				}
				else
				{
					if (current->childLess)
					{
						KDNode * next_subtree = nearest_leaf_node(goal, current->childLess);
						if (next_subtree && (next_subtree != current))
						{
							n2p.emplace(next_subtree);
							end_nodes.emplace(current);
							continue;
						}
					}
				}
			}
		}

	}
	if (best_node)
	{
		return best_node->dimensions;
	}
	assert(false);
	return item;
}

template <typename DIM, typename DIMR>
std::vector<DIM> KDTree<DIM, DIMR>::within(DIM item, double dist, bool ordered)
{
    std::vector<DIM> rv;
    DIMR goal = _accessor(item);
    KDNode * startNode = nearest_leaf_node(goal,  root);
    if(startNode)
    {
        DIMR temp = _accessor(startNode->dimensions);
        std::stack<KDNode *> n2p;
        std::stack<KDNode *> end_nodes;
        n2p.emplace(startNode);
        //end_nodes.push(root);
        // see wikipedia on nn in KDTree
        while(!n2p.empty())
        {
            KDNode * current = n2p.top();
            n2p.pop();
            temp = _accessor( current->dimensions);
            double d = _distancef(goal,  temp);
            if(d < dist)
            {
                rv.emplace_back(current->dimensions);
            }

            if((!end_nodes.empty()) && (current == end_nodes.top()))
            {
                end_nodes.pop();
                        continue;
            }

            //if points on the other side of spliting plane that could be closer
            // goto the leaf of that tree and continue
            DIMR zeros;
            zeros.resize(ndim,  0.0);
            int cdim = current->dimchoice;
            zeros[cdim] = _units[cdim];
            double d_pt_pl = _distancetoplane(goal,  temp,  zeros);
            if(current->parent)
            n2p.emplace(current->parent);
            if(d_pt_pl <= dist)
            {
                // there may be points on the other side of the plane
                if (goal[cdim] < temp[cdim] )
                {
                    if (current->childGreater)
                    {
                        KDNode * next_subtree = nearest_leaf_node(goal, current->childGreater);
                        if(next_subtree && (next_subtree !=  current))
                        {
                            n2p.emplace(next_subtree);
                            end_nodes.emplace(current);
                        continue;
                        }
                    } 
                }
                else
                {
                    if (current->childLess)
                    {
                        KDNode * next_subtree = nearest_leaf_node(goal, current->childLess);
                        if(next_subtree && (next_subtree != current))
                        {
                            n2p.emplace(next_subtree);
                            end_nodes.emplace(current);
                        continue;
                        }
                    }
                }
            }
        }

    }

	if (ordered)
	{
		std::sort(rv.begin(), rv.end(), [&goal, this](DIM one, DIM two)
		{
			auto temp1 = _accessor(one);
			double d1 = _distancef(goal, temp1);

			auto temp2 = _accessor(two);
			double d2 = _distancef(goal, temp2);
			return d1 < d2;
		});
	}
    return rv;
}

#if !defined(NDEBUG) && defined(KDDEBUG)
template <typename DIM, typename DIMR>
void KDTree<DIM, DIMR>::verifyTree()
{
    if(!root)
    {
        return;
    }
    std::function<void(KDNode *)> func = [this]( KDNode * current){
        if(current)
        {
            if(current->parent)
            {
                assert(((current->parent->dimchoice + 1) % ndim) ==  current->dimchoice);
                assert((current->parent->childLess == current) || (current->parent->childGreater == current));
            }
            if(current->childLess == nullptr && current->childGreater == nullptr)
            {
                assert(current->depth == 1);
            }
            assert(current->depth == (std::max(current->greaterChildDepth(), current->lessChildDepth()) +1));
            assert(current->childGreater != current);
            assert(current->childLess != current);
            
            //check that the tree is actually split by the current node
            std::function<void(KDNode *)> verSplitLess = [this, current]( KDNode * c)
            {
                int cdim = current->dimchoice;
                DIMR pcut = _accessor(current->dimensions);//parent cut
                DIMR ccut = _accessor(c->dimensions);//current cut
                //check the nodes exist before comparing in that direction
                assert (ccut[cdim] <= pcut[cdim]);
            };
            if(current->childLess)
            {
                DFT(verSplitLess, current->childLess);
            }
            
            //check that the tree is actually split by the current node
            std::function<void(KDNode *)> verSplitGreater = [this, current]( KDNode * c)
            {
                int cdim = current->dimchoice;
                DIMR pcut = _accessor(current->dimensions);//parent cut
                DIMR ccut = _accessor(c->dimensions);//current cut
                //check the nodes exist before comparing in that direction
                assert (ccut[cdim] >= pcut[cdim]);
            };
            if(current->childGreater)
            {
                DFT(verSplitGreater, current->childGreater);
            }
        }
    };
    DFT(func, root);    
}
#endif
//currently the order is in the order of DFT
template <typename DIM, typename DIMR>
void KDTree<DIM, DIMR>::PrintTree(std::function<std::string(DIMR)> printer)
{
    if(!root)
    {
        return;
    }
    int d=0;
    std::function<void(KDNode *)> prefunc = [this,&d, &printer]( KDNode * current){
        if(current)
        {
            std::string dir = "b";//b for base/root
            if(current->parent)
            {
                if(current->parent->childGreater == current)
                {
                    dir = "r";
                }
                else if(current->parent->childLess == current)
                {
                    dir = "l";
                }
                else
                    assert(false);
            }
            for(int i =0; i < d;i++)
            {
                std::cout << "  ";
            }
            std::cout << "de: " << current-> depth << " di: " << current->dimchoice << " data: " << printer(_accessor(current->dimensions)) << " " << dir << std::endl;
            d++;
        }
    };
    std::function<void(KDNode *)> postfunc = [this,&d, &printer]( KDNode * current){
        if(current)
        {
            d--;
        }
    };
    DFTPrePost(prefunc, postfunc, root); 
}
