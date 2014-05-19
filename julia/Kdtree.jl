#!/usr/bin/env julia
# Fredrik Boulund 2014
# Simple octree implementation based on wikipedia article:
#  http://en.wikipedia.org/wiki/K-d_tree

module Kdtree

immutable Atom
	pos::(Float64,Float64,Float64)
	radius::Float64
end

type Node
	atom::Atom
	left
	right
end

function kdtree(atoms, depth=0)
	# Constructs a kdtree containing Atom types.
	# Input:
	#  atoms    An array of Atoms to construct tree from.
	#  depth    Current depth in the tree (determines what
	#		    axis to split on).
	# Returns:
	#  Node		The root node of the complete kd-tree.

	# If less than three atoms remain; special case
	if length(atoms) == 2
		return Node(atoms[1],
					kdtree(atoms[2:2], depth+1),
					nothing)
	elseif length(atoms) == 1
		return Node(atoms[1],
					nothing,
					nothing)
	end
		
	# Select axis based on depth so that axis cycles through all valid values
	axis = mod(depth,3)+1

	# Sort point list and choose median as pivot element
	atoms = atoms[sortperm([p.pos[axis] for p in atoms])]
	median = div(length(atoms),2)+1

	# Create node and construct subtrees
	return Node(atoms[median], 
				 kdtree(atoms[1:median-1], depth+1), 
				 kdtree(atoms[median+1:end], depth+1))
end


function sq_dist(a, b)
	# Computes squared distance between to points in Euclidian space.
	# Input:
	#  a	Point a, coordinates in tuple (x,y,z)
	#  b	Point b, coordinates in tuple (x,y,z)
	# Returns:
	#  dist 	Squared distance between the points.
	return (a[1]-b[1])^2 + (a[2]-b[2])^2 + (a[3]-b[3])^2
end


function nn_search(node, qpos, best_node=nothing, best_dist=Inf)
	# Searches a kd-tree for the nearest neighbor to a given point.
	# Input:
	#  node		  A Node in the tree. To query the entire tree,
	#             input the root node.
	#  qpos	 	  The position from which to find the nearest neighbor 
	#			  in the tree. A tuple of floats (x,y,z).
	#  best_node  Best Node found.
	#  best_dist  Shortest distance found so far. 
	# Returns:
	#  best_dist
	#  best_node
	#

	if node.left == node.right == nothing
		dist = sq_dist(node.atom.pos, qpos)
		if dist < best_dist
			best_dist = dist
			best_node = node
			return (best_dist, node)
		end
	elseif node.left != nothing
		best_dist, node = nn_search(node.left, qpos, best_node, best_dist)
	end
	if node.right != nothing
		best_dist, node = nn_search(node.right, qpos, best_node, best_dist)
	end
	return (best_dist, node)
end
			
	

end
