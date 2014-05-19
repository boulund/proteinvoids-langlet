#!/usr/bin/env julia
# Fredrik Boulund 2014
# Fredrik Erik Viktor Protein Voids
# Reimplementation of 2009 Bachelor thesis project in Julia

push!(LOAD_PATH, ".") # Include current dir when searching for modules.

using ArgParse
import Kdtree   # Custom Kd-tree implementation for Atoms
import ReadPDB 

settings = ArgParseSettings()
@add_arg_table settings begin
	"file"
		help = "A PDB file from which to read atom coordinates."
		required = true
	"--types"
		help = "A tab separated file with atom information (Residue, Atom, Group)"
	"--print"
		help = "Print the data structures and results after each step. Not recommended."
		action = :store_true
	"--opt"
		help = "an optional argument"
end
args = parse_args(ARGS, settings)


atoms = ReadPDB.readPDB(args["file"])
if args["print"] println(atoms) end
tree = Kdtree.kdtree(atoms)
if args["print"] println(tree) end
dist, node = Kdtree.nn_search(tree, (0,0,0))
if args["print"] println("$(sqrt(dist)) $node") end
