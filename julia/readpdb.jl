#!/usr/bin/env julia
# Fredrik Boulund 2014-05-12
# First attempts at Julia 

require("ArgParse")
using ArgParse

settings = ArgParseSettings()
@add_arg_table settings begin
	"file"
		help = "A PDB file from which to read atom coordinates."
		required = true
	"--opt"
		help = "an optional argument"
end
args = parse_args(ARGS, settings)

immutable Atom
	x::Float64
	y::Float64
	z::Float64
	r::Float64
end

global atoms
count = 1
file = open(args["file"])
for line=eachline(file)
	#println(chomp(line))
	if beginswith(line, "ATOM")
		_, _, name, residue, _, _, x, y, z, _, _, _, = split(line)
		x, y, z = parsefloat(x), parsefloat(y), parsefloat(z)
		if count == 1
			atoms = [Atom(x,y,z,1)]
		else
			atoms = [atoms Atom(x,y,z,1)]
		end
		count += 1
	end
end

println("Read $count atoms from '$(args["file"])'")
