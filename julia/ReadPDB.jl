#!/usr/bin/env julia
# Fredrik Boulund 2014

module ReadPDB
require("Kdtree")
import Kdtree
Atom = Kdtree.Atom

function readPDB(filename)
	atoms = []
	count = 1
	file = open(filename)
	for line=eachline(file)
		if beginswith(line, "ATOM") 
			_, _, name, residue, _, _, x, y, z, _, _, _, = split(line)
			x, y, z = parsefloat(x), parsefloat(y), parsefloat(z)
			if count == 1
				atoms = [Atom((x,y,z),1)]
			else
				atoms = [atoms Atom((x,y,z),1)]
			end
			count += 1
		end
	end
	println("Read $count atoms from '$filename'")
	return atoms
end

end
