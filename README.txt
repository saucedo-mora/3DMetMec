The code is in a single file 3DMetMec, and the geometrical and mechanical data 
needs to be in the same folder. The code is authomatic, but the BC for the BDF file 
needs to be inserted manually in the code for each case.

geo_data.stl is a stl file with the geometrical external geometry.

mech_data.csv is a file with the mechanical properties (Young modulus)
in different parts of the component. The file is component by coordinates x, y, z of
points with known stiffness.

The units of the files geo_data.stl and mech_data.csv needs to be consistent, and the
name can not be changed.

The 3DMetMEc.jl needs to be open ina compiler and just run it.
