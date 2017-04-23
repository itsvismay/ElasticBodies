from pymeshfix import _meshfix

# Read mesh from infile and output cleaned mesh to outfile
file = "/home/vismay/ElasticBodies/3dUnion"
_meshfix.CleanFromFile(file+"/BrokenDesign/out.off", file+"/BrokenDesign/out.off")