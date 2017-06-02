import sys

pathToConfig = sys.argv[1]

file_write = open(pathToConfig, 'w')

file_write.write('.00001\n')
file_write.write('t\n')
file_write.write('n\n')
file_write.write('1500\n')
file_write.write('0.000150\n')
file_write.write('-10\n')
file_write.write('1\n')
file_write.write('8300000000\n')
file_write.write('.35\n')
file_write.write('pR\n')
file_write.write('neo\n')
file_write.write('newton\n')
file_write.write(sys.argv[2]+'\n')

file_write.close()

