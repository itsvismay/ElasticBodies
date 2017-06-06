import sys, getopt, os

fileToFix = 'blah.off'

try:
	opts, args = getopt.getopt(sys.argv[1:], '', ["input="])
except getopt.GetoptError:
	print 'ERROR BAD INPUT'
	sys.exit(-2)
for opt, arg in opts:
	if opt == '--input':
		fileToFix = arg

fixthis = open(fileToFix, 'r')
lines = fixthis.readlines()

fixthis.close()

fixthis = open(fileToFix, 'w')

for line in lines:
	if line.isspace() == False:
		fixthis.write(line)

fixthis.close()
