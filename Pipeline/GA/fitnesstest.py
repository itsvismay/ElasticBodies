experimentDir = ''

try:
	opts, args = getopt.getopt(sys.argv[1:], 'r', ["experimentDir="])
except getopt.GetoptError:
	print 'Error Bad Input'
	sys.exit(-2)
for opt, arg in opts:
	if opt == "--experimentDir":
		experimentDir = arg

fitFile = open(experimentDir+'fitness.txt', 'w')
fitFile.write(str(10))
fitFile.close()
