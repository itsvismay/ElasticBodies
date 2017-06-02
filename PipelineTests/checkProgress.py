import sys, os, subprocess

#name = "Individual_"
name = "TwoVarCurveFinalSectionsShort_"
ind = 0
done = 0
almostDone15 = 0
almostDone14 = 0
almostDone13 = 0
almostDone12 = 0
almostDone11 = 0
almostDone10 = 0
almostDone9 = 0
almostDone8 = 0
almostDone7 = 0
almostDone6 = 0
almostDone5 = 0
almostDone4 = 0
almostDone3 = 0
almostDone2 = 0
almostDone1 = 0
almostDone0 = 0
almostDonehalf = 0
broke = 0

completedInd = []
brokeInd = []

for ind in range(0, 500):
  try:
    f = open(name+str(ind)+'/curvePrepedRemesh/position.txt')
    lines = f.readlines()
    print ind, len(lines)
    if len(lines) == 5002:
      done = done + 1
    elif len(lines) >= 1500:
      almostDone15 = almostDone15 + 1
    elif len(lines) >= 1400:
      almostDone14 = almostDone14 + 1
    elif len(lines) >= 1300:
      almostDone13 = almostDone13 + 1
    elif len(lines) >= 1200:
      almostDone12 = almostDone12 + 1
    elif len(lines) >= 1100:
      almostDone11 = almostDone11 + 1
    elif len(lines) >= 1000:
      almostDone10 = almostDone10 + 1
    elif len(lines) >= 900:
      almostDone9 = almostDone9 + 1
    elif len(lines) >= 800:
      almostDone8 = almostDone8 + 1
    elif len(lines) >= 700:
      almostDone7 = almostDone7 + 1
    elif len(lines) >= 600:
      almostDone6 = almostDone6 + 1
    elif len(lines) >= 500:
      almostDone5 = almostDone5 + 1
    elif len(lines) >= 400:
      almostDone4 = almostDone4 + 1
    elif len(lines) >= 300:
      almostDone3 = almostDone3 + 1
    elif len(lines) >= 200:
      almostDone2 = almostDone2 + 1
    elif len(lines) >= 100:
      almostDone1 = almostDone1 + 1
    elif len(lines) >= 80:
      almostDonehalf = almostDonehalf + 1
    else:
      almostDone0 = almostDone0 + 1
  except:
    blah = 0
    broke = broke + 1
    brokeInd.append(ind)

print brokeInd

print "Finished:", done
print "1500: ", almostDone15
print "1400: ", almostDone14
print "1300: ", almostDone13
print "1200: ", almostDone12
print "1100: ", almostDone11
print "1000: ", almostDone10
print "900: ", almostDone9
print "800: ", almostDone8
print "700: ", almostDone7
print "600: ", almostDone6
print "500: ", almostDone5
print "400: ", almostDone4
print "300: ", almostDone3
print "200: ", almostDone2
print "100: ", almostDone1
print "50: ", almostDonehalf
print "000: ", almostDone0
print "Broke:", broke
