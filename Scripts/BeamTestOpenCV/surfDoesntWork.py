import cv2
from matplotlib import pyplot as plt

img1 = cv2.imread('tracker.png',0)          # queryImage

# Initiate SIFT detector
sift = cv2.xfeatures2d.SIFT_create()
kp1, des1 = sift.detectAndCompute(img1,None)

def surf_method(kp1, des1, img1, kp2, des2, img2):
	# FLANN parameters
	FLANN_INDEX_KDTREE = 0
	index_params = dict(algorithm = FLANN_INDEX_KDTREE, trees = 5)
	search_params = dict(checks=50)   # or pass empty dictionary

	flann = cv2.FlannBasedMatcher(index_params,search_params)

	matches = flann.knnMatch(des1,des2,k=2)

	# Need to draw only good matches, so create a mask
	matchesMask = [[0,0] for i in xrange(len(matches))]

	# ratio test as per Lowe's paper
	for i,(m,n) in enumerate(matches):
	    if m.distance < 0.7*n.distance:
		matchesMask[i]=[1,0]

	draw_params = dict(matchColor = (0,255,0),
		           singlePointColor = (255,0,0),
		           matchesMask = matchesMask,
		           flags = 0)

	img3 = cv2.drawMatchesKnn(img1,kp1,img2,kp2,matches,None,**draw_params)
	return img3

#----------------------------
cap = cv2.VideoCapture("./beam_vid.mp4")

if not cap.isOpened():
	print("not opened")
else:
	num = 0
	while(True):
		flag, frame = cap.read()
		#gray = cv2.cvtColor(frame[400:, 300:1800], cv2.COLOR_BGR2HSV)
		if flag:
			img2 = frame[400:, 300:1800]
			kp2, des2 = sift.detectAndCompute(img2, None)
			img3 = surf_method(kp1, des1, img1, kp2, des2, img2)
			cv2.imshow('video', img3)
			print(num)
			num+=1
			cv2.waitKey()
		
cap.release()


