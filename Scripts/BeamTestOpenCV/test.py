import cv2


y_ref = 350
num = 0
def mouse_cap(event, x, y, flags, param):
	if event == cv2.EVENT_LBUTTONDOWN:
		cv2.circle(img, (x,y), 10, (255, 255, 255), 1)
		print("time ", num, y-y_ref)
		ix, iy = x, y

cap = cv2.VideoCapture("./beam_vid.mp4")

if not cap.isOpened():
	print("not opened")
else:
	
	while(cap.isOpened()):
		flag, frame = cap.read()
		cv2.namedWindow('video')
		cv2.setMouseCallback('video', mouse_cap)
		if flag and num>300:
			img = frame[100:, 300:1800]
			cv2.line(img, (0, y_ref), (1150, y_ref), (255, 255, 255), 2)
			cv2.imshow('video', img)
			cv2.waitKey()
		num+=1
		
cap.release()


