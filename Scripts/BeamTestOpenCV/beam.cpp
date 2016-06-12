// #include <iostream>
// #include <opencv2/opencv.hpp>
// #include <fstream>
 
// using namespace std;
// using namespace cv;
// ofstream dataFile;
 
// int main (int argc, const char * argv[])
// {
//     dataFile.open("data.txt");
//     VideoCapture cap("../beam1.mp4");
//     cap.set(CV_CAP_PROP_FRAME_WIDTH, 640);
//     cap.set(CV_CAP_PROP_FRAME_HEIGHT, 480);    
//     if (!cap.isOpened())
//         return -1;
 
//     Mat img;
//     HOGDescriptor hog;
//     hog.setSVMDetector(HOGDescriptor::getDefaultPeopleDetector());
    
//     namedWindow("video capture", CV_WINDOW_AUTOSIZE);
//     int t=0;
//     while (true)
//     {
//         t++;
//         cap >> img;
//         if (!img.data)
//             break;
 
//         vector<Rect> found, found_filtered;
//         hog.detectMultiScale(img, found, 0, Size(8, 8), Size(16, 16), 1.2, 2, false);
 
//         size_t i, j;
//         for (i=0; i<found.size(); i++)
//         {
//             Rect r = found[i];
//             for (j=0; j<found.size(); j++)
//                 if (j!=i && (r & found[j])==r)
//                     break;
//             if (j==found.size())
//                 found_filtered.push_back(r);
//         }
//         double avgD = 0;
//         for (i=0; i<found_filtered.size(); i++)
//         {
//     	    Rect r = found_filtered[i];
//             r.x += cvRound(r.width*0.1);
//             r.width = cvRound(r.width*0.8);
//             r.y += cvRound(r.height*0.06);
//     	    r.height = cvRound(r.height*0.9);
//             if(r.br().y>=img.rows/2){
//                 avgD+=r.br().x;  
//     	       rectangle(img, r.tl(), r.br(), cv::Scalar(0,255,0), 2);
//             }
// 	    }
//         if(found_filtered.size()>0){
//             avgD = avgD/found_filtered.size();
//         }
//         dataFile<<t<<","<<(1-(1.0/avgD))<<"\n";
//         dataFile.flush();
//         imshow("video capture", img);
//         if (waitKey(20) >= 0)
//             break;
//     }
//     dataFile.close();
//     return 0;
// }

#include <stdio.h>
#include <iostream>
#include "opencv2/core.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/highgui.hpp"

using namespace cv;
using namespace cv::xfeatures2d;

void readme();

/** @function main */
int main( int argc, char** argv )
{

  // Mat img_1 = imread( "../scsht.png", IMREAD_GRAYSCALE );
  Mat img_2 = imread( "../siftimg.png", IMREAD_GRAYSCALE );
  
  VideoCapture cap("../test.mp4");
  cap.set(CV_CAP_PROP_FRAME_WIDTH, 640);
  cap.set(CV_CAP_PROP_FRAME_HEIGHT, 480);    
  if (!cap.isOpened())
    return -1;

  Mat img_1;
  namedWindow("video", CV_WINDOW_AUTOSIZE);

    //-- Step 1: Detect the keypoints using SURF Detector
    int minHessian = 6000;
    Ptr<SURF> detector = SURF::create( minHessian );

  while(true){
    cap>>img_1;
    if( !img_1.data || !img_2.data )
    { std::cout<< " --(!) Error reading images " << std::endl; return -1; }

    std::vector<KeyPoint> keypoints_1, keypoints_2;

    detector->detect( img_1, keypoints_1 );
    detector->detect( img_2, keypoints_2 );

    // //-- Draw keypoints
    Mat img_keypoints_1; Mat img_keypoints_2;

    drawKeypoints( img_1, keypoints_1, img_keypoints_1, Scalar::all(-1), DrawMatchesFlags::DEFAULT );
    drawKeypoints( img_2, keypoints_2, img_keypoints_2, Scalar::all(-1), DrawMatchesFlags::DEFAULT );

    //-- Show detected (drawn) keypoints
    // imshow("video", img_keypoints_1 );
    // imshow("Keypoints 2", img_keypoints_2 );
    imshow("video", img_keypoints_1);
    waitKey(1);
  }
  


  return 0;
  }

  /** @function readme */
  void readme()
  { std::cout << " Usage: ./SURF_detector <img1> <img2>" << std::endl; }