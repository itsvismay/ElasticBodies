use <BezierScad.scad>;
//bezierSpring(x3=-40.0,y3=20.0);
torsionSpring();
maxTheta = 2 * PI;
numCtrlPoints = maxTheta / (0.25 * PI);
numKnots = 0;
width = [2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0];
height = [10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0];
displacement = [0.0,0.5,0.5,0.5,0.5,0.5,0.5,0.5];

numCtrlPoints = 8;
numKnots = numCtrlPoints + 4;

module torsionSpringBSpline() {
    
}

module bSpline(ctrlPts, width = [1], centered = false, showCtrls = true, resolution = 5)
    {
    if (showCtrls)
    {
        for (pt = ctrlPts)
        {
            % translate([pt[0], pt[1], 0]) circle(1);
        }
    }
    knot = [
        //0, 1/8, 2/8, 3/8, 5/8, 5/8, 6/8, 7/8, 8/8
        0/11, 1/11, 2/11, 3/11, 4/11, 5/11, 6/11, 7/11, 8/11, 9/11, 10/11, 11/11
    ];
    
    samples = len(knot) * resolution;
    change = 1 / samples;
    currentJ = 3;
    lastPointX = 0;
    lastPointY = 0;
    pointX = 0;
    pointY = 0;
    
    // test
    for (j = [0 : numCtrlPoints-1])
    {
        basisVal = Basis(j, 3, knot, change*16);
        pointX = pointX + (ctrlPts[j][0] * basisVal);
        pointY = pointY + (ctrlPts[j][1] * basisVal);
        echo(ctrlPts[j]);
        echo(basisVal);
        echo(pointX);
        echo(pointY);
    }
    
    echo(-5);
    echo(pointX);
    echo(pointY);
    
    pointX = 0;
    pointY = 0;
    
    for (i = [0 : samples-1])
    {
        // draw line in bspline
        currentT = change * i;
        currentI = i / resolution;
        pointX = 0;
        pointY = 0;
        
        for (j = [0 : numCtrlPoints-1])
        {
            basisVal = Basis(j, currentJ, knot, currentT);
            pointX = pointX + ctrlPts[j][0] * basisVal;
            pointY = pointY + ctrlPts[j][1] * basisVal;
        }
        
        drawBasicShape(lastPointX, lastPointY, pointX, pointY, width[currentI]);
        
        lastPointX = pointX;
        lastPointY = pointY;
    }
    
    
    function Basis (i, j, knot, t) =
        j == 0 && t < knot[i+1] && t >= knot[i]
        ?
            1
        :
            j == 0
            ?
                0
            :
                ((t - knot[i]) / (knot[i+j] - knot[i])) * Basis(i,j-1,knot,t) +
                ((knot[i+j+1] - t) / (knot[i+j+1] - knot[i+1])) * Basis(i+1,j-1,knot,t); 

}

module drawBasicShape(x1, y1, x2, y2, width)
{
    // width: width of extrude
    length = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
    if (width!=0.0) {
        square(size = [width, length], center = false);
        //translate([width/2, 0, 0]) circle(d=width);
        //translate([width/2, length, 0]) circle(d=width);
    }
}

module translateAndRotate(x1, y1, x2, y2, width)
{    
    if (width != 0.0) {
        angle = atan((y2-y1)/(x2-x1)) + 270;
        if (x2>=x1) {
            translate([x1, y1, 0.0]) rotate([0.0, 0.0, angle]) translate([-width/2, 0.0, 0.0]) scale(1.0, 1.0, 1.0001) drawBasicShape(x1, y1, x2, y2, width);
        } 
        else if (x2<x1) { 
            translate([x1, y1, 0.0]) rotate([0.0, 0.0, angle+180]) translate([-width/2, 0.0, 0.0]) scale(1.0, 1.0, 1.0001) drawBasicShape(x1, y1, x2, y2, width);}
    }      
}

module torsionSpring() {
    cpts = [ 
        [0.0 * cos(0 * (0.25 * PI)), 0.0 * sin(0 * (0.25*PI))],
        [5.0 * cos(1 * (0.25 * PI)), 5.0 * sin(1 * (0.25*PI))],
        [10.0 * cos(2 * (0.25 * PI)), 10.0 * sin(2 * (0.25*PI))],
        [15.0 * cos(3 * (0.25 * PI)), 15.0 * sin(3 * (0.25*PI))],
        [20.0 * cos(4 * (0.25 * PI)), 20.0 * sin(4 * (0.25*PI))],
        [25.0 * cos(5 * (0.25 * PI)), 25.0 * sin(5 * (0.25*PI))],
        [30.0 * cos(6 * (0.25 * PI)), 30.0 * sin(6 * (0.25*PI))],
        [35.0 * cos(7 * (0.25 * PI)), 35.0 * sin(7 * (0.25*PI))],
    ];
    
    cpts2 = [ 
        [0.0 * cos(0 * (45)), 0.0 * sin(0 * (45))],
        [5.0 * cos(1 * (45)), 5.0 * sin(1 * (45))],
        [10.0 * cos(2 * (45)), 10.0 * sin(2 * (45))],
        [15.0 * cos(3 * (45)), 15.0 * sin(3 * (45))],
        [20.0 * cos(4 * (45)), 20.0 * sin(4 * (45))],
        [25.0 * cos(5 * (45)), 25.0 * sin(5 * (45))],
        [30.0 * cos(6 * (45)), 30.0 * sin(6 * (45))],
        [35.0 * cos(7 * (45)), 35.0 * sin(7 * (45))],
    ];
    
    cpts3 = [ 
        [0.0 * cos(0 * (45)), 0.0 * sin(0 * (45))],
        [5.0 * cos(1 * (45)), 5.0 * sin(1 * (45))],
        [10.0 * cos(2 * (45)), 10.0 * sin(2 * (45))],
        [15.0 * cos(3 * (45)), 15.0 * sin(3 * (45))],
        [17.0 * cos(4 * (45)), 17.0 * sin(4 * (45))],
        [19.0 * cos(5 * (45)), 19.0 * sin(5 * (45))],
        [21.0 * cos(6 * (45)), 21.0 * sin(6 * (45))],
        [24.0 * cos(7 * (45)), 24.0 * sin(7 * (45))],
        //[25.0 * cos(8 * (45)), 25.0 * sin(8 * (45))],
        //[26.0 * cos(9 * (45)), 26.0 * sin(9 * (45))],
        //[27.0 * cos(10 * (45)), 27.0 * sin(10 * (45))],
        //[28.0 * cos(11 * (45)), 28.0 * sin(11 * (45))],
        //[29.0 * cos(12 * (45)), 29.0 * sin(12 * (45))],
        //[30.0 * cos(13 * (45)), 30.0 * sin(13 * (45))],
        //[31.0 * cos(14 * (45)), 31.0 * sin(14 * (45))],
        //[32.0 * cos(15 * (45)), 32.0 * sin(15 * (45))],
        //[33.0 * cos(16 * (45)), 33.0 * sin(16 * (45))],
        //[34.0 * cos(17 * (45)), 34.0 * sin(17 * (45))],
        //[35.0 * cos(18 * (45)), 35.0 * sin(18 * (45))],
        //[36.0 * cos(19 * (45)), 36.0 * sin(19 * (45))],
        //[37.0 * cos(20 * (45)), 37.0 * sin(20 * (45))],
        //[38.0 * cos(21 * (45)), 38.0 * sin(21 * (45))],
    ];
    
    cpts4 = [
        [24.0 * cos(7 * (45)), 24.0 * sin(7 * (45))],
        [25.0 * cos(8 * (45)), 25.0 * sin(8 * (45))],
        [26.0 * cos(9 * (45)), 26.0 * sin(9 * (45))],
        [27.0 * cos(10 * (45)), 27.0 * sin(10 * (45))],
        [28.0 * cos(11 * (45)), 28.0 * sin(11 * (45))],
        [29.0 * cos(12 * (45)), 29.0 * sin(12 * (45))],
        [30.0 * cos(13 * (45)), 30.0 * sin(13 * (45))],
        [31.0 * cos(14 * (45)), 31.0 * sin(14 * (45))]
    ];
    //radius = 0;
    //theta = 0;
    //num = 0;
    //for (i = [0 : numCtrlPoints])
    //{
    //    echo(len(cpts));
    //    echo(cpts[0][0]);
    //    // = radius * cos(theta);
    //    //cpts[1] = radius * sin(theta);
    //    theta = theta + 0.25 * PI;
    //    radius = radius + displacement[i];
    //}
    
    //BezLine(
    //    cpts3,
    //    width,
    //    resolution = 6,
    //    centered = true
    //);
    
    BezLine(
        cpts4,
        width,
        resolution = 6,
        centered = true
    );
    
    bSpline(
        cpts3,
        width = width,
        centered = true,
        showCtrls = true,
        resolution = 5
    );
    // to be implemented
}

module bezierSpring(x3,y3) {
	x1=0.0;
	x2=20.0;
	y2=10.0;
	x4=10.0;
	y4=30.0;
	x5=0.0;
	width2=1*1.30;
	width3=4*1.30;
	width4=1*1.30;

	thk=1.30;
	dia1=4.30;
	dia2=dia1-2*thk;
	height=40.0;
	union() {
		bezierCurve();
		//lowerPart();
		//translate([0,height+dia1-thk/2.0]) upperPart();
	}

	module bezierCurve() {
		linear_extrude(height=10)
		BezLine([
			[x1,-thk/2],
			[x2, y2],
			[x3,y3],
			[x4,y4],
			[x5,height]],
			width=[thk,width2,width3,width4,thk],resolution=6,centered=true);
	}

	//module lowerPart() {
	//	translate([-10.0,-dia1,0.0])
	//	%cube([20, thk, 10], center=false);
	//	%translate([10,-dia1/2.0])
	//		difference() {
	//			cylinder(h=10,d1=dia1,d2=dia1,center=false,$fn=40);
	//			translate([0,0,-1])
	//			cylinder(h=12, d1=dia2, d2=dia2,center=false,$fn=40);
	//			translate([-2*dia1,-dia1,-1])
	//			cube([2*dia1,2*dia1,12], center=false);			}
	//	translate([-10.0,-thk,0])
	//	cube([20,thk,10],center=false);
	//}

	//module upperPart() {
	//	translate([-10.0,-dia1,0.0])
	//	cube([20, thk, 10], center=false);
	//	%translate([10,-dia1/2,0])
	//		difference() {
	//			cylinder(h=10, d1=dia1, d2=dia1 , center=false, $fn=40);
	//			translate([0,0,-1])
	//			cylinder(h=12, d1=dia2, d2=dia2 , center=false, $fn=40);
	//			translate([-2*dia1,-dia1,-1])
	//			cube([2*dia1,2*dia1,12], center=false);
	//		}
	//	%translate([-10.0,-thk,0])
	//	cube([20, thk, 10], center=false);
	//}
}
