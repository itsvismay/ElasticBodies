use <BezierScad.scad>;
bezierSpring(x3=-40.0, y3=20.0);
             
// height will be fixed and will be set before the optimization run.

$fn = 60;

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
    
    thk = 1.30; 
    dia1 = 4.30;
    dia2 = dia1-2*thk;
    height = 40.0;

    // Whole spring
    union() {
        bezierCurve();
        lowerPart();
        translate([0,height+dia1-thk/2,0]) upperPart();
        }

    module bezierCurve() {
        linear_extrude(height = 10) 
        BezLine([
                [x1,-thk/2], 
                [x2, y2], 
                [x3,y3],
                [x4,y4],
                [x5,height]
                ], width = [thk, width2, width3, width4, thk], 
                   resolution = 6, centered = true);
        }

    // LOWER FLAT PART
    module lowerPart() {
        translate([-10.0,-dia1,0.0]) 
        %cube([20, thk, 10], center=false);
        %translate([10,-dia1/2,0])
            difference() {
                cylinder(h=10, d1=dia1, d2=dia1 , center=false, $fn=40);
                translate([0,0,-1]) 
                cylinder(h=12, d1=dia2, d2=dia2 , center=false, $fn=40);
                translate([-2*dia1,-dia1,-1]) 
                cube([2*dia1,2*dia1,12], center=false);
                }
        translate([-10.0,-thk,0]) 
                cube([20, thk, 10], center=false);
                }
                
    // UPPER FLAT PART
    module upperPart() {
        translate([-10.0,-dia1,0.0]) 
        cube([20, thk, 10], center=false);
        %translate([10,-dia1/2,0])
            difference() {
                cylinder(h=10, d1=dia1, d2=dia1 , center=false, $fn=40);
                translate([0,0,-1]) 
                cylinder(h=12, d1=dia2, d2=dia2 , center=false, $fn=40);
                translate([-2*dia1,-dia1,-1]) 
                cube([2*dia1,2*dia1,12], center=false);
                }
        %translate([-10.0,-thk,0]) 
                cube([20, thk, 10], center=false);
                }
               
                }