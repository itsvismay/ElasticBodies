module drawBasicShape(x1, y1, x2, y2, width)
{
    // width: width of extrude
    length = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
    if (length == 0)
        echo(5);
    if (width!=0.0) {
        square(size = [width, length], center = false);
        translate([width/2, 0, 0]) circle(d=width);
        translate([width/2, length, 0]) circle(d=width);
    }
}

module translateAndRotate(x1, y1, x2, y2, width)
{    
    //angle= atan((y2-y1)/(x2-x1)) < 0 ? atan((y2-y1)/(x2-x1)) + 270 : atan((y2-y1)/(x2-x1));
    //if (y2 == y1 && x2 == x1)
    //    echo(0);
    if (width != 0.0) {
        
        angle = atan((y2-y1)/(x2-x1)) + 270;
        //translate([x1, y1, 0.0]) rotate([0.0, 0.0, angle]) translate([-width/2, 0.0, 0.0]) %drawBasicShape(x1, y1, x2, y2, width);
        if (x2>=x1) {
            translate([x1, y1, 0.0]) rotate([0.0, 0.0, angle]) translate([-width/2, 0.0, 0.0]) drawBasicShape(x1, y1, x2, y2, width);
        } 
        else if (x2<x1) { 
            translate([x1, y1, 0.0]) rotate([0.0, 0.0, angle+180]) translate([-width/2, 0.0, 0.0]) drawBasicShape(x1, y1, x2, y2, width);}
    }        
}

/////////////////////////////
//         POINTS
/////////////////////////////

$fn = 200;


union() {
// Read File Containing [X, Y, Width] Info
include <gcode_for_loop.scad>;
}