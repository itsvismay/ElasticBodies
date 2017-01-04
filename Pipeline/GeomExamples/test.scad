///////////////////////////////////////////////////////
// PARAMETERS
thkInPlane = 2.0; //mm Thickness of the circular section before scaling
widthInPlane = 40; //width of spring in mm
thkOutOfPlane = 10; //mm

$fn=100;
diaExt = 30;
scaleX1 = 2.0;
scaleY1 = 1.0;
scaleX2 = 1.85;
scaleY2 = 1.0;
overlap = 0.5;

diaInt = diaExt - 2*thkInPlane;

//To calculate the length of the bottom flat piece (cube), we need to use "scaleX" information:
lengthOfCube = widthInPlane/2 - scaleX1*diaExt/4; //echo(lengthOfCube);

module scaled_half_circle(thkOut, dia1, dia2, scale_x1, scale_y1, scale_x2, scale_y2) {
translate([0,dia1/2,0])
intersection(){
difference() {
	scale([scale_x1, scale_y1, 1]) cylinder(h=thkOut, d1=dia1, d2=dia1, center=false);
	scale([scale_x2, scale_y2, 1]) translate([0,0,-2.5]) cylinder(h=thkOut+5, d1=dia2, d2=dia2, center=false);
	}
	translate([dia1/2,0,thkOut/2]) cube([dia1,dia1,thkOut*2], center=true);
	}
}


// 1.SECTION
cube([lengthOfCube,2,thkOutOfPlane], center=false);
translate([lengthOfCube,0,0]) scaled_half_circle(thkOutOfPlane, diaExt, diaInt, scaleX1, scaleY1, scaleX2, scaleY2);
translate([0,diaExt-2,0]) cube([lengthOfCube,2,thkOutOfPlane], center=false);

// 2.SECTION
translate([0,1*(diaExt-overlap),0]){
cube([lengthOfCube,2,thkOutOfPlane], center=false);
translate([lengthOfCube,0,0]) scaled_half_circle(thkOutOfPlane, diaExt, diaInt, scaleX1, scaleY1, scaleX2, scaleY2);
translate([0,diaExt-2,0]) cube([lengthOfCube,2,thkOutOfPlane], center=false);
}

// 3.SECTION
translate([0,2*(diaExt-overlap),0]){
cube([lengthOfCube,2,thkOutOfPlane], center=false);
translate([lengthOfCube,0,0]) scaled_half_circle(thkOutOfPlane, diaExt, diaInt, scaleX1, scaleY1, scaleX2, scaleY2);
translate([0,diaExt-2,0]) cube([lengthOfCube,2,thkOutOfPlane], center=false);
}

// MIRROR

mirror([1,0,0]){
// 1.SECTION
cube([lengthOfCube,2,thkOutOfPlane], center=false);
translate([lengthOfCube,0,0]) scaled_half_circle(thkOutOfPlane, diaExt, diaInt, scaleX1, scaleY1, scaleX2, scaleY2);
translate([0,diaExt-2,0]) cube([lengthOfCube,2,thkOutOfPlane], center=false);

// 2.SECTION
translate([0,1*(diaExt-overlap),0]){
cube([lengthOfCube,2,thkOutOfPlane], center=false);
translate([lengthOfCube,0,0]) scaled_half_circle(thkOutOfPlane, diaExt, diaInt, scaleX1, scaleY1, scaleX2, scaleY2);
translate([0,diaExt-2,0]) cube([lengthOfCube,2,thkOutOfPlane], center=false);
}

// 3.SECTION
translate([0,2*(diaExt-overlap),0]){
cube([lengthOfCube,2,thkOutOfPlane], center=false);
translate([lengthOfCube,0,0]) scaled_half_circle(thkOutOfPlane, diaExt, diaInt, scaleX1, scaleY1, scaleX2, scaleY2);
translate([0,diaExt-2,0]) cube([lengthOfCube,2,thkOutOfPlane], center=false);
}

}
