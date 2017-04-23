thkInPlane = 2.95391200;
widthInPlane = 11.63328000;
heightInPlane = 33.25864000;
numberOfSectionsSet = 7.91986000;
thkOutOfPlane = 9.75150000;
numberOfSections = floor(numberOfSectionsSet);
radius = 0.5 * (heightInPlane/numberOfSections + thkInPlane/2);
$fn=75;

spring();

module halfCylinder(rad, thk)
{
	difference() {
		difference() {
			scale([1,1,1]) cylinder(r=radius, h=thk, center=true);
			scale([0.9,1,1]) cylinder(r=radius-thkInPlane, h=thk+2, center=true);
		};
	translate([-rad,0,0]) cube([2*rad,2*rad,thk+2], center=true);
	};
}

module spring() {
cube([widthInPlane-radius, thkInPlane, thkOutOfPlane], center=false);
for (i=[1:numberOfSections]) {
	echo(i);
	if (i%2 != 0) {
		echo("in the first if-loop");
		translate([widthInPlane-radius,2*(i-1)*radius+radius-(i-1)*thkInPlane,thkOutOfPlane/2])
		halfCylinder(radius, thkOutOfPlane);
		translate([radius,2*i*radius-i*thkInPlane,0])
		cube([widthInPlane-2*radius,thkInPlane,thkOutOfPlane], center=false);
	}
	else {
		echo("in the else-loop");
		translate([radius,(2*i-1)*radius-(i-1)*thkInPlane,thkOutOfPlane/2])
		rotate([0,0,180])
		halfCylinder(radius, thkOutOfPlane);
		translate([radius,2*i*radius-i*thkInPlane,0])
		cube([widthInPlane-2*radius,thkInPlane,thkOutOfPlane], center=false);
	}
}
if (numberOfSections%2 != 0) {
translate([0,numberOfSections*(2*radius-thkInPlane),0])
cube([radius, thkInPlane, thkOutOfPlane], center=false);
}
else {
translate([widthInPlane-radius,numberOfSections*(2*radius-thkInPlane),0])
	cube([radius, thkInPlane, thkOutOfPlane], center=false);
}
}

