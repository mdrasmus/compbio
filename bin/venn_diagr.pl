#!/usr/bin/perl

$step = .001;


$A1 = $ARGV[0]; #Area 1
$A2 = $ARGV[1]; #Area 2
$O = $ARGV[2];  #Overlap
$pi = atan2(1,1) * 4;

$r1 = sqrt($A1/$pi);
$r2 = sqrt($A2/$pi);

$x = $r1+$r2; #start condition, circles touch each other

if($O > 0){
    do {
	$x-=$step;
	
	$x1 = .5 * ($r1**2 - $r2**2 + $x**2) / $x;
	
	$A = (($pi*($r1**2)*acos($x1/$r1)/$pi) + ($pi*($r2**2)*acos(($x-$x1)/$r2)/$pi) - $x*$r1*sin(acos($x1/$r1)));
	
    }while((100*abs($A-$O)/$O) > .1);
    print stderr "Converged: $r1\t$r2\t$x\t$A\n";
}

print "%!PS\n";
print "/centershow { dup stringwidth pop -2 div -0.5 rmoveto show } bind def\n";
print "25 200 translate\n";
print "/Helvetica findfont 8 scalefont setfont\n";
print "1 setlinewidth\n";
print "1 0 0 setrgbcolor\n";
print "200 200 $r1 0 360 arc stroke\n";
printf("%.3f 200 moveto ($A1) centershow stroke\n", 200-$r1/2);
print "0 0 1 setrgbcolor\n";
printf("%.3f 200 $r2 0 360 arc stroke\n", 200+$x);
printf("%.3f 200 moveto ($A2) centershow stroke\n", 200+$x+$r2/2);
print "0 0 0 setrgbcolor\n";
printf("%.3f 200 moveto ($O) centershow stroke\n", 200+$x/2);
printf "showpage\n";

sub asin { atan2($_[0], sqrt(1 - $_[0] * $_[0])) }

sub acos { atan2( sqrt(1 - $_[0] * $_[0]), $_[0] ) }
