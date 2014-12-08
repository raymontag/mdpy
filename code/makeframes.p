#!/usr/bin/perl
#
# makeframes.p   
#
# Eric Jeckelmann -- University of Mainz, October 2003 
#                 -- Leibniz Universitaet Hannover, July 2011
#
# make PostScript frames showing hard-core particle motion
#

$radius=1;   # hard-core particle radius

# frame size
$xbound=400;
$ybound=400;
$wallwidth=4;

$frame=0;
$file = $ARGV[0];
open(FILE,$file);
while (<FILE>) 
    {
# read in system size and draw box frame
    if(/System size = (\d+) x (\d+)/)
        {
        $nx = $1;
        $ny = $2;
        print STDERR "nx, ny =  $nx, $ny \n";
		if($ny > $nx) { $ybound *= $ny/$nx; }
		if($nx > $ny) { $xbound *= $nx/$ny; }
        $xscale = ($xbound-2*$wallwidth)/$nx;
        $yscale = ($ybound-2*$wallwidth)/$ny;
        $xcorner=$wallwidth/$xscale;
        $ycorner=$wallwidth/$yscale;
        $wallwidth=sqrt($xcorner*$ycorner);
        }
# read in number of particles
    elsif(/Number of particles = (\d+)/)
        {
        $np = $1;
        print STDERR "np =  $np \n";
        }
# read in the movie time step  
    elsif(/Movie time step = (.*)/)
        {
        $movie_time_step = $1;
        print STDERR "movie time step =  $movie_time_step \n";
        }
# read in time
    elsif(/Time = (.*)/)
        {
        print STDERR "Frame time = $1 \n";

# making PS file        
        print STDERR "    saving frame =  $frame \n";
        if($frame < 10) { $filename="frame-000".$frame; }
        elsif($frame < 100) { $filename="frame-00".$frame; }
        elsif($frame < 1000) { $filename="frame-0".$frame; }
        else { $filename="frame-".$frame; }
        open(OUT,">$filename\.ps"); 
        $frame++;
        print OUT "%!PS-Adobe-3.0 \n";
        print OUT "%%BoundingBox: 0 0 $xbound $ybound \n"; 
        print OUT "/circle { currentpoint 1 0 360 gsave newpath 0.9 0.3 0.3 setrgbcolor 0.02 setlinewidth arc fill stroke grestore } def \n"; 
        print OUT "/bluecircle { currentpoint 1 0 360 gsave newpath 0.0 0.3 1.0 setrgbcolor 0.02 setlinewidth arc fill stroke grestore } def \n"; 
        print OUT "/lightgreencircle { currentpoint 1 0 360 gsave newpath 0.3 1.0 0.3 setrgbcolor 0.02 setlinewidth arc fill stroke grestore } def \n"; 
        print OUT "/deletecircle { currentpoint 1 0 360 gsave newpath 1.0 1.0 1.0 setrgbcolor 0.02 setlinewidth arc fill stroke grestore } def \n"; 
        print OUT "$xscale $yscale scale \n";
        print OUT "$xcorner $ycorner translate \n";
        print OUT 0, " " , 0," moveto \n";
        print OUT 0, " ", $ny, " lineto \n";
        print OUT $nx, " " , $ny, " lineto \n"; 
        print OUT $nx, " " , 0, " lineto \n";
        print OUT "closepath \n";
        print OUT "1.0 1.0 1.0 setrgbcolor fill \n";
#        print OUT "0.0 0.0 0.5 setrgbcolor \n";
#        print OUT "$wallwidth setlinewidth \n";
#        print OUT "1 setlinecap \n";
#        print OUT -$xcorner/2, " " , -$ycorner/2," moveto \n";
#        print OUT -$xcorner/2, " ", $ny+$ycorner/2, " lineto \n";
#        print OUT $nx+$xcorner/2, " " , $ny+$ycorner/2, " lineto \n"; 
#        print OUT $nx+$xcorner/2, " " , -$ycorner/2, " lineto \n";
        print OUT "closepath \n";
        print OUT "stroke \n";
        }
    elsif(/Particle (\d+) : (.*) (.*) (.*) (.*)/)
        {
        if($np > 3 && $np%2 == 1 && $1==$np)
            {
            print OUT "$2 $3 moveto bluecircle\n"; 
            }
        else
            {
            print OUT "$2 $3 moveto circle\n"; 
            }
        if($1==$np) 
            { 
            print OUT "showpage \n";
            close(OUT);
            }
        }
    }
close(FILE);





