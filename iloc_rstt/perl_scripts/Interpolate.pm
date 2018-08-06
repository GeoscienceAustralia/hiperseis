#
#
# Interpolate.pm: interpolation functions
#
# Istvan Bondar, 2009/07/29
#
# $d2y = spline($n, $x, $y)
#     Calculates interpolating coefficients for a natural spline 
# $yp = spline_int($xp, $n, $x, $y, $d2y)
#     Returns interpolated function value f(xp) by cubic spline interpolation
# $yp = bilinear_int($xp1, $xp2, $nx1, $nx2, $x1, $x2, $y)
#     Returns interpolated function value f(xp1,xp2) by bilinear interpolation
# ($klo, $khi) = bracket($xp, $n, $x) 
#     For a vector x, ordered in ascending order, returns indices klo and khi
#     such that x[klo] <= xp < x[khi] 
#
package Interpolate;

require Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(spline spline_int bilinear_int bracket);

use strict;

#
#  $d2y = spline($n, $x, $y)
#     Calculates interpolating coefficients for a natural spline 
#     Assumes that $y is ordered by increasing $x
#  Input Arguments:                                                  
#     $n   - number of points          
#     $x   - x array          
#     $y   - y array                     
#  Return:                                                           
#     $d2y - second derivatives of the natural spline interpolating function
#
sub spline {
    my ($n, $x, $y) = @_;
    my (@p, @d2y, $temp, $d, $i);    
    @d2y = @p = (); 
    $d2y[0] = $p[0] = $d2y[$n-1] = 0.;
    for ($i = 1; $i < $n - 1; $i++) {
        $d = ($x->[$i] - $x->[$i-1]) / ($x->[$i+1] - $x->[$i-1]);
        $temp = $d * $d2y[$i-1] + 2.;
        $d2y[$i] = ($d - 1.) / $temp;
        $p[$i] = ($y->[$i+1] - $y->[$i])   / ($x->[$i+1] - $x->[$i]) -
                 ($y->[$i]   - $y->[$i-1]) / ($x->[$i]   - $x->[$i-1]);
        $p[$i] = (6 * $p[$i] / ($x->[$i+1] - $x->[$i-1]) -$d * $p[$i-1]) / $temp;
    }
    for ($i = $n - 2; $i >= 0; $i--) {
        $d2y[$i] = $d2y[$i] * $d2y[$i+1] + $p[$i];
    }
    return \@d2y;
}

#
# $yp = spline_int($xp, $n, $x, $y, $d2y)
#     Returns interpolated function value f(xp) by cubic spline interpolation
#  Input Arguments:                                                  
#     $xp  - x point to be interpolated          
#     $n   - number of points          
#     $x   - x array          
#     $y   - y array                     
#     $d2y - second derivatives of the natural spline interpolating function
#  Return:                                                           
#     interpolated function value $yp = f($xp)
#  Calls:                                                            
#     bracket                                                  
#
sub spline_int {
    my ($xp, $n, $x, $y, $d2y) = @_; 
    my ($h, $g, $a, $b, $c, $d, $yp, $klo, $khi, $k);
#
#   bracket $xp
#
    ($klo, $khi) = bracket($xp, $n, $x);
#
#   interpolate $yp
#
    $h = $x->[$khi] - $x->[$klo];
    $g = $y->[$khi] - $y->[$klo];
    $a = ($x->[$khi] - $xp) / $h;
    $b = ($xp - $x->[$klo]) / $h;
    $c = ($a * $a * $a - $a) * $h * $h / 6.;
    $d = ($b * $b * $b - $b) * $h * $h / 6.;
    $yp = $a *   $y->[$klo] + $b *   $y->[$khi] + 
          $c * $d2y->[$klo] + $d * $d2y->[$khi];
    return $yp;
}   

#
# $yp = bilinear_int($xp1, $xp2, $nx1, $nx2, $x1, $x2, $y)
#     Returns interpolated function value f(xp1,xp2) by bilinear interpolation
#  Input Arguments:                                                  
#     xp1  - x1 point to be interpolated          
#     xp2  - x2 point to be interpolated          
#     nx1  - number of points in x1         
#     nx2  - number of points in x2         
#     x1   - x1 vector          
#     x2   - x2 vector          
#     y    - y matrix over x1 and x2                     
#  Return:                                                           
#     interpolated function value yp = f(xp1, xp2)
#  Calls:                                                            
#     bracket 
#
sub bilinear_int {
    my ($xp1, $xp2, $nx1, $nx2, $x1, $x2, $y) = @_;
    my ($ilo, $ihi, $jlo, $jhi, $f1, $f2, $yp);
#
#   bracket xp1 and xp2
#
    ($ilo, $ihi) = bracket($xp1, $nx1, $x1);
    ($jlo, $jhi) = bracket($xp2, $nx2, $x2);
#
#   scalers
#
    $f1 = ($xp1 - $x1->[$ilo]) / ($x1->[$ihi] - $x1->[$ilo]);
    $f2 = ($xp2 - $x2->[$jlo]) / ($x2->[$jhi] - $x2->[$jlo]);
#
#   interpolate
#
    $yp = (1. - $f1) * (1. - $f2) * $y->[$ilo][$jlo] + 
                $f1  * (1. - $f2) * $y->[$ihi][$jlo] + 
                $f1  *       $f2  * $y->[$ihi][$jhi] + 
          (1. - $f1) *       $f2  * $y->[$ilo][$jhi]; 
    return $yp;
}
               
#
# ($klo, $khi) = bracket($xp, $n, $x) 
#     For a vector x, ordered in ascending order, returns indices klo and khi
#     such that x[klo] <= xp < x[khi] 
#  Input Arguments:                                                  
#     $xp  - x point to be bracketed          
#     $n   - number of points in x         
#     $x   - x array          
#  Return:                                                           
#     $klo - lower index
#     $khi - upper index
#
sub bracket {
    my ($xp, $n, $x) = @_;
    my ($klo, $khi, $k);
    $klo = 0; 
    $khi = $n - 1;
    return ($klo, $khi) if ($n < 2);
    while (($khi - $klo) > 1) {
        $k = ($khi + $klo) >> 1;
        if ($x->[$k] > $xp) {
            $khi = $k;
        }
        else {
            $klo = $k;
        }
    }
    return ($klo, $khi);
}

1;

