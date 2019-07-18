#!/usr/bin/perl
#--------------------------------------------------------#
#           plot time-history from PM energy file        #
#                                                        #
# read in time-history from an .pm energy file and       #
# create a plot using gnuplot. the plot is written as a  #
# jpeg image.                                            #
#                                                        #
# the mean and standard deviation computed, printed out  #
# and displayed on the plot as a solid black line.       #
#                                                        #
# the data from the .pm  file is written to a file in a  #
# simple format that can be read in by gnuplot.          #
#                                                        #
# input: file name                                       #
#                                                        #
# output: mean and standard deviation, png image file,   #
#         .dat file.                                     #
#                                                        #
# usage: eplot.pl  <file>                                #
#                                                        #
#--------------------------------------------------------#
use strict;

#-------------------#
# set pm file name  #
#-------------------#

if (1+$#ARGV != 1) {
  die "Usage: eplot <file> \n";
  }

my $pm_file = shift(@ARGV);

#-------------------#
# read pm file name #
#-------------------#
my ($n, $xaxis, $yaxis, $legend, $tdata, $vdata) = readEnergyFile ($pm_file);
my @time = @$tdata;
my @data = @$vdata;

my $emax = $data[0];

for (my $i = 1; $i < $n; $i++ ) {
  #print $time[$i] . ", " . $data[$i] . "\n";

  if ($data[$i] > $emax) {
    $emax = $data[$i];
    }
  }

my $mean; 
my $sd; 
my ($mean, $sd) = compStandardDeviation(@data);
print "\n\n";
print ">>> max=" . $emax . "\n";
print ">>> mean=" . $mean . "\n";
print ">>> sd=" . $sd . "\n";

#-------------------------#
# write data to .dat file #
#-------------------------#
my $name = (split('\.', $pm_file))[0];
#print ">>> name = " . $name . "\n";
my $dat_file = $name . ".dat";
open(DAT_FP,">$dat_file") || die "cant open file $dat_file: $!\n";

for (my $i = 0; $i < $n; $i++) {
  print DAT_FP $time[$i] . " " . $data[$i] . "\n";
  }

close(DAT_FP);

#--------------------#
# create plot        #
#--------------------#
my $bt = $time[0];
my $et = $time[$n-1];
createGnuPlot($xaxis, $yaxis, $legend, $mean, $sd, \@time, \@data);

#============================================================================#
#-------------------------  s u b r o u t i n e s ---------------------------#
#============================================================================#
#--------------------------------------------------------------#
#                                                              #
#                          readEnergyFile                      #
#                                                              #
#--------------------------------------------------------------#

sub
readEnergyFile {

  return undef unless(scalar(@_));
 
  my $pm_file = $_[0];

  open(PM_FP,"<$pm_file") || die "cant open file $pm_file: $!\n";

  #-------------------#
  # parse the pm file # 
  #-------------------#

  my $n = 0;
  my @time = ();
  my @vals= ();
  my $xaxis = "time";
  my $yaxis = "energy";
  my $legend = "total energy";

  while ($_ = <PM_FP>) {
    my $line = $_;
    my $i = index($line, "#");
    my $j = index($line, "@");
    #print ">>> line= " . $line . "\n";

    if ($j != -1) {
      $line =~ s/\s+/ /g;
      my @tok = split(' ', $line); 

      if ($tok[1] eq "interaction") {
        my @tok = split('\"', $line); 
        my $tok1 = $tok[1];
        }
      }
    elsif ($i == -1) {
      #my ($t, $val) = split(' ',$line);
      my @toks = split(' ',$line);

      if ($toks[0] eq "interaction") {
        }
      else {
        my $t = $toks[0];
        my $num_vals = scalar(@toks);
        $time[$n] = $t;
        $vals[$n] = $toks[$num_vals-1];
        $n += 1;
        }
      }
    }

  return ($n, $xaxis, $yaxis, $legend, \@time, \@vals);
  }

#--------------------------------------------------------------#
#                                                              #
#                       compStandardDeviation                  #
#                                                              #
#--------------------------------------------------------------#

sub 
compStandardDeviation {

  my(@data) = @_;
  my $n = scalar(@data);

  #-------------------#
  # check for no data #
  #-------------------#

  return undef unless($n);
   
  my $sum1 = 0;
  my $sum2 = 0;

  #---------------#
  # compute mean  #
  #---------------#

  foreach my $x (@data) {
    $sum1 += $x;
    }

  my $mean = $sum1 / $n; 

  #------------#
  # compute sd #
  #------------#

  foreach my $x (@data) {
    $sum2 += ($x - $mean)*($x - $mean);
    }

  my $std_dev = sqrt($sum2 / $n);
  return ($mean, $std_dev);
  }

#--------------------------------------------------------------#
#                                                              #
#                          createGnuPlot                       #
#                                                              #
# create a gnuplot and write it as an png image file.          #
#--------------------------------------------------------------#

sub
createGnuPlot {

  my ($xaxis, $yaxis, $legend, $mean, $sd, $time, $values) = @_;

  my $num_data = scalar(@{$time});
  my $start = 0;
  my $end = $num_data-1;
  my $conv = 1;
  print "\n---------- plot data ---------- \n";
  print ">>> num data = $num_data \n";

  #------------------------#
  # create pipe to gnuplot #
  #------------------------#

  open (GP, "|/usr/bin/gnuplot -persist") or die "no gnuplot";

  # force buffer to flush after each write
  use FileHandle;
  GP->autoflush(1);

  print GP "set term x11; \n";

  print GP "set nokey; \n";
  print GP "set xlabel \"Time\" \n";
  print GP "set ylabel \"Energy\" \n"; 
  print GP "set pointsize 0.75 \n";

  print GP "plot \"-\" using 1:2 with lines \n";

  for (my $i = $start; $i <= $end; $i++) {
    my $val = ${$values}[$i] / $conv;
    print GP "${$time}[$i]  $val \n";
    #print "${$time}[$i]  ${$values}[$i] \n";
    }

  print GP "end \n";
  close GP;


  }


