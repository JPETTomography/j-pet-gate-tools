#!/usr/bin/perl
use strict;
use warnings;
use Forks::Super;
use File::Path qw(make_path remove_tree);

# A subroutine to return the index of a string in an array.
# Usage: indexArray(string_to_search_for, array)
sub indexArray
{
 1 while $_[0] ne pop;
 @_-1;
}

print "Perl script to execute GATE scripts in batch.\n";
print "Written by Geron Bindseil, July 2011.\n";

# Determine the number of locical CPUs available
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
my $info = Sys::Info->new;
my $cpu  = $info->device( CPU =>  );
my $cpucount = $cpu->count;

# Set the maximum number of processes to allow
$Forks::Super::MAX_PROC = 8;
$Forks::Super::ON_BUSY = 'block';

# Set the main macro file path and the output directory
my $MainMacro = "insert_main_macro_name_here.mac";
my $outputpath = '/Volumes/Data/Work/Gate/Output/' . substr($MainMacro,0, - 4);

make_path($outputpath);
print "\nOutput will be written to: $outputpath \n";

# Set the array of parameters for the gate script
# root file number to start on
my $startnumber = 0;

# Populate the Times
my $numruns = 54;
my $count = $startnumber;
my $TimeIncrement = 50;

my @TimeSlice;
my @TimeStart;
my @TimeStop;
my @OutputFileName;

for($count=$startnumber;$count<($numruns+$startnumber);$count++){
  push(@TimeSlice, $TimeIncrement);
  push(@TimeStart, $count*$TimeIncrement);
  push(@TimeStop, ($count+1)*$TimeIncrement);
  push(@OutputFileName, $outputpath . '/' . substr($MainMacro,0, - 4) . '_' . $count);
}

# Loop over the parameters and fork the GATE simulations
my @childs;
$count = 0;


for($count=0;$count<$numruns;$count++) {
  $childs[$count] = fork { cmd => "Gate -a TimeSlice $TimeSlice[$count] -a TimeStart $TimeStart[$count] -a TimeStop $TimeStop[$count] -a FileName $OutputFileName[$count] $MainMacro > $OutputFileName[$count].out" };
}

print "End of perl batch script.\n";
