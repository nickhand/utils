#!/usr/bin/perl
#
# submit_prepare.pl -- Marshall Perrin <mperrin@astro.berkeley.edu>
# 
#     Package up an article for submission to ApJ or ArXiV.org (astro-ph)
#     This involves renaming things to meet the ApJ file name standards,
#     and/or re-compressing image files smaller to meet ArXiV's file
#     size limit, and then compressing everything into a tar file.
#
#     It also strips out comments from the LaTeX file, both to reduce
#     file size further and in case you left any embarassing comments in that
#     you don't want the journal editor to see. ;-)
#
# Usage:    submit_prepare.pl [-x] [-r #] [-t] <filename>
#
# Options:
#           -x   prepare for ArXiV submission. This involves recompressing
#                all figures to a smaller size in order to meet the < 1MB limit.
#           -r # What resolution to use for JPEGs during ArXiV preparation.
#                Default is 150 dpi; make it smaller to shrink files more.
#           -t   automatically tar things up into submit.tar.gz
#
# Outputs: A directory called "submit" which contains all the relevant stuff.
#       **WARNING** This script will overwrite any existing "submit"
#       subdirectory of the current directory. Any files therein will be lost!
#
# Requirements:
#       You must have jpeg2ps and gs in your path. Get gs from Fink and
#       jpeg2ps from http://www.pdflib.com/products/more/jpeg2ps.html
#       or DarwinPorts.
#
# History:
#       2005-10-07      Began. M. Perrin
#       2006-04-24      Added option handling, checks for file name
#                       extensions, and general cleanup of the code.
#
#
# TODO: Make it check for includegraphics commands too, and optionally
#       convert them to plotones, stripping out any arguments

use Getopt::Std;
# See http://aplawrence.com/Unix/perlgetopts.html
# TODO make this better. 
#        add an option to automatically latex and xdvi it for testing.
%options = ();
getopts("xr:t",\%options);

if ($#ARGV < 0) { die "Usage: submit_prepare.pl [-x] [-r #] [-t] <filename>\n"}

if (defined $options{r}) {
    $resolution = $options{r};
    print "  Using resolution = $resolution\n";
} else { $resolution = 150; }

# are we preparing for submission to a journal (default; use full-res figures)
# or to ArXiv.org (aka astro-ph, in which case we need to compress the figures);
$arxiv_mode = defined $options{x};

$inputtex = $ARGV[0];
unless ($inputtex =~ /\.tex$/) {$inputtex .= ".tex";};


print "Input file is $inputtex\n";

unless (-e $inputtex) {die "Couldn't find that input file!\n"}

mkdir "submit";
`rm submit/*`; # clean out any leftover files there

open IN,$inputtex or die "Error:$!";
open OUT,">submit/ms.tex" or die;

$figcounter = 1;
while (<IN>) {
	
	#Remove any comments in the LaTeX file before submitting
	s/^%.*\n|([^\\])%.*\n/$1/go;
	
	# Look to see if this line references a figure.
	if (/\\plotone{([^\}]*)}/) {
	    $fn = $1;
            # output a new version of the input line which links to
            # 'f1', 'f2', etc instead.
	    s/$fn/f$figcounter.eps/;
            prepare_figure($fn);
	}
        # repeat the above for includegraphics instead of plotone.
	if (/\\includegraphics(.*){([^\}]*)}/) {
            $args = $1;
	    $fn = $2;
            prepare_figure($fn);
	}
	
	print OUT;
	
}

# copy the BBL file, too (necessary for ApJ submission since they
# don't do BiBTeX directly)
$bblname = $inputtex;
$bblname =~ s/\.tex/\.bbl/;
`cp $bblname submit/ms.bbl`;

print "\nOutput written to submit/ms.tex\n";

if ($arxiv_mode) {
    $size = (split(/\s+/,`du -sk submit`))[0];
    print "      Total size of submit directory:  $size kB\n\n";
    if ($size >= 1024) {print "**** Warning! Size > 1 MB! ****\n"; }
}


if (defined $options{t}) {
    system('tar cvzf submit.tar.gz submit');
    $size = (split(/\s+/,`ls -l submit.tar.gz`))[4];
    $size = int($size/1024);
    print "       Created submit.tar.gz    $size kB\n\n";
}



######################################################################

sub prepare_figure {
    my ($fn);
    $fn = shift;
    print "figure: $fn\t=> f$figcounter.eps\n\t";
            # if $fn does not have an extension, assume it's .eps
            unless ($fn =~ /\.eps$/) { $fn .= ".eps"};
            #TODO error checking for the case where $fn is not an EPS file!

            unless (-e $fn) { die "I can't find the file $fn!"; }

            
       if ($arxiv_mode) {
            # Now compress the figure for ArXiV's draconian requirements
            `gs -r$resolution -dEPSCrop -dTextAlphaBits=4 -sDEVICE=jpeg -sOutputFile=submit/f$figcounter.jpg -dBATCH -dNOPAUSE $fn`;
           `jpeg2ps submit/f$figcounter.jpg > submit/f$figcounter.eps`;
           `rm submit/f$figcounter.jpg`;
            # and use whichever is smaller
            $oldsize = -s "$fn";
            $newsize = -s "submit/f$figcounter.eps";
            print "\tOld size:\t$oldsize\tNew size:\t$newsize\n";
 #           if ($oldsize < $newsize) {`cp $fn.eps submit/f$figcounter.eps"
        
        } else  {
            # We're just submitting to a journal. This is easy.
	    `cp $fn submit/f$figcounter.eps`;
        }

     $figcounter++;


}
