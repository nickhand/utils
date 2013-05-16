
2013-02-05
----------

Contact: Rui Pereira <rui.pereira@in2p3.fr>

This distribution contains 32 spectra, in either ascii table (*.dat) or
1D FITS (*.fit) formats. All files can be checked for consistency using:

   sha1sum -c checksum.sha1

The ascii tables contain 3 space delimited columns, respectively
wavelength, flux and variance (photon noise). The 1D FITS files contain
the flux in the primary extension and the variance as a separate
extension.  The spectra have not been redshift-corrected, but have been
corrected for Milky Way extinction as described in the paper.  The
FLUXERR keyword represents the (achromatic) absolute flux calibration
accuracy [mag] as estimated by the calibration pipeline, and should be
added to the error of any synthetic magnitudes (not colors) computed
using the spectra. The filename notation is 11fe[PM]xxx.*, where [PM]
stands for Plus or Minus and xxx is the phase of observation, in days
relative to B-band maximum, times 10, eg.:

   11feM083.fit = - 8.3 d
   11feP217.fit = +21.7 d

The data contained in this distribution complement a publication from
the Nearby Supernova Factory (SNfactory).  The acquisition, reduction,
and analysis of these data are fully described in that work.  SNfactory
is providing these data to be used by other scientists for the purpose
of furthering their analyses.  Work benefitting from the use of these
data should include the following citation, conveniently included here
as an aastex bibliography entry that may be cut and pasted into a
manuscript:

\bibitem[Pereira et al.(2013)]{Pereira13} Pereira, R., et al.\ 2013,
\aap, XXX, XXX

Alternatively, the following bibtex entry may be used:

@article{2013A&A...XXX..XXXX,
}

*** PLEASE BE ADVISED *** PLEASE BE ADVISED *** PLEASE BE ADVISED ***

Some or all of the data included here will eventually be superceded by
improved reductions in a subsequent SNfactory data release.  For that
reason, we suggest that users occasionally confirm that they have the
most recent version of the data by checking the SNfactory website.  We
would also recommend that people looking to get the data for the first
time get it from the website directly.

Please do not publicly re-host the data on the web somewhere else
without checking with us first.  We would like it to be clear to
everyone the authoritative reductions are always available at the
SNfactory project web site. 
