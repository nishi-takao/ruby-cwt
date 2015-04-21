# CWT test data set generating script

## How to run
Just run cwt_test.rb.

    $ ./cwt_test.rb

## outputs
+ i.dat : source data (Sin curve, two perioda cycles/64 points)

    sample_no. value
     :

+ o.dat : CWT output values (noctave=8, nvoice=4)

    sample_00 Frequency_Division(octave*voice) real imag. abs.
     :
    (empty line)
    sample_01 Frequency_Division(octave*voice) real imag. abs.
     :

A empty line gives for each a source sampling point (noctave*nvoice=32 output lines)

## How to visualize

We recommend gnuplot to show results.

    $ gnuplot
    gnuplot> set pm3d map
    gnuplot> splot "o.dat" using 1:2:3 # real part
    gnuplot> splot "o.dat" using 1:2:4 # imaginary part
    gnuplot> splot "o.dat" using 1:2:5 # abs. values

