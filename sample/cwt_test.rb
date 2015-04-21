#!/usr/bin/ruby
#
# Testing data generator for CWT
#
# % ./cwt_test.rb
#
# Time-stamp: <2015-04-21 18:27:56 zophos>
#
require 'narray'
require 'cwt'

SRCDATA_FILE='i.dat'
DSTDATA_FILE='o.dat'

########################################################################
#
# create source data
#
# 2 cycles Sin curve par 64 samples
#
# v[x]=sin(x/64*4*PI)
#
v=NVector.float(64).indgen!
v.mul!(NMath::PI/16.0)
v=NMath.sin(v)


########################################################################
#
# continuous wavelet transform
# octave = 8
# voice par octave = 4
#
# => result[64,32] (complex)
#
w=v.cwt(8,4)

w_real=w.real # extract real part
w_imag=w.imag # extract imaginary part
w_abs=w.abs   # absolute value of complex number


(xr,yr)=w.shape


########################################################################
#
# dump source data
#
# data format:
#
# sample#00 source_value\n
#   :
# \n
# sample#01 source_value\n
#   :
#
File.open(SRCDATA_FILE,'w'){|f|
    xr.times{|x|
        f<<"#{x} #{v[x]}\n" #
    }
}

########################################################################
#
# dump result
#
# data format:
#
# sample# octave*voice real imagenary absolute\n
#   :
#
File.open(DSTDATA_FILE,'w'){|f|
    xr.times{|x|
        yr.times{|y|
            f<<"#{x} #{y} #{w_real[x,y]} #{w_imag[x,y]} #{w_abs[x,y]}\n"
        }
        f<<"\n"
    }
}
