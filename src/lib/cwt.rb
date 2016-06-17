#
# Continuous wavelet transform using Morlet wavelet for NArray
#
# Time-stamp: <2016-06-17 15:17:30 zophos>
#
# Copyright (c) 2013 NISHI, Takao <zophos@ni.aist.go.jp> AIST
# All rights reserved.
#
# This is free software with ABSOLUTELY NO WARRANTY.
#
# You can redistribute it and/or modify it
# under the terms of GPLv2 or later.
#
require 'narray'
require 'cwt_morlet'

unless(defined?(Math::log2))
    Math::LOG2=Math.log(2.0)
    def Math.log2(x);Math.log(x)/Math::LOG2;end
end

class NVector
    def cwt(noctave,nvoice=1,w0=Math::PI*2.0)
        oct=Math::log2(noctave)
        unless oct.to_i.to_f==oct
            raise ArgumentError,'noctave must be equal to 2**n.'
        end

        o=_adjust_shape_for_cwt.cwt_morlet(noctave,
                                           nvoice,
                                           w0)[0...self.size,
                                               true]

        #
        # Wavelet scale and Fourier frequency trasformation factor
        #
        # Torrence, C. and G.P. Compo
        # "A Practical Guide to Wavelet Analysis",
        # Bull. An. Meteorol. Soc., 79(1), pp.61--78, 1998.
        # http://paos.colorado.edu/research/wavelets/bams_79_01_0061.pdf
        #
        # 3.h Wavelet scale and Fourier frequency
        #
        freq_fact=(w0+Math::sqrt(2+w0*w0))/(4*Math::PI)

        o.instance_variable_set(:@noctave,noctave)
        o.instance_variable_set(:@nvoice,nvoice)
        o.instance_variable_set(:@w0,w0)
        o.instance_variable_set(:@freq_trans_factor,freq_fact)

        o.instance_eval(<<_EOS_
def noctave
    @noctave
end
def nvoice
    @nvoice
end
def w0
    @w0
end
def freq_trans_factor
    @freq_trans_factor
end
def scale2freq(s)
    (2.0**(s.to_f/@nvoice+1.0))/@freq_trans_factor
end
def freq2scale(f)
    (Math.log2(f*@freq_trans_factor)-1)*@nvoice
end
_EOS_
                       )
        o
    end

    private
    def _adjust_shape_for_cwt
        sz=2**Math::log2(self.total.to_f).ceil

        buf=NVector.complex(sz)
        buf[0...self.total]=self.reshape(self.total)

        buf
    end
end
