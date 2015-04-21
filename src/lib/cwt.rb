#
# Continuous wavelet transform using Morlet wavelet for NArray
#
# Time-stamp: <2013-02-18 10:19:18 zophos>
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

        _adjust_shape_for_cwt.cwt_morlet(noctave,nvoice,w0)[0...self.size,true]
    end

    private
    def _adjust_shape_for_cwt
        sz=2**Math::log2(self.total.to_f).ceil

        buf=NVector.complex(sz)
        buf[0...self.total]=self.reshape(self.total)

        buf
    end
end
