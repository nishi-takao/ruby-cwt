# ruby-cwt

NISHI, Takao <zophos@ni.aist.go.jp>

## What's this?

ruby-cwt is a library of continuous wavelet transform using Morlet
wavelet for Ruby/NArray.

It's ported from "Rwave: Time-Frequency analysis of 1-D signals".
http://cran.r-project.org/web/packages/Rwave/

## How to build

### Requirements

+ Ruby
+ NArray (with headers)
+ make
+ gcc

### Building and Installing

1. run extconf.rb
2. make
3. make install

eg;
     $ ruby extconf.rb
     $ make
     $ sudo make install

If extconf.rb miss to found narray.h, use --with-narray-include option.
eg;
     $ ruby extconf.rb --with-narray-include=/usr/lib/ruby/vendor_ruby/1.9.1/x86_64-linux


## License
Copyright (c) 2013 NISHI, Takao <zophos@ni.aist.go.jp> 
AIST
All rights reserved.

Original copyright:

(c) Copyright 1997 by
Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang
Princeton University
All right reserved.


This is free software with ABSOLUTELY NO WARRANTY.

You can redistribute it and/or modify it under the terms of GPLv2 or
later.

