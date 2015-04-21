#!/usr/bin/ruby
#
# Time-stamp: <2013-02-12 16:42:44 zophos>
#
require 'mkmf'

#
# quick hack for preheaders
#
def have_header(header,*preheaders,&b)
  checking_for header do
    headers=preheaders.flatten
    headers+=[header]
    if try_cpp(cpp_include(headers), &b)
      $defs.push(format("-DHAVE_%s", header.tr("a-z./\055", "A-Z___")))
      true
    else
      false
    end
  end
end


unless(have_library('m','nan','math.h'))
   print <<-EOS
   ** configure error **  
   libm.a is not found. There seems to be something wrong with your environs.

   EOS
   exit(-1)
end

dir_config('narray',$archdir,$archdir)
unless(have_header('narray.h') &&
       have_header('narray_config.h'))
   print <<-EOS
   ** configure error **  
   Header narray.h and/or narray_config.h is not found.
   If you have these files in /narray/dir/include, try the following:

   % ruby extconf.rb --with-narray-include=/narray/dir/include

   EOS
   exit(-1)
end

$objs=['cwt_morlet.o']
create_makefile('cwt_morlet')
