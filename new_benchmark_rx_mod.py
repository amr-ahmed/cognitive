#!/usr/bin/python
#!/usr/bin/env python
#
# Copyright 2010,2011 Free Software Foundation, Inc.
# 
# This file is part of GNU Radio
# 
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.


from gnuradio import gr, gru
from gnuradio import eng_notation
from gnuradio.eng_option import eng_option
from optparse import OptionParser
import xmlrpclib
from time import sleep
import time
global ddf,r

# From gr-digital
from gnuradio import digital

# from current dir
from receive_path import receive_path
from uhd_interface import uhd_receiver
from uhd_interface import uhd_interface

import xmlrpclib

import struct
import sys

#import os
#print os.getpid()

fname = raw_input('\n Type a file name: ')
if fname =='':
    print "\n You must specify a file name\n"
else:
    print "\nReady to go, file named: " + str(fname) + " will be saved to disk\n"
    output_file = open(fname, 'a')

xmlrpc_client_set_freq = xmlrpclib.Server('http://localhost:1222')

class my_top_block(gr.top_block):
    def __init__(self, demodulator, rx_callback, options):
        gr.top_block.__init__(self)

  

        if(options.rx_freq is not None):
            # Work-around to get the modulation's bits_per_symbol
            args = demodulator.extract_kwargs_from_options(options)
            symbol_rate = options.bitrate / demodulator(**args).bits_per_symbol()

            self.source = uhd_receiver(options.args, symbol_rate,
                                       options.samples_per_symbol,
                                       options.rx_freq, options.rx_gain,
                                       options.spec, options.antenna,
                                       options.verbose)
            options.samples_per_symbol = self.source._sps

        elif(options.from_file is not None):
            sys.stderr.write(("Reading samples from '%s'.\n\n" % (options.from_file)))
            self.source = gr.file_source(gr.sizeof_gr_complex, options.from_file)
        else:
            sys.stderr.write("No source defined, pulling samples from null source.\n\n")
            self.source = gr.null_source(gr.sizeof_gr_complex)
        
        # Set up receive path
        # do this after for any adjustments to the options that may
        # occur in the sinks (specifically the UHD sink)
        self.rxpath = receive_path(demodulator, rx_callback, options) 
        
        self.connect(self.source, self.rxpath)


# /////////////////////////////////////////////////////////////////////////////
#                                   main
# /////////////////////////////////////////////////////////////////////////////

global n_rcvd, n_right, received,ddf

def main():
    global n_rcvd, n_right, received ,my_access,ddf,r
      
    n_rcvd = 0
    n_right = 0
    received = 0
    my_access=10101
    freqlst=[474000000,482000000,490000000,498000000,506000000,514000000,522000000,530000000,538000000,546000000,554000000]
    r=0
    
    def rx_callback(ok, payload):
        ddf=1 
        global n_rcvd, n_right, received,ddf,r
	 
	if payload[0:2] is not '':
            (pktno,) = struct.unpack('!H', payload[0:2])
	    (access,) = struct.unpack('!H', payload[2:4])
	else:
	    print "Document received"
	    received = 1
        c = 0
	# !Q means: ! a network stream big-endian, and Q means a unsigned long long variable
	# Code only admits a payload length of 8 bits.
	#(content,) = struct.unpack('!Q', payload[0:8])
	
	if received == 0:
            n_rcvd += 1
            if ok:
                n_right += 1	    	
         
            print "ok = %5s  pktno = %4d  n_rcvd = %4d  n_right = %4d  access=%4d" % (
                ok, pktno, n_rcvd, n_right , access )
           			
	    if access != my_access or payload == None : 
			 print "changing freq to >>>>>>> %4d " % freqlst[r]
               		 xmlrpc_client_set_freq.set_freq(freqlst[r])
	        	 r=(r+1)%11		
    	    # Saving samples into a file	
	    output_file.write(payload[2:])	
	    #print "pktno = " + str(pktno) + " with content = " + str(content)
	    #output_file.write(str(content))
	else:
	    output_file.write('#EOF')
    demods = digital.modulation_utils.type_1_demods()
    
    # Create Options Parser:
    parser = OptionParser (option_class=eng_option, conflict_handler="resolve")
    expert_grp = parser.add_option_group("Expert")

    parser.add_option("-m", "--modulation", type="choice", choices=demods.keys(), 
                      default='psk',
                      help="Select modulation from: %s [default=%%default]"
                            % (', '.join(demods.keys()),))
    parser.add_option("","--from-file", default=None,
                      help="input file of samples to demod")

    receive_path.add_options(parser, expert_grp)
    uhd_receiver.add_options(parser)

    for mod in demods.values():
        mod.add_options(expert_grp)

    (options, args) = parser.parse_args ()

    if len(args) != 0:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if options.from_file is None:
        if options.rx_freq is None:
            sys.stderr.write("You must specify -f FREQ or --freq FREQ\n")
            parser.print_help(sys.stderr)
            sys.exit(1)

    # build the graph
    tb = my_top_block(demods[options.modulation], rx_callback, options)
    r = gr.enable_realtime_scheduling()
    if r != gr.RT_OK:
        print "Warning: Failed to enable realtime scheduling."
    while(1):
     ddf=0
     tb.start()        # start flow graph
     time.sleep(4)     # receive packets for 4sec 
     while ddf==1:               
           tb.stop()           
    	   tb.wait()
           ddf=0            
           tb.start()
           time.sleep(7)
           
     print ddf
     xmlrpc_client_set_freq.set_freq(freqlst[r])
     print"the current freq:%d"%freqlst[r]
     r= ( r+1 ) % 11 
     tb.stop()
     tb.wait()                 # wait for it to finish
                       
    if (n_rcvd > 0):
        PER = (n_rcvd - n_right) / float(n_rcvd)
        print "\n Packet Error Rate is: %.5f \n\n" % (PER)
    else:
	print "\n No packets received for determining the PER.\n"

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
