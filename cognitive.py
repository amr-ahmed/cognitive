#!/usr/bin/env python
# Copyright 2005,2007,2011 Free Software Foundation, Inc.
# This file is part of GNU Radio
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
#
import numpy as np
from numpy import fft
from gnuradio import gr, eng_notation, window
from gnuradio import audio
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from optparse import OptionParser
import sys
import math
import struct
import threading
import new_benchmark_tx_mod

sys.stderr.write("Warning: this may have issues on some machines+Python version combinations to seg fault due to the callback in bin_statitics.\n\n")

class ThreadClass(threading.Thread):
    def run(self):
        return

class tune(gr.feval_dd):
    """
    This class allows C++ code to callback into python.
    """
    def __init__(self, tb):
        gr.feval_dd.__init__(self)
        self.tb = tb

    def eval(self, ignore):
        """
        This method is called from gr.bin_statistics_f when it wants
        to change the center frequency.  This method tunes the front
        end to the new center frequency, and returns the new frequency
        as its result.
        """
  global new_freq
        try:
            # We use this try block so that if something goes wrong
            # from here down, at least we'll have a prayer of knowing
            # what went wrong.  Without this, you get a very
            # mysterious:
            #
            #   terminate called after throwing an instance of
            #   'Swig::DirectorMethodException' Aborted
            #
            # message on stderr.  Not exactly helpful ;)

            self.new_freq = self.tb.set_next_freq()
            return self.new_freq

        except Exception, e:
            print "tune: Exception: ", e

class parse_msg(object):
    def __init__(self, msg):
        self.center_freq = msg.arg1()
        self.vlen = int(msg.arg2())
        assert(msg.length() == self.vlen * gr.sizeof_float)

        # FIXME consider using NumPy array
        t = msg.to_string()
        self.raw_data = t
        self.data = struct.unpack('%df' % (self.vlen,), t)

class my_top_block(gr.top_block):
    global target_freq
    def __init__(self):
        gr.top_block.__init__(self)

        usage = "usage: %prog [options] min_freq max_freq"
        parser = OptionParser(option_class=eng_option, usage=usage)
        parser.add_option("-a", "--args", type="string", default="fpga=usrp1_fpga_4rx.rbf",
                          help="UHD device device address args [default=%default]")
        parser.add_option("", "--spec", type="string", default=None,
	                  help="Subdevice of UHD device where appropriate")
        parser.add_option("-A", "--antenna", type="string", default=None,
                          help="select Rx Antenna where appropriate")
        parser.add_option("-s", "--samp-rate", type="eng_float", default=10.66e6,
                          help="set sample rate [default=%default]")
        parser.add_option("-g", "--gain", type="eng_float", default=35,
                          help="set gain in dB (default is midpoint)")
        parser.add_option("", "--tune-delay", type="eng_float",
                          default=10e-3, metavar="SECS",
                          help="time to delay (in seconds) after changing frequency [default=%default]")
        parser.add_option("", "--dwell-delay", type="eng_float",
                          default=10e-3, metavar="SECS",
                          help="time to dwell (in seconds) at a given frequncy [default=%default]")
        parser.add_option("-F", "--fft-size", type="int", default=256,
                          help="specify number of FFT bins [default=%default]")
        parser.add_option("", "--real-time", action="store_true", default=False,
                          help="Attempt to enable real-time scheduling")

        (options, args) = parser.parse_args()
       # if len(args) != 2:
        #    parser.print_help()
          # sys.exit(1)

        self.min_freq = 470000000
        self.max_freq = 558000000

        if self.min_freq > self.max_freq:
            # swap them
            self.min_freq, self.max_freq = self.max_freq, self.min_freq

	self.fft_size = options.fft_size

        if not options.real_time:
            realtime = False
        else:
            # Attempt to enable realtime scheduling
            r = gr.enable_realtime_scheduling()
            if r == gr.RT_OK:
                realtime = True
            else:
                realtime = False
                print "Note: failed to enable realtime scheduling"

        # build graph
        self.u = uhd.usrp_source(device_addr=options.args,
                                 stream_args=uhd.stream_args('fc32'))

        # Set the subdevice spec
        if(options.spec):
            self.u.set_subdev_spec(options.spec, 0)

        # Set the antenna
        if(options.antenna):
            self.u.set_antenna(options.antenna, 0)

        usrp_rate = options.samp_rate
        self.u.set_samp_rate(usrp_rate)
        dev_rate = self.u.get_samp_rate()

	s2v = gr.stream_to_vector(gr.sizeof_gr_complex, self.fft_size)

        mywindow = window.blackmanharris(self.fft_size)
        fft = gr.fft_vcc(self.fft_size, True, mywindow)
        power = 0
        for tap in mywindow:
            power += tap*tap
	#3power = sum(map(lambda x: x*x, mywindow))
        c2mag = gr.complex_to_mag_squared(self.fft_size)
	#ref_scale=13490.0
        # FIXME the log10 primitive is dog slow
        #log = gr.nlog10_ff(20,self.fft_size, 
        #                    -20*math.log10(self.fft_size)-10*math.log10(power/self.fft_size)- 20*math.log10(ref_scale/2))
	#self.fft_size,
        # Set the freq_step to 75% of the actual data throughput.
        # This allows us to discard the bins on both ends of the spectrum.
        
        self.freq_step = 0.75 * usrp_rate
        self.min_center_freq = self.min_freq + self.freq_step/2
        self.nsteps = math.ceil((self.max_freq - self.min_freq) / self.freq_step)
        self.max_center_freq = self.min_center_freq + ((self.nsteps-1) * self.freq_step)

        self.next_freq = self.min_center_freq

        tune_delay  = max(0, int(round(options.tune_delay * usrp_rate / self.fft_size)))  # in fft_frames
        dwell_delay = max(1, int(round(options.dwell_delay * usrp_rate / self.fft_size))) # in fft_frames

        self.msgq = gr.msg_queue(16)
        self._tune_callback = tune(self)        # hang on to this to keep it from being GC'd
        stats = gr.bin_statistics_f(self.fft_size, self.msgq,
                                    self._tune_callback, tune_delay,
					dwell_delay)

        # FIXME leave out the log10 until we speed it up
	#self.connect(self.u, s2v, fft, c2mag, log, stats)
	self.connect(self.u, s2v, fft, c2mag, stats)

        if options.gain is None:
            # if no gain was specified, use the mid-point in dB
            g = self.u.get_gain_range()
            options.gain = float(g.start()+g.stop())/2.0

        self.set_gain(options.gain)
	print "gain =", options.gain

    def set_next_freq(self):
        target_freq = self.next_freq
        self.next_freq = self.next_freq + self.freq_step
        if self.next_freq >= self.max_center_freq:
            self.next_freq = self.min_center_freq

        if not self.set_freq(target_freq):
            print "Failed to set frequency to", target_freq
            sys.exit(1)

        return target_freq


    def set_freq(self, target_freq):
        """
        Set the center frequency we're interested in.

        @param target_freq: frequency in Hz
        @rypte: bool
        """
        r = self.u.set_center_freq(target_freq)
        if r:
            return True

        return False

    def set_gain(self, gain):
        self.u.set_gain(gain)

    def ab(self):
	return options.args


def main_loop(tb):
 t=0
 while(1):
    max_sweeps = 1
    num_sweeps = 0
    st = tb.nsteps * tb.fft_size
    j=0
    recmu1=[]
    rec3=[]
    rec4=[]
    freqlst = [] 
    binlst = []
    rectemp = []
    u = []
    freelst = []
    power=open("power.dat", "w")
  
    while 1:

        # Get the next message sent from the C++ code (bin_statistics)
        # It contains the center frequency and the mag squared of the fft

        m = parse_msg(tb.msgq.delete_head()) # This is a blocking call.
        fft_new = fft.fftshift(m.data)

        # - compute the square root for each element -
        # this results in the RF signal magnitude
        # since the fft outputs magnitude squared.
        # i need to test and make sure this is not too slow
        fft_new = np.sqrt(fft_new)
        fft_step = tb.freq_step / len(fft_new)
        myfreq = m.center_freq - (tb.freq_step / 2)
        
        for mypoint in fft_new:
	   	  freqlst.append(myfreq)
	    	  binlst.append(mypoint)
                  print "%s    %s" % (myfreq, mypoint)
                  power=open("power.dat", "a")
                  hh=str(myfreq)
		  bi= str(mypoint)
		  todo= hh + "     " + bi +'\n'
                  power.write(todo)
                  myfreq += fft_step	  

        # if we have done enough complete sweeps accross
        # the entire frequency range then stop/quit...
	if myfreq >= 557418945 :            
		break  
      
	#if m.center_freq == tb.max_center_freq:
         #   sys.exit(1)
	    #num_sweeps += 1
            #if num_sweeps == max_sweeps:
             #   sys.exit(0)
                # what is the proper way to stop the signal graph and quit?
    mu1=0
    mu1=np.mean(binlst)
    z=0
    while z <= 2815:
        if binlst[z] <= mu1:
            recmu1.append(binlst[z])   
        z=z+1
    mu2=np.mean(recmu1)
    std1=np.std(recmu1)
    j3=0
    i=0
    j2=0
    while i <= 2814:	
        if ( binlst[i+1]-binlst[i] ) >= (3*std1):            
             while j2 <= 256:
                 if i <= 2814:
                      rec3.append(binlst[i])
                      j3=j3+1
                      i=i+1
    	              j2=j2+1	            
                 else:
                      i=i+1
                      j2=j2+1
	else: 
		i=i+1	
        j2=0
    j4=0
    j5=j3	
    p=0
    n=0
    while p <= j5:
          p8=0
          oo=len(rec3)        
          u=rec3[p8:p8+256]
          p8=p8+256
          ll=u[0]
          pp=0
          while pp <= 255:
	        if ll < u[pp]:
		   ll=u[pp]               
                pp=pp+1
          rec4.append(ll)
          j4=j4+1
	  p= p+256
          n = n+1	
    print "-------------------------"
    ll=rec4[0]
    pp=0
    while pp < len(rec4):
	if ll >rec4[pp]:
		ll=rec4[pp]
        pp=pp+1
    min1=ll
    print "THRESHOLD = %s"%min1
    print "--------------------------"
    
    channel=21
    fp=68
    while channel < 32:
        flag=0
	y=fp
        z=fp+20
	while y < z:
		                
    		if binlst[y]>=min1:
                   flag=1
                   break
                y=y+1
        if flag==1:
            
            if channel==21:
	       print "channel %d  (from 470 to 478) "%(channel) + "occupied"
            elif channel==22:
	       print "channel %d  (from 478 to 486) "%(channel) + "occupied"
            elif channel==23:
	       print "channel %d  (from 486 to 494) "%(channel) + "occupied"
            elif channel==24:
	       print "channel %d  (from 494 to 502) "%(channel) + "occupied"
            elif channel==25:
	       print "channel %d  (from 502 to 510) "%(channel) + "occupied"
            elif channel==26:
	       print "channel %d  (from 510 to 518) "%(channel) + "occupied"    
            elif channel==27:
	       print "channel %d  (from 518 to 526) "%(channel) + "occupied"
            elif channel==28:
	       print "channel %d  (from 526 to 534) "%(channel) + "occupied"
            elif channel==29:
	       print "channel %d  (from 534 to 542) "%(channel) + "occupied"
            elif channel==30:
	       print "channel %d  (from 542 to 550) "%(channel) + "occupied"
            elif channel==31:
	       print "channel %d  (from 550 to 558) "%(channel) + "occupied"
        else:
            if channel==21:
	       print "channel %d  (from 470 to 478) "%(channel) + "free"
	       freelst.append(474000000)
            elif channel==22:
	       print "channel %d  (from 478 to 486) "%(channel) + "free"
	       freelst.append(482000000)
            elif channel==23:
	       print "channel %d  (from 486 to 494) "%(channel) + "free"
	       freelst.append(490000000)
            elif channel==24:
	       print "channel %d  (from 494 to 502) "%(channel) + "free"
	       freelst.append(498000000)
            elif channel==25:
	       print "channel %d  (from 502 to 510) "%(channel) + "free"
	       freelst.append(506000000)
            elif channel==26:
	       print "channel %d  (from 510 to 518) "%(channel) + "free"    
	       freelst.append(514000000)
            elif channel==27:
	       print "channel %d  (from 518 to 526) "%(channel) + "free"
	       freelst.append(522000000)
            elif channel==28:
	       print "channel %d  (from 526 to 534) "%(channel) + "free"
	       freelst.append(530000000)
            elif channel==29:
	       print "channel %d  (from 534 to 542) "%(channel) + "free"
	       freelst.append(538000000)
            elif channel==30:
	       print "channel %d  (from 542 to 550) "%(channel) + "free"
	       freelst.append(546000000)
            elif channel==31:
	       print "channel %d  (from 550 to 558) "%(channel) + "free"
	       freelst.append(554000000)
        fp=fp+256
	channel=channel+1
       
    new_benchmark_tx_mod.main('%4d',"fpga=usrp1_fpga_4rx.rbf") %freelst[t]
    del binlst[:]
    print "done"
    ts = len(freelst)
    if t >= ts:
        t=0
    else:
        t=t+1
    del freelst[:]
    
if __name__ == '__main__':
    t = ThreadClass()
    t.start()

    tb = my_top_block()
    try:
        tb.start()
        main_loop(tb)
        sys.exit(1)
    except KeyboardInterrupt:
        pass
