#!/usr/bin/env python3
from matplotlib import pyplot as plt
import numpy as np
import PyPore

from PyPore.DataTypes import *
import seaborn as sns

import sys
import os

from scipy import signal

class TraceFile( File ):
    '''
    A modified version of PyPore's File designed to analyze protein translocation events.
    
    '''
    def __init__( self, filename=None, current=None, timestep=None, **kwargs ):
        super(TraceFile,self).__init__(filename,current,timestep,**kwargs)
        self.filename=filename

    def parse( self, parser = lambda_event_parser( threshold=90 ) ):
        '''
        Applies one of the plug-n-play event parsers for event detection. The parser must have a .parse method
        which returns a tuple corresponding to the start of each event, and the ionic current in them. 
        '''
        
        self.events = [ Event( current=seg.current,
                               start=seg.start / self.second,
                               end=(seg.start+seg.duration) / self.second,
                               duration=seg.duration / self.second,
                               second=self.second,
                               file=self,
                               I_0=seg.I_0) for seg in parser.parse( self.current ) ]
        if len(self.events)>0:
            self.I_0 = self.events[0].I_0

        self.event_parser = parser
    def splitEvents(self,segmentCount=100):
        for event in self.events:
            if len(event.current) >= segmentCount*2:
                event.segmentCount=segmentCount
                #print(event.start*file.second,event.end*file.second)
                thisEventCurrent=np.array(self.current[int(event.start*self.second):int(event.end*self.second)])
                #print(event.current[0],event.current[-1])
                #print(thisEventCurrent[0],thisEventCurrent[-1])
                mask=np.zeros(thisEventCurrent.shape[0])
                i=0
                step=thisEventCurrent.shape[0]/segmentCount
                while(i<thisEventCurrent.shape[0]-1):
                    mask[int(i)]=1
                    i+=step
                tics = np.concatenate( ( np.where(mask==1)[0], [mask.shape[0]] ) )
                #print(np.count_nonzero(mask))
                #print(tics.shape[0])
                #if tics.shape[0] is 51:
                #    print (tics)
                del mask
                event.segments = [ Segment(current=current, copy=True, 
                               start=tics[i],
                               duration=current.shape[0] ) for i, current in enumerate( np.split( thisEventCurrent, tics[1:-1]) ) ]
                event.means=[]
                for segment in event.segments:
                    segment.event = event
                    segment.scale( float(self.second) )
                    event.means.append(segment.mean)
                    
                    #################################################
                    
                    #################################################

        pass
    
    def fft(self,n=-1,cutoff_freq=None):#number of coefficients, set to -1 for maximum
        if n==-1:
            for event in self.events:
                event.fourierCoeffs=np.fft.fft(event.current)
        else:
            for event in self.events:
                if n>len(event.current):
                    event.fourierCoeffs=[]
                    continue
                ##########################
                if cutoff_freq is not None and n<=len(event.current):
                    nyquist = self.second / 2.
                    b, a = signal.bessel(1, cutoff_freq / nyquist, btype='low', analog=0, output='ba' )
                    filtered_current = signal.filtfilt( b, a, np.array( event.current ).copy() )
                    cfft=np.fft.fft(filtered_current) #filtered cfft
                    np.put(cfft,range (n,len(event.current)),0) #filtered np.put command
                    event.fourierCoeffs=np.array(cfft,copy=True)
                ##########################
                else:
                    cfft=np.fft.fft(event.current) #original cfft
                    np.put(cfft,range (n,len(event.current)),0) #original np.put command
                    event.fourierCoeffs=np.array(cfft,copy=True)
                    
    def fft_real(self, n=-1, cutoff_freq=None):#number of coefficients, set to -1 for maximum
        if n==-1:
            for event in self.events:
                event.fourierCoeffs=np.fft.fft(event.current)
        else:
            for event in self.events:
                if (len(event.current) % 2 == 0 and (len(event.current)/2) + 1 >= n) or (len(event.current) % 2 != 0 and (len(event.current)+1)/2 >= n):
                    if cutoff_freq is not None:
                        nyquist = self.second / 2.
                        b, a = signal.bessel(1, cutoff_freq / nyquist, btype='low', analog=0, output='ba' )
                        filtered_current = signal.filtfilt( b, a, np.array( event.current ).copy() )
                        #event.fourierCoeffs=[] #holds real Fourier coefficients
                        #event.fourierPhasors=[] #holds phase angles
                        cfft=np.fft.rfft(filtered_current)/len(filtered_current)
                        #fft_phase = np.angle(cfft) #phase angles
                        cfft=np.real(cfft) #keeping the sign of the magnitude
                        #event_fourierCoeffs=np.array(cfft, copy=True)
                        #event_fourierPhasors=np.array(fft_phase, copy=True)    
                        #allfourierCoeffs.append(event_fourierCoeffs)
                        #allfourierPhasors.append(event_fourierPhasors) #not using currently
                        event.fourierCoeffs=np.array(cfft,copy=True)
                    else:
                        #event.fourierCoeffs=[] #holds real Fourier coefficients
                        #event.fourierPhasors=[] #holds phase angles
                        cfft=np.fft.rfft(event.current)/len(event.current)
                        #fft_phase = np.angle(cfft) #phase angles
                        cfft=np.real(cfft) #keeping the sign of the magnitude
                        #event_fourierCoeffs=np.array(cfft, copy=True)
                        #event_fourierPhasors=np.array(fft_phase, copy=True)    
                        #allfourierCoeffs.append(event_fourierCoeffs)
                        #allfourierPhasors.append(event_fourierPhasors) #not using currently
                        event.fourierCoeffs=np.array(cfft,copy=True)     
#########################################################################             



    def fft_segments_real(self, n=-1, cutoff_freq=None):#, n_segments = 10):#number of coefficients, set to -1 for maximum
        if n==-1:
            for event in self.events:
                event.fourierCoeffs=np.fft.fft(event.current)
        else:
            for event in self.events:
                event.segFCs=[]
                for segment in event.segments:
                    if (len(segment.current) % 2 == 0 and (len(segment.current)/2) + 1 >= n) or (len(segment.current) % 2 != 0 and (len(segment.current)+1)/2 >= n):
                        
                        if cutoff_freq is not None:
                            nyquist = self.second / 2.
                            b, a = signal.bessel(1, cutoff_freq / nyquist, btype='low', analog=0, output='ba' )
                            filtered_current = signal.filtfilt( b, a, np.array( segment.current ).copy() )
                            #event.fourierCoeffs=[] #holds real Fourier coefficients
                            #event.fourierPhasors=[] #holds phase angles
                            cfft=np.fft.rfft(filtered_current)/len(filtered_current)
                            #fft_phase = np.angle(cfft) #phase angles
                            cfft=np.real(cfft) #keeping the sign of the magnitude
                            #event_fourierCoeffs=np.array(cfft, copy=True)
                            #event_fourierPhasors=np.array(fft_phase, copy=True)    
                            #allfourierCoeffs.append(event_fourierCoeffs)
                            #allfourierPhasors.append(event_fourierPhasors) #not using currently
                            segment.fourierCoeffs=np.array(cfft,copy=True)
                            event.segFCs.append(segment.fourierCoeffs)
                        else:
                            #event.fourierCoeffs=[] #holds real Fourier coefficients
                            #event.fourierPhasors=[] #holds phase angles
                            cfft=np.fft.rfft(segment.current)/len(segment.current)
                            #fft_phase = np.angle(cfft) #phase angles
                            cfft=np.real(cfft) #keeping the sign of the magnitude
                            #event_fourierCoeffs=np.array(cfft, copy=True)
                            #event_fourierPhasors=np.array(fft_phase, copy=True)    
                            #allfourierCoeffs.append(event_fourierCoeffs)
                            #allfourierPhasors.append(event_fourierPhasors) #not using currently
                            segment.fourierCoeffs=np.array(cfft,copy=True)
                            event.segFCs.append(segment.fourierCoeffs)
                    else:
                        event.segFCs.append(None)

#########################################################################          

#########################################################################             

    def fft_correspond(self, top_FCs=-1, cutoff_freq=None, low_freq=1e3, high_freq=9e3):#, n_segments = 10):#number of coefficients, set to -1 for maximum
        sample_rate = self.second
        if top_FCs==-1:
            for event in self.events:
                event.fourierCoeffs=np.fft.fft(event.current)
        else:
            for event in self.events:
                event.segFCs=[]
                for segment in event.segments:
                    if (sample_rate/len(segment.current))*top_FCs <= high_freq:
                        cfft = np.fft.rfft(event.current)
                        #cfft = np.fft.rfft(seg_event_9.current)
                        N = len(segment.current)
                        n = np.arange(N)
                        T = N/sample_rate
                        freq = n/T
    
                        freq_array = np.array(freq)
                        #freq_range = range(1000, 10000, 1000)
            
                        freq_iter = freq_array[np.where((freq_array >= low_freq) & (freq_array <= high_freq))]
                        cfft_result = np.abs(cfft[np.where((freq_array >= low_freq) & (freq_array <= high_freq))]/N) #cfft_result = np.abs(cfft[np.where((freq_array >= 1000) & (freq_array <= 9000))]/N) 
                    
                        top_fcs = np.sort(cfft_result)[::-1]
                        top_fcs_freqs = np.sort((cfft_result*freq_iter))[::-1]
                        
                        event.segFCs.append(top_fcs_freqs[:top_FCs])
                    else:
                        event.segFCs.append(None)
                        
                        
    #########################################################################          

#########################################################################             

    def fft_welch_event(self, nsegs=5, low_freq=1e3, high_freq=9e3):#, n_segments = 10):#number of coefficients, set to -1 for maximum
        sample_rate = self.second
        for event in self.events:
            N = len(event.current)
            if round(N/nsegs) >= 1:
                f, Pxx_den = signal.welch(event.current, sample_rate, nperseg=round(N/nsegs))###welch
                freq_iter = f[np.where((f >= low_freq) & (f <= high_freq))]##welch
                cfft_result = Pxx_den[np.where((f >= low_freq) & (f <= high_freq))]##welch
            
                if cfft_result.size > 0:
                    welch_mean = np.mean(cfft_result)
                    welch_std = np.std(cfft_result)
                    event.welch_mean = welch_mean
                    event.welch_std = welch_std
                else:
                    event.welch_mean = None
                    event.welch_std = None
            else:
                    event.welch_mean = None
                    event.welch_std = None
                
    def fft_welch_segments(self, nsegs=5, low_freq=1e3, high_freq=9e3):#, n_segments = 10):#number of coefficients, set to -1 for maximum
        sample_rate = self.second
        for event in self.events:
            event.seg_welch_means=[]
            event.seg_welch_stds=[]
            for segment in event.segments:
                N = len(segment.current)
                if round(N/nsegs) >= 1:
                    f, Pxx_den = signal.welch(segment.current, sample_rate, nperseg=round(N/nsegs))###welch
                    freq_iter = f[np.where((f >= low_freq) & (f <= high_freq))]##welch
                    cfft_result = Pxx_den[np.where((f >= low_freq) & (f <= high_freq))]##welch
            
                    if cfft_result.size != 0:
                        welch_mean = np.mean(cfft_result)
                        welch_std = np.std(cfft_result)
                        event.seg_welch_means.append(welch_mean)
                        event.seg_welch_stds.append(welch_std)
                    else:
                        event.seg_welch_means.append(None)
                        event.seg_welch_stds.append(None)
                else:
                    event.seg_welch_means.append(None)
                    event.seg_welch_stds.append(None)
            
            
        
        
            
    def crossOver(self):#, n_segments = 10):#number of coefficients, set to -1 for maximum
        #sample_rate = self.second
        for event in self.events:
            event.crmeans=[]
            for segment in event.segments:
                #N = len(segment.current)
                #mean  = 51.5
                crossCount = np.sum((segment.current[:-1]>segment.mean) != (segment.current[1:]>segment.mean))
                event.crmeans.append(crossCount)
                
    def derivSignal(self):#, n_segments = 10):#number of coefficients, set to -1 for maximum
        #sample_rate = self.second
        for event in self.events:
            event.derivmeans=[]
            event.derivstds=[]
            for segment in event.segments:
                #N = len(segment.current)
                #mean  = 51.5
                deriv = np.diff(segment.current)
                event.derivmeans.append(np.mean(deriv))
                event.derivstds.append(np.std(deriv))
                
    def quantileSignal(self):#, n_segments = 10):#number of coefficients, set to -1 for maximum
        #sample_rate = self.second
        for event in self.events:
            event.first_quartiles=[]
            event.medians=[]
            event.third_quartiles=[]
            for segment in event.segments:
                #N = len(segment.current)
                #mean  = 51.5
                quant_array = np.quantile(segment.current, [0.25, 0.5, 0.75], interpolation="midpoint")
                event.first_quartiles.append(quant_array[0])
                event.medians.append(quant_array[1])
                event.third_quartiles.append(quant_array[2])
                
    
class FilterDerivativeParser( parser ):
    '''
    This parser will segment a file using a filter-derivative method. It will
    first apply a bessel filter at a certain cutoff to the current, then it will
    take the derivative of that, and segment when the derivative passes a
    threshold. It will then identify which segments are baseline current and which correspond to an event.
    '''

    def __init__( self, low_threshold=0.025, high_threshold=0.12, cutoff_freq=10000.,
        sampling_freq=1.e5 ,openholder=None,openlimits=[0.8,1.2]):
        self.low_threshold = low_threshold
        self.high_threshold = high_threshold
        self.cutoff_freq = cutoff_freq
        self.sampling_freq = sampling_freq
        self.openholder=openholder
        self.openlimits=openlimits

    def parse( self, current ):
        '''
        Apply the filter-derivative method to filter the ionic current.
        '''
        binCount=100
        edges=np.histogram(current,np.histogram_bin_edges(current,binCount))
        currentCpy=np.array(current,copy=True)
        currentCpy=currentCpy[((edges[1][-1]/2)<currentCpy)]
        upperEdges=np.histogram(currentCpy,np.histogram_bin_edges(currentCpy,int(binCount/2)))
        i_peak=np.argmax(upperEdges[0])
        I_0_MLE=np.mean(upperEdges[1][i_peak-1:i_peak+3])
        # Filter the current using a first order Bessel filter twice, one in
        # both directions to preserve phase
        from scipy import signal
        nyquist = self.sampling_freq / 2.
        b, a = signal.bessel( 1, self.cutoff_freq / nyquist, btype='low', analog=0, output='ba' )
        filtered_current = signal.filtfilt( b, a, np.array( current ).copy() )
        
        # Take the derivative
        deriv = np.abs( np.diff( filtered_current ) )
        
        # Find the edges of the blocks which fulfill pass the lower threshold
        blocks = np.where( deriv > self.low_threshold*I_0_MLE, 1, 0 )
        block_edges = np.abs( np.diff( blocks ) )
        tics = np.where( block_edges == 1 )[0] + 1 

        # Split points are points in the each block which pass the high
        # threshold, with a maximum of one per block 
        split_points = [0] 

        for start, end in zip( tics[:-1:2], tics[1::2] ): # For all pairs of edges for a block..
            segment = deriv[ start:end ] # Save all derivatives in that block to a segment
            segmentCurrent= filtered_current[start:end]
            
            if (np.argmax( segment ) > self.high_threshold*I_0_MLE) and (np.amax(segmentCurrent)>I_0_MLE*0.8): # If the maximum derivative in that block is above a threshold..
                split_points = np.concatenate( ( split_points, [ start, end ] ) ) # Save the edges of the segment 
                # Now you have the edges of all transitions saved, and so the states are the current between these transitions
        tics = np.concatenate( ( split_points, [ current.shape[0] ] ) )
        
        events=[ Segment( current=current[ tics[i]:tics[i+1] ], start=tics[i],duration=tics[i+1]-tics[i] ) 
                   for i in range( 0, len(tics)-1, 2 ) ]
        
        currents=[]
        currents = [event.current for event in events if event.min>I_0_MLE*self.openlimits[0] and event.max<I_0_MLE*self.openlimits[1]]
        conCurrents=np.hstack(currents)
        I_0 = np.mean(conCurrents) # find open pore current I_0 for the trace
        for event in events:
            event.I_0=I_0
        if self.openholder is not None: 
            self.openholder.append(conCurrents)
        return [event for event in events if event.min< I_0_MLE*0.30 and event.max<I_0_MLE*0.5]

    def set_params( self ):
        self.low_thresh = float( self.lowThreshInput.text() )
        self.high_thresh = float( self.highThreshInput.text() )
         
                    
#                     if (len(segment.current) % 2 == 0 and (len(segment.current)/2) + 1 >= n) or (len(segment.current) % 2 != 0 and (len(segment.current)+1)/2 >= n):
                        
#                         if cutoff_freq is not None:
#                             nyquist = self.second / 2.
#                             b, a = signal.bessel(1, cutoff_freq / nyquist, btype='low', analog=0, output='ba' )
#                             filtered_current = signal.filtfilt( b, a, np.array( segment.current ).copy() )
#                             #event.fourierCoeffs=[] #holds real Fourier coefficients
#                             #event.fourierPhasors=[] #holds phase angles
#                             cfft=np.fft.rfft(filtered_current)/len(filtered_current)
#                             #fft_phase = np.angle(cfft) #phase angles
#                             cfft=np.real(cfft) #keeping the sign of the magnitude
#                             #event_fourierCoeffs=np.array(cfft, copy=True)
#                             #event_fourierPhasors=np.array(fft_phase, copy=True)    
#                             #allfourierCoeffs.append(event_fourierCoeffs)
#                             #allfourierPhasors.append(event_fourierPhasors) #not using currently
#                             segment.fourierCoeffs=np.array(cfft,copy=True)
#                             event.segFCs.append(segment.fourierCoeffs)
#                         else:
#                             #event.fourierCoeffs=[] #holds real Fourier coefficients
#                             #event.fourierPhasors=[] #holds phase angles
#                             cfft=np.fft.rfft(segment.current)/len(segment.current)
#                             #fft_phase = np.angle(cfft) #phase angles
#                             cfft=np.real(cfft) #keeping the sign of the magnitude
#                             #event_fourierCoeffs=np.array(cfft, copy=True)
#                             #event_fourierPhasors=np.array(fft_phase, copy=True)    
#                             #allfourierCoeffs.append(event_fourierCoeffs)
#                             #allfourierPhasors.append(event_fourierPhasors) #not using currently
#                             segment.fourierCoeffs=np.array(cfft,copy=True)
#                             event.segFCs.append(segment.fourierCoeffs)
#                     else:
#                         event.segFCs.append(None)

#########################################################################          
             

            





