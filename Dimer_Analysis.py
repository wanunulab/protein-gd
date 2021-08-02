#!/usr/bin/env python3
from matplotlib import pyplot as plt
import numpy as np
import PyPore

from PyPore.DataTypes import *
import seaborn as sns

import sys
import os
# def prepareFiles(directory,metaOutput="metadata.json",verbose=False ,acceptedExts={".bin",".opt",".abf"},meta_format="json"):
#     for root,dirs,files in os.walk(directory,topdown=True):
#         print (root)
        
#         for name in files:
#             if os.path.splitext(name)[1]==".bin":
#                 print (" ",name)
        
            

# prepareFiles("F:\\ProteinData\\Meni\\")


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

        pass
    def fft(self,n=-1):#number of coefficients, set to -1 for maximum
        if n==-1:
            for event in self.events:
                event.fourierCoeffs=np.fft.fft(event.current)
        else:
            for event in self.events:
                if n>len(event.current):
                    event.fourierCoeffs=[]
                    continue
                cfft=np.fft.fft(event.current)
                np.put(cfft,range (n,len(event.current)),0)
                event.fourierCoeffs=np.array(cfft,copy=True)
            
                
class Experiment( object ):
    '''
    An experiment represents a series of files which all are to be analyzed together, and have
    their results bundled. This may be many files run on the same nanopore experiment, or the
    same conditions being tested over multiple days. It attempts to construct a fast and efficient
    way to store the data. Functions are as close to the file functions as possible, as you are
    simply applying these functions over multiple files.
    '''

    def __init__( self, filenames, name=None ):
        '''
        Take in the filenames and store them for later analysis
        '''

        self.filenames = filenames
        self.name = name or "Experiment"
        self.files = []

    def parse( self, event_detector=lambda_event_parser( threshold=90 ), 
        segmenter=SpeedyStatSplit( prior_segments_per_second=10, cutoff_freq=2000. ),
        filter_params=(1,2000),
        verbose=True, meta=False  ):
        '''
        Go through each of the files and parse them appropriately. If the segmenter
        is set to None, then do not segment the events. If you want to filter the
        events, pass in filter params of (order, cutoff), otherwise None.
        '''

        # Go through each file one at a time as a generator to ensure many files
        # are not open at the same time.
        for file in map( TraceFile, self.filenames ):
            if verbose:
                print ("Opening {}".format( file.filename ))

            file.parse( parser=event_detector )

            if verbose:
                print ("\tDetected {} Events".format( file.n ))
            
            # If using a segmenter, then segment all of the events in this file
            for i, event in enumerate( file.events ):
                if filter_params is not None:
                    event.filter( *filter_params )
                if segmenter is not None:
                    event.parse( parser=segmenter )
                    if verbose:
                        print ("\t\tEvent {} has {} segments".format( i+1, event.n ))

            if meta:
                file.to_meta()
            self.files.append( file )

    def apply_hmm( self, hmm, filter=None, indices=None ):
        segments = []
        for event in self.get( "events", filter=filter, indices=indices ):
            _, segs = hmm.viterbi( np.array( seg.mean for seg in event.segments ) )
            segments = np.concatenate( ( segments, segs ) )
        return segments

    def delete( self ):
        with ignored( AttributeError ):
            del self.events

        with ignored( AttributeError ):
            del self.segments

        for file in self.files:
            file.delete()
        del self

    @property
    def n( self ):
        return len( self.files )

    @property
    def events( self ):
        '''
        Return all the events in all files.
        '''
        
        try:
            return reduce( list.__add__, [ file.events for file in self.files ] )
        except:
            return []

    @property
    def segments( self ):
        '''
        Return all segments from all events in an unordered manner.
        '''

        try:
            return reduce( list.__add__, [ event.segments for event in self.events ] )
        except:
            return []
 
class Sample( object ):
    '''A container for events all suggested to be from the same substrate.'''
    def __init__( self, events=[], files=[], label=None ):
        self.events = events
        self.files = files
        self.label = label

    def delete( self ):
        with ignored( AttributeError ):
            for file in self.files:
                file.delete()

        for event in self.events:
            event.delete()
        del self.events
        del self.files
        del self

# =============================================================================
#     self.events = [ Event( current=seg.current,
#                                start=seg.start / self.second,
#                                end=(seg.start+seg.duration) / self.second,
#                                duration=seg.duration / self.second,
#                                second=self.second,
#                                file=self ) for seg in parser.parse( self.current ) ]
#     
# =============================================================================

# def _lambda_select( self, events ):
#         '''
#         From all of the events, filter based on whatever set of rules has been initiated with.
#         ''' 
#         return [ event for event in events if np.all( [ rule( event ) for rule in self.rules ] ) ]
    
#     def parse( self, current ):
#         '''
#         Perform a large capture of events by creating a boolean mask for when the current is below a threshold,
#         then detecting the edges in those masks, and using the edges to partitition the sample. The events are
#         then filtered before being returned. 
#         '''
#         mask = np.where( current < self.threshold, 1, 0 ) # Find where the current is below a threshold, replace with 1's
#         mask = np.abs( np.diff( mask ) )                  # Find the edges, marking them with a 1, by derivative
#         tics = np.concatenate( ( [0], np.where(mask ==1)[0]+1, [current.shape[0]] ) )
#         del mask
#         events = [ Segment(current=np.array(current), copy=True, 
#                             start=tics[i],
#                             duration=current.shape[0] ) for i, current in enumerate( np.split( current, tics[1:-1]) ) ]
#         return [ event for event in self._lambda_select( events ) ]


# #mean_list = ' '.join( str( round( segment.mean, 8 ) ) for segment in Event.segments )









