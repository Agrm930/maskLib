#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 13:52:03 2021

@author: sasha
"""

#Run in terminal cd /Users/eddiemarici/Desktop/Masklib/maskLib && PYTHONPATH=/Users/eddiemarici/Desktop/Masklib python3 "DXF/3DMM_Rutgers1_Weak_Coupled_Transmon (1).py" 2>&1

import numpy as np

import maskLib.MaskLib as m
from maskLib.microwaveLib import Strip_straight
from maskLib.markerLib import MarkerSquare, MarkerCross
from maskLib.utilities import doMirrored
from dxfwrite import DXFEngine as dxf

# ===============================================================================
# Design Parameters - Change these to modify the design
# ===============================================================================
CHIP_SIZE_x = 21000          # Chip dimensions in microns (500x500)
CHIP_SIZE_y = 21000          # Chip dimensions in microns (500x500)
PAD_SEPARATION = 100     # Distance between pad centers in microns
PAD_WIDTH = 200          # Length of each pad in microns
PAD_HEIGHT = 200         # Width of each pad in microns
FRAME_LAYER = '703/0'    # Layer for chip frame boundary
METAL_LAYER = 'BASEMETAL'  # Layer for pads

# ===============================================================================
# wafer setup
# ===============================================================================

w = m.Wafer('chiplet array','DXF/',CHIP_SIZE_x,CHIP_SIZE_y,padding=2500,waferDiameter=m.waferDiameters['4in'],sawWidth=300,singleChipColumn=False,centerChip=False,frame=True,markers=False)
#set wafer properties
# w.frame: draw frame layer?
# w.solid: draw things solid?
# w.multiLayer: draw in multiple layers?
# w.singleChipColumn: only make one column of chips?

w.SetupLayers([
    ['BASEMETAL',4],
    ['DICEBORDER',5],
    ['Opt_Mark',3],
   # ['GOLDMARKERS',7],
    ['703/0', 9]
    ])

#setup junction layers
#setupJunctionLayers(w)

#initialize the wafer (remember to finalize any wafer properties like layers before initializing!)
w.init()


#do dicing border (by default located on layer 'MARKERS', so let's put it on layer 'DICEBORDER' instead)
w.DicingBorder(layer='DICEBORDER')

#do optical markers
#(note: mirrorX and mirrorY are true by default, but I've exposed them here to demonstrate how they work)
doMirrored(MarkerCross, w, (30000,30000),(500,500), 20,layer='Opt_Mark',mirrorX=True,mirrorY=True)
doMirrored(MarkerCross, w, (0,45000),(500,500), 20,layer='Opt_Mark',mirrorX=True,mirrorY=True)
doMirrored(MarkerCross, w, (45000,0),(500,500), 20,layer='Opt_Mark',mirrorX=True,mirrorY=True)

# ===============================================================================
# chip class definition
# ===============================================================================
class SimpleChip(m.Chip):
    def __init__(self, wafer, chipID, layer, defaults=None, **kwargs):
        m.Chip.__init__(self, wafer, chipID, layer, defaults={'w':200, 'r_out':10, 'r_ins':0})
        if defaults is not None:
            for d in defaults:
                self.defaults[d] = defaults[d]
        
        # Left pad
        Strip_straight(self, self.centered((-PAD_SEPARATION/2 - PAD_WIDTH/2, 0)), PAD_WIDTH, w=PAD_HEIGHT)
        # Right pad
        Strip_straight(self, self.centered((PAD_SEPARATION/2 + PAD_WIDTH/2, 0)), PAD_WIDTH, w=PAD_HEIGHT)
     


        
# ===============================================================================
# generate chips
# ===============================================================================
# Create and set the default chip
#import time
#cache_bust = str(int(time.time() * 1000) % 100000)
w.setDefaultChip(SimpleChip(w, f'3DMM2_CHIP_DEFAULT', w.defaultLayer))


# Populate remaining chips in buffer
for i in range(1, len(w.chips)):
    w.setChipBuffer(SimpleChip(w, f'3DMM2_CHIP{i}', w.defaultLayer).save(w), i)
    
# Save a sample chip DXF if we have enough chips
if len(w.chips) > 1:
    w.chips[10].save(w, drawCopyDXF=True, dicingBorder=False)
    

# Now that all chips are saved in the blocks section, write instances of the chips at the right spots on the wafer
w.populate()
w.save()