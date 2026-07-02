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
from maskLib.junctionLib import setupJunctionLayers
from maskLib.markerLib import MarkerSquare, MarkerCross
from maskLib.utilities import doMirrored
from dxfwrite import DXFEngine as dxf

# ===============================================================================
# Design Parameters - Change these to modify the design
# ===============================================================================
FIELD_PADDING = 500     # Padding inside each field (keeps pads away from field edge)
FIELD_SIZE = 500       # Field size in microns
FIELD_SAW = 0         # Spacing between fields (saw width) in microns
CHIPLET_SIZE_x = 21000 + FIELD_SAW      # Large chiplet dimensions in microns
CHIPLET_SIZE_y = 21000 + FIELD_SAW      # Large chiplet dimensions in microns

PAD_SEPARATION = 130     # Distance between pad centers in microns
PAD_WIDTH = 180          # Length of each pad in microns
PAD_HEIGHT = 200         # Width of each pad in microns
FRAME_LAYER = '703/0'    # Layer for chip frame boundary
METAL_LAYER = 'BASEMETAL'  # Layer for pads

# ===============================================================================
# wafer setup
# ===============================================================================

w = m.Wafer('chiplet array','DXF/',CHIPLET_SIZE_x,CHIPLET_SIZE_y,padding=1500,waferDiameter=m.waferDiameters['4in'],sawWidth=500,singleChipColumn=False, centerChip=False, frame=True, markers=False
)
#set wafer properties
# w.frame: draw frame layer?
# w.solid: draw things solid?
# w.multiLayer: draw in multiple layers?
# w.singleChipColumn: only make one column of chips?

w.SetupLayers([
    ['BASEMETAL',4],
    ['DICEBORDER',5],
    ['Opt_Mark',3],
    ['TiW_Mark',7],
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
#doMirrored(MarkerCross, w, (30000,30000),(500,500), 20,layer='Opt_Mark',mirrorX=True,mirrorY=True)
doMirrored(MarkerCross, w, (0,45000),(500,500), 20,layer='Opt_Mark',mirrorX=True,mirrorY=True)
doMirrored(MarkerCross, w, (45000,0),(500,500), 20,layer='Opt_Mark',mirrorX=True,mirrorY=True)

# ===============================================================================
# chiplet class definition
# ===============================================================================
class SimpleChiplet(m.Chip):
    def __init__(self, wafer, chipID, layer, defaults=None, **kwargs):
        m.Chip.__init__(self, wafer, chipID, layer, defaults={'w':200, 'r_out':10, 'r_ins':0})
        if defaults is not None:
            for d in defaults:
                self.defaults[d] = defaults[d]
        
        # Tile fields inside the large chiplet
        step = FIELD_SIZE + FIELD_SAW
        nx = int((CHIPLET_SIZE_x - 2 * FIELD_PADDING) / step)
        ny = int((CHIPLET_SIZE_y - 2 * FIELD_PADDING) / step)
        print(f"Fields per chiplet: nx={nx}, ny={ny} ({nx*ny} total)")
        
        # Center the field grid in the large chiplet
        offset_x = (CHIPLET_SIZE_x - nx * step) / 2 + step / 2
        offset_y = (CHIPLET_SIZE_y - ny * step) / 2 + step / 2
        
        for ix in range(nx):
            for iy in range(ny):
                cx = offset_x + ix * step
                cy = offset_y + iy * step
                
                # Field frame on 703/0 layer
                self.add(dxf.rectangle(
                    (cx - FIELD_SIZE/2, cy - FIELD_SIZE/2),
                    FIELD_SIZE, FIELD_SIZE,
                    layer=wafer.lyr(FRAME_LAYER)
                ))
                
                # Left pad
                Strip_straight(self,
                    (cx - PAD_SEPARATION/2 - PAD_WIDTH, cy),
                    PAD_WIDTH, w=PAD_HEIGHT)
                # Right pad
                Strip_straight(self,
                    (cx + PAD_SEPARATION/2, cy),
                    PAD_WIDTH, w=PAD_HEIGHT)
        
        # Markers at four corners of the chiplet (placed once, outside field loop)
        marker_size = 100  # marker size in microns
        
        # Bottom-left
        MarkerCross(self, (FIELD_PADDING/2, FIELD_PADDING/2), (marker_size, marker_size), 5, layer='TiW_Mark')
        # Bottom-right
        MarkerCross(self, (CHIPLET_SIZE_x - FIELD_PADDING/2, FIELD_PADDING/2), (marker_size, marker_size), 5, layer='TiW_Mark')
        # Top-left
        MarkerCross(self, (FIELD_PADDING/2, CHIPLET_SIZE_y - FIELD_PADDING/2), (marker_size, marker_size), 5, layer='TiW_Mark')
        # Top-right
        MarkerCross(self, (CHIPLET_SIZE_x - FIELD_PADDING/2, CHIPLET_SIZE_y - FIELD_PADDING/2), (marker_size, marker_size), 5, layer='TiW_Mark')

# ===============================================================================
# Corner chip class definition (smaller chips for wafer corners)
# ===============================================================================
CORNER_CHIP_SIZE = int((CHIPLET_SIZE_x-FIELD_SAW-2*FIELD_PADDING)//2 + 2*FIELD_PADDING) # Size of corner chips in microns

class CornerChip(m.Chip):
    def __init__(self, wafer, chipID, layer, defaults=None, **kwargs):
        # Temporarily disable wafer frame to prevent parent from adding full-size frame
        original_frame = wafer.frame
        wafer.frame = False
        m.Chip.__init__(self, wafer, chipID, layer, defaults={'w':200, 'r_out':10, 'r_ins':0})
        wafer.frame = original_frame  # Restore original setting
        
        if defaults is not None:
            for d in defaults:
                self.defaults[d] = defaults[d]

        step = FIELD_SIZE + FIELD_SAW
        mx = int((CORNER_CHIP_SIZE - 2 * FIELD_PADDING) / step)
        my = int((CORNER_CHIP_SIZE - 2 * FIELD_PADDING) / step)
        print(f"Corner chip size: {CORNER_CHIP_SIZE} microns")
        print(f"Fields per corner chip: mx={mx}, my={my} ({mx*my} total)")
        
        # Center the field grid in the corner chip with step offset (same as SimpleChiplet)
        offset_x = (CORNER_CHIP_SIZE - mx * step) / 2 + step / 2
        offset_y = (CORNER_CHIP_SIZE - my * step) / 2 + step / 2
        
        for ix in range(mx):
            for iy in range(my):
                cx = offset_x + ix * step
                cy = offset_y + iy * step
                
                # Field frame on 703/0 layer
                self.add(dxf.rectangle(
                    (cx - FIELD_SIZE/2 + step/2, cy - FIELD_SIZE/2 + step/2),
                    FIELD_SIZE, FIELD_SIZE,
                    layer=wafer.lyr(FRAME_LAYER)
                ))
                
                # Left pad (with same step/2 offset as field frame)
                Strip_straight(self,
                    (cx - PAD_SEPARATION/2 - PAD_WIDTH + step/2, cy + step/2),
                    PAD_WIDTH, w=PAD_HEIGHT)
                # Right pad (with same step/2 offset as field frame)
                Strip_straight(self,
                    (cx + PAD_SEPARATION/2 + step/2, cy + step/2),
                    PAD_WIDTH, w=PAD_HEIGHT)
        
        # Override chip dimensions for corner chip
        self.width = CORNER_CHIP_SIZE
        self.height = CORNER_CHIP_SIZE
        self.center = (self.width / 2, self.height / 2)
        
        # Add frame with same offset alignment as fields
        self.add(dxf.rectangle(
            (step/2, step/2),
            CORNER_CHIP_SIZE, CORNER_CHIP_SIZE,
            layer=wafer.lyr(FRAME_LAYER)
        ))

        # Add a test structure or marker to identify it
        #MarkerCross(self, (0, 0), (500, 500), 10, layer='Opt_Mark')
     


        
# ===============================================================================
# generate chiplets
# ===============================================================================
# Create and set the default chiplet
#import time
#cache_bust = str(int(time.time() * 1000) % 100000)
w.setDefaultChip(SimpleChiplet(w, f'3DMM2_CHIPLET_DEFAULT', w.defaultLayer))


# Populate remaining chiplets in buffer
for i in range(1, len(w.chips)):
    w.setChipBuffer(SimpleChiplet(w, f'3DMM2_CHIPLET{i}', w.defaultLayer).save(w), i)
    
# Save a sample chip DXF if we have enough chips
if len(w.chips) > 1:
    w.chips[10].save(w, drawCopyDXF=True, dicingBorder=False)

# Add four corner chips manually
# Calculate corner positions (at ~70% of wafer radius, roughly at 45 degrees)
wafer_radius = m.waferDiameters['4in'] / 2
corner_distance = wafer_radius * 0.55

corner_positions = [
    (-corner_distance, -corner_distance),  # Bottom-left
    (corner_distance, -corner_distance),   # Bottom-right
    (-corner_distance, corner_distance),   # Top-left
    (corner_distance, corner_distance),    # Top-right
]

for idx, (cx, cy) in enumerate(corner_positions):
    corner_chip = CornerChip(w, f'CORNER_CHIP_{idx}', w.defaultLayer)
    corner_chip.save(w)
    # Adjust insertion point to compensate for CORNER_CHIP_SIZE/2 shift
    adj_x = cx - CORNER_CHIP_SIZE / 2 - FIELD_SIZE
    adj_y = cy - CORNER_CHIP_SIZE / 2 - FIELD_SIZE
    insert_pt = w.chipSpace((adj_x, adj_y))
    w.drawing.add(dxf.insert(corner_chip.ID, insert=insert_pt, layer=w.lyr(corner_chip.layer)))
    
# Now that all chips are saved in the blocks section, write instances of the chips at the right spots on the wafer
w.populate()
w.save()