#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 13:52:03 2021

@author: sasha
Edited by Agrim, 2026 (junction stack, centered tabs, corner-based coordinates)
"""

#Run in terminal cd /Users/eddiemarici/Desktop/Masklib/maskLib && PYTHONPATH=/Users/eddiemarici/Desktop/Masklib python3 "DXF/3DMM_Rutgers1_Weak_Coupled_Transmon (1).py" 2>&1

import numpy as np

import maskLib.MaskLib as m
from maskLib.junctionLib import setupJunctionLayers, JcalcTabDims, JContact_slot, Transmon3DWithShunt
from maskLib.fluxoniumLib import smallJJ, leads_for_tmon_dosearray_custom
from maskLib.Entities import SolidPline, RoundRect
from maskLib.markerLib import MarkerSquare, MarkerCross
from maskLib.utilities import doMirrored, cornerRound
from dxfwrite import DXFEngine as dxf
from dxfwrite import const

def JunctionWithLeads(chip, pos):
    '''
    Draw the junction stack in the gap between the transmon pads, following
    TmonDoseArrayPrathu.py: a lead from each pad (thin line + contact pad +
    wedge taper) meeting a Dolan junction (small finger, bridge, big finger,
    with undercut and shift layers) centered in the gap.

    pos is the same point passed to Transmon3DWithShunt (bottom-left corner
    of the top pad). All dimensions come from the JJ_* design parameters.
    '''
    total_JJ_length = 2*JJ_FINGER_LENGTH + JJ_BRIDGE_LENGTH
    wedge_to_JJ = JJ_WEDGE_TO_WEDGE/2 - total_JJ_length/2
    lead_pad_to_contact = TMON_SEPARATION/2 - JJ_WEDGE_TO_WEDGE/2 - JJ_WEDGE_LENGTH - JJ_CONTACT_LENGTH + JJ_PAD_OVERLAP
    assert lead_pad_to_contact > 0, 'Leads do not fit: shrink JJ_WEDGE_TO_WEDGE or contact/wedge lengths, or increase TMON_SEPARATION'

    # chip.add() shifts objects with a .points attribute (like the transmon's
    # SolidPline pads) by chip.origin_offset, but not the plain
    # polylines/rectangles drawn here -- apply the offset manually so both
    # always land in the same place. With centerChip=False this is (0,0).
    x0 = pos[0] + chip.origin_offset[0]
    y0 = pos[1] + chip.origin_offset[1]

    xc = x0 + TMON_PAD_WIDTH/2  # center line of the pads

    for side, ystart in (('top', y0 + JJ_PAD_OVERLAP),
                         ('bottom', y0 - TMON_SEPARATION - JJ_PAD_OVERLAP)):
        leads_for_tmon_dosearray_custom(chip, m.Structure(chip, start=(xc, ystart), direction=0),
                                        toporbottom=side, layer='LEADS',
                                        leadLpadtocontact=lead_pad_to_contact, leadLcontacttoJJ=wedge_to_JJ,
                                        leadW=JJ_LEAD_WIDTH, contactW=JJ_CONTACT_WIDTH,
                                        contactL=JJ_CONTACT_LENGTH, wedgeL=JJ_WEDGE_LENGTH)

    # junction centered vertically in the gap (x nudged by leadW/2 to counter
    # the perpendicular lead offset smallJJ applies internally)
    smallJJ(chip, m.Structure(chip, start=(xc + JJ_LEAD_WIDTH/2, y0 - TMON_SEPARATION/2), direction=90),
            smallfingerlayer='SMALLFINGER', bigfingerlayer='BIGFINGER', bridgelayer='BRIDGE',
            Ulayer='UNDERCUT', Slayer='SHIFT',
            gap=JJ_BRIDGE_WIDTH, leadW=JJ_LEAD_WIDTH, fingerL=JJ_FINGER_LENGTH,
            bigfingerW=JJ_BIGFINGER_WIDTH, smallfingerW=JJ_SMALLFINGER_WIDTH,
            bridgeW=JJ_BRIDGE_WIDTH, bridgeL=JJ_BRIDGE_LENGTH, undercut=JJ_UNDERCUT)

# ===============================================================================
# Design Parameters - Change these to modify the design
# ===============================================================================
FIELD_PADDING = 500     # Padding inside each field (keeps pads away from field edge)
FIELD_SIZE = 500       # Field size in microns
FIELD_SAW = 0         # Spacing between fields (saw width) in microns
CHIPLET_SIZE_x = 21000 + FIELD_SAW      # Large chiplet dimensions in microns
CHIPLET_SIZE_y = 21000 + FIELD_SAW      # Large chiplet dimensions in microns

# Transmon3DWithShunt parameters (scaled to fit 500 micron fields)
TMON_PAD_WIDTH = 200       # Width of transmon pads in microns
TMON_PAD_HEIGHT = 100      # Height of transmon pads in microns
TMON_LEAD_WIDTH = 20       # Width of leads in microns
TMON_LEAD_HEIGHT = 100     # Height of leads in microns
TMON_PAD_RADIUS = 10       # Radius of pad corners in microns
TMON_SEPARATION = 130      # Separation between pads in microns
TMON_SHUNT_WIDTH = 5       # Width of shunt line in microns
TMON_SHUNT_DIST = 50       # Distance from pads to shunt in microns
TMON_SHUNT_LENGTH = TMON_SEPARATION + TMON_PAD_HEIGHT  # Total shunt length in microns

# Junction parameters (leads + Dolan junction drawn in the gap between the pads)
# Lead/contact geometry follows TmonDoseArrayPrathu.py, scaled down to fit the
# 130 micron pad separation (Prathu's wedge-to-wedge span of 180 assumed a 400 micron gap)
JJ_LEAD_WIDTH = 1            # Width of thin leads in microns
JJ_CONTACT_WIDTH = 20        # Width of contact pads on the leads in microns
JJ_CONTACT_LENGTH = 10       # Length of contact pads in microns
JJ_WEDGE_LENGTH = 10         # Length of wedge taper from contact pad to thin lead in microns
JJ_WEDGE_TO_WEDGE = 40       # Distance between the two wedge tips in microns (must fit inside TMON_SEPARATION)
JJ_PAD_OVERLAP = 20          # How far leads overlap into the big pads in microns
JJ_FINGER_LENGTH = 1.5       # Length of small/big fingers in microns
JJ_SMALLFINGER_WIDTH = 0.140 # Width of small finger in microns
JJ_BIGFINGER_WIDTH = 0.340   # Width of big finger in microns (small finger + 0.2)
JJ_BRIDGE_WIDTH = 0.840      # Width of bridge in microns (small finger + 0.7)
JJ_BRIDGE_LENGTH = 0.250     # Length of bridge in microns
JJ_UNDERCUT = 0.2            # Undercut width in microns

FRAME_LAYER = '703/0'      # Layer for chip frame boundary
METAL_LAYER = 'BASEMETAL'  # Layer for metal structures

RENDER_FULL_WAFER = True          # Set False to skip writing the full wafer DXF
EXPORT_SAMPLE_CHIPLET_DXF = True  # Set False to skip the standalone chiplet DXF
EXPORT_CORNER_CHIP_DXF = True     # Set False to skip the standalone corner-chip DXF
REUSE_IDENTICAL_CHIPS = True      # Reuse one generated block for repeated chiplets/corner chips

# ===============================================================================
# wafer setup
# ===============================================================================

w = m.Wafer('shunted probe array','DXF/',CHIPLET_SIZE_x,CHIPLET_SIZE_y,padding=1500,waferDiameter=m.waferDiameters['4in'],sawWidth=500,singleChipColumn=False, centerChip=False, frame=True, markers=False
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
    ['703/0', 9],
    ['LEADS',6],
    ['SMALLFINGER',1],
    ['BIGFINGER',8],
    ['BRIDGE',2],
    ['UNDERCUT',30],
    ['SHIFT',40]
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
        # centerChip=False: chip.add() applies no origin shift, so all entity
        # types share one corner-based coordinate system (no +chipsize/2 hacks)
        m.Chip.__init__(self, wafer, chipID, layer, defaults={'w':200, 'r_out':10, 'r_ins':0}, centerChip=False)
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
                
                # Transmon3D with shunt (scaled to fit 500 micron field)
                tpos = (cx - TMON_PAD_WIDTH/2, cy + TMON_SEPARATION/2)
                Transmon3DWithShunt(self, tpos,
                                    padw=TMON_PAD_WIDTH, padh=TMON_PAD_HEIGHT,
                                    leadw=TMON_LEAD_WIDTH, leadh=TMON_LEAD_HEIGHT,
                                    padradius=TMON_PAD_RADIUS,
                                    separation=TMON_SEPARATION, shunt=True,
                                    shunt_width=TMON_SHUNT_WIDTH, shunt_dist=TMON_SHUNT_DIST,
                                    shunt_length=TMON_SHUNT_LENGTH, shunt_side='left', flipped=True,
                                    tab=True, tab_shift_x=TMON_PAD_WIDTH/2 - TMON_LEAD_WIDTH/2,
                                    layer=wafer.lyr(METAL_LAYER))

                # Leads and Dolan junction in the pad gap
                JunctionWithLeads(self, tpos)
        
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
        # Temporarily override wafer chip dimensions so Chip gets the correct width/height.
        original_frame = wafer.frame
        original_chip_x = wafer.chipX
        original_chip_y = wafer.chipY
        wafer.frame = False
        wafer.chipX = CORNER_CHIP_SIZE + wafer.sawWidth
        wafer.chipY = CORNER_CHIP_SIZE + wafer.sawWidth
        # centerChip=False: chip.add() applies no origin shift, so all entity
        # types share one corner-based coordinate system (no +chipsize/2 hacks)
        m.Chip.__init__(self, wafer, chipID, layer, defaults={'w':200, 'r_out':10, 'r_ins':0}, centerChip=False)
        wafer.chipX = original_chip_x
        wafer.chipY = original_chip_y
        wafer.frame = original_frame
        
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
                    (cx - FIELD_SIZE/2, cy - FIELD_SIZE/2),
                    FIELD_SIZE, FIELD_SIZE,
                    layer=wafer.lyr(FRAME_LAYER)
                ))
                
                # Transmon3D with shunt (scaled to fit 500 micron field)
                tpos = (cx - TMON_PAD_WIDTH/2, cy + TMON_SEPARATION/2)
                Transmon3DWithShunt(self, tpos,
                                    padw=TMON_PAD_WIDTH, padh=TMON_PAD_HEIGHT,
                                    leadw=TMON_LEAD_WIDTH, leadh=TMON_LEAD_HEIGHT,
                                    padradius=TMON_PAD_RADIUS,
                                    separation=TMON_SEPARATION, shunt=True,
                                    shunt_width=TMON_SHUNT_WIDTH, shunt_dist=TMON_SHUNT_DIST,
                                    shunt_length=TMON_SHUNT_LENGTH, shunt_side='left', flipped=True,
                                    tab=True, tab_shift_x=TMON_PAD_WIDTH/2 - TMON_LEAD_WIDTH/2,
                                    layer=wafer.lyr(METAL_LAYER))

                # Leads and Dolan junction in the pad gap
                JunctionWithLeads(self, tpos)
        

         # Markers at four corners of the chiplet (placed once, outside field loop)
        marker_size = 100  # marker size in microns
        
        # Bottom-left
        MarkerCross(self, (FIELD_PADDING/2, FIELD_PADDING/2), (marker_size, marker_size), 5, layer='TiW_Mark')
        # Bottom-right
        MarkerCross(self, (CORNER_CHIP_SIZE - FIELD_PADDING/2, FIELD_PADDING/2), (marker_size, marker_size), 5, layer='TiW_Mark')
        # Top-left
        MarkerCross(self, (FIELD_PADDING/2, CORNER_CHIP_SIZE - FIELD_PADDING/2), (marker_size, marker_size), 5, layer='TiW_Mark')
        # Top-right
        MarkerCross(self, (CORNER_CHIP_SIZE - FIELD_PADDING/2, CORNER_CHIP_SIZE - FIELD_PADDING/2), (marker_size, marker_size), 5, layer='TiW_Mark')
        # Add frame for the actual corner-chip boundary.
        self.add(dxf.rectangle(
            (0, 0),
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
default_chiplet = SimpleChiplet(w, '3DMM2_CHIPLET_DEFAULT', w.defaultLayer)
w.setDefaultChip(default_chiplet)


# Populate remaining chiplets in buffer
if REUSE_IDENTICAL_CHIPS:
    for i in range(1, len(w.chips)):
        w.setChipBuffer(default_chiplet, i)
else:
    for i in range(1, len(w.chips)):
        w.setChipBuffer(SimpleChiplet(w, f'3DMM2_CHIPLET{i}', w.defaultLayer).save(w), i)
    
# Save a sample chip DXF if we have enough chips
if EXPORT_SAMPLE_CHIPLET_DXF and len(w.chips) > 1:
    w.chips[1].save(w, drawCopyDXF=True, dicingBorder=False)

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

if REUSE_IDENTICAL_CHIPS:
    corner_chip = CornerChip(w, 'CORNER_CHIP_TEMPLATE', w.defaultLayer)
    corner_chip.save(w, drawCopyDXF=EXPORT_CORNER_CHIP_DXF, dicingBorder=False)

    if RENDER_FULL_WAFER:
        for cx, cy in corner_positions:
            # Adjust insertion point to compensate for CORNER_CHIP_SIZE/2 shift
            adj_x = cx - CORNER_CHIP_SIZE / 2 - FIELD_SIZE
            adj_y = cy - CORNER_CHIP_SIZE / 2 - FIELD_SIZE
            insert_pt = w.chipSpace((adj_x, adj_y))
            w.drawing.add(dxf.insert(corner_chip.ID, insert=insert_pt, layer=w.lyr(corner_chip.layer)))
else:
    for idx, (cx, cy) in enumerate(corner_positions):
        corner_chip = CornerChip(w, f'CORNER_CHIP_{idx}', w.defaultLayer)
        corner_chip.save(w, drawCopyDXF=EXPORT_CORNER_CHIP_DXF and idx == 0, dicingBorder=False)
        # Adjust insertion point to compensate for CORNER_CHIP_SIZE/2 shift
        adj_x = cx - CORNER_CHIP_SIZE / 2 - FIELD_SIZE
        adj_y = cy - CORNER_CHIP_SIZE / 2 - FIELD_SIZE
        insert_pt = w.chipSpace((adj_x, adj_y))
        if RENDER_FULL_WAFER:
            w.drawing.add(dxf.insert(corner_chip.ID, insert=insert_pt, layer=w.lyr(corner_chip.layer)))
    
# Now that all chips are saved in the blocks section, optionally write the full wafer DXF
if RENDER_FULL_WAFER:
    w.populate()
    w.save()