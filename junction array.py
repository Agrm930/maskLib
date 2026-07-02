#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 13:52:03 2021

@author: sasha
"""

#Run in terminal cd /Users/eddiemarici/Desktop/Masklib/maskLib && PYTHONPATH=/Users/eddiemarici/Desktop/Masklib python3 "DXF/3DMM_Rutgers1_Weak_Coupled_Transmon (1).py" 2>&1

import numpy as np

import maskLib.MaskLib as m
from maskLib.junctionLib import setupJunctionLayers, JcalcTabDims, JContact_slot
from maskLib.Entities import SolidPline, RoundRect
from maskLib.markerLib import MarkerSquare, MarkerCross
from maskLib.utilities import doMirrored, cornerRound
from dxfwrite import DXFEngine as dxf
from dxfwrite import const


def Transmon3DWithShunt(chip, pos, padw=1500, padh=750, padw2=None, padh2=None, leadw=100, leadw2=None, leadh=2000, leadh2=None, separation=200, padradius=20,
                        tab=False, tab_gapw=3, tab_gapl=0.5, tab_tabw=2, tab_tabl=0.5, tab_taboffs=-0.5, tab_r_out=1.5, tab_r_ins=1.5,
                        tab_offset_x=0, tab_offset_y=0, tab_shift_x=0,
                        tabShoulder=False, tabShoulderWidth=30, tabShoulderLength=80, tabShoulderRadius=None,
                        flipped=False, rotation=0, bgcolor=None, shunt=False, shunt_width=10, shunt_dist=150, shunt_length=400, shunt_side='left', **kwargs):
    '''
    Local copy from junctionLib so this script can be edited independently.
    '''

    def struct():
        if isinstance(pos, m.Structure):
            return pos
        elif isinstance(pos, tuple):
            return m.Structure(chip, start=pos, direction=rotation)
        else:
            return chip.structure(pos)

    padstart = struct().start

    if padw2 is None:
        padw2 = padw
    if padh2 is None:
        padh2 = padh
    if leadw2 is None:
        leadw2 = leadw
    if leadh2 is None:
        leadh2 = leadh

    if flipped:
        tab_shift_x = tab_shift_x - 10
        toppad = [
            padstart,
            (padstart[0] + padw, padstart[1]),
            (padstart[0] + padw, padstart[1] + padh),
            (padstart[0], padstart[1] + padh),
            padstart
        ]

        botpad = [
            (padstart[0], padstart[1] - separation),
            (padstart[0], padstart[1] - separation - padh2),
            (padstart[0] + padw2, padstart[1] - separation - padh2),
            (padstart[0] + padw2, padstart[1] - separation),
            (padstart[0], padstart[1] - separation)
        ]

        radius = padradius
        p1_top = cornerRound(toppad[0], 3, radius)
        p2_top = cornerRound(toppad[1], 4, radius)
        p3_top = cornerRound(toppad[2], 1, radius)
        p4_top = cornerRound(toppad[3], 2, radius)
        toppad_fillet = p4_top + p3_top + p2_top + p1_top

        p1_bot = cornerRound(botpad[0], 2, radius)
        p2_bot = cornerRound(botpad[1], 3, radius)
        p3_bot = cornerRound(botpad[2], 4, radius)
        p4_bot = cornerRound(botpad[3], 1, radius)
        botpad_fillet = p4_bot + p3_bot + p2_bot + p1_bot

        if tab:
            slotposadjust = JcalcTabDims(chip, struct(), gapw=tab_gapw, gapl=tab_gapl, tabw=tab_tabw, tabl=tab_tabl, taboffs=tab_taboffs, r_out=tab_r_out, r_ins=tab_r_ins, absoluteDimensions=False, stemw=None, steml=None, **kwargs)
            topslotx = padstart[0] + leadw/2 + tab_shift_x
            topslotpos = (topslotx, padstart[1] + slotposadjust[1])
            slot_points_top = [
                (topslotx + slotposadjust[0], padstart[1]),
                (topslotx + slotposadjust[0], padstart[1] + slotposadjust[1]),
                (topslotx - slotposadjust[0], padstart[1] + slotposadjust[1]),
                (topslotx - slotposadjust[0], padstart[1])
            ]
            botslotx = padstart[0] + leadw2/2 + tab_shift_x
            botslotpos = (botslotx, padstart[1] - separation)
            slot_points_bot = [
                (botslotx - slotposadjust[0], padstart[1] - separation),
                (botslotx - slotposadjust[0], padstart[1] - separation - slotposadjust[1]),
                (botslotx + slotposadjust[0], padstart[1] - separation - slotposadjust[1]),
                (botslotx + slotposadjust[0], padstart[1] - separation)
            ]

        if shunt:
            if shunt_side == 'left':
                shunt_start_x = padstart[0]
            else:
                shunt_start_x = padstart[0] + padw

            shunt_start = (shunt_start_x, padstart[1] - separation / 2 - shunt_length / 2)
            shunt_end = (shunt_start_x, padstart[1] - separation / 2 + shunt_length / 2)

            shunt_points_inner = [
                shunt_end,
                (shunt_start[0], shunt_end[1]),
                (shunt_start[0]-shunt_dist, shunt_end[1]),
                (shunt_start[0]-shunt_dist, shunt_start[1]),
                shunt_start
            ]

            shunt_points_outer = [
                (shunt_start[0], shunt_start[1]-shunt_width),
                (shunt_start[0]-shunt_dist-shunt_width, shunt_start[1]-shunt_width),
                (shunt_end[0]-shunt_dist-shunt_width, shunt_end[1]+shunt_width),
                (shunt_end[0], shunt_end[1]+shunt_width)
            ]

            combined_pad = p1_top + shunt_points_inner + p1_bot + p4_bot + p3_bot + p2_bot + shunt_points_outer + p4_top + p3_top + p2_top

            if tab:
                combined_pad = p1_top + shunt_points_inner + p1_bot + slot_points_bot + p4_bot + p3_bot + p2_bot + shunt_points_outer + p4_top + p3_top + p2_top + slot_points_top

            chip.add(SolidPline((0, 0), points=combined_pad, **kwargs))
        else:
            if tab:
                top_pad_points = p4_top + p3_top + p2_top + slot_points_top + p1_top
                bot_pad_points = p4_bot + p3_bot + p2_bot + p1_bot + slot_points_bot
            else:
                top_pad_points = toppad_fillet
                bot_pad_points = botpad_fillet

            chip.add(SolidPline((0, 0), points=top_pad_points, **kwargs))
            chip.add(SolidPline((0, 0), points=bot_pad_points, **kwargs))

        if tab:
            slotposadjust = JcalcTabDims(chip, struct(), gapw=tab_gapw, gapl=tab_gapl, tabw=tab_tabw, tabl=tab_tabl, taboffs=tab_taboffs, r_out=tab_r_out, r_ins=tab_r_ins, absoluteDimensions=False, stemw=None, steml=None, **kwargs)
            topslotx = padstart[0] + leadw/2 - tab_offset_x + tab_shift_x
            topslotpos = (topslotx, padstart[1] + slotposadjust[1])
            botslotx = padstart[0] + leadw2/2 + tab_shift_x
            botslotpos = (botslotx, padstart[1] - separation - tab_offset_y)

            slot_kwargs = {k: v for k, v in kwargs.items() if k not in ['gapw', 'gapl', 'tabw', 'tabl', 'taboffs', 'r_out', 'r_ins']}

            JContact_slot(chip, m.Structure(chip, start=topslotpos, direction=rotation-90), hflip=flipped, rotation=0,
                         gapw=tab_gapw, gapl=tab_gapl, tabw=tab_tabw, tabl=tab_tabl, taboffs=tab_taboffs, r_out=tab_r_out, r_ins=tab_r_ins,
                         **slot_kwargs)
            JContact_slot(chip, m.Structure(chip, start=botslotpos, direction=rotation-90), hflip=not flipped, rotation=0,
                         gapw=tab_gapw, gapl=tab_gapl, tabw=tab_tabw, tabl=tab_tabl, taboffs=tab_taboffs, r_out=tab_r_out, r_ins=tab_r_ins,
                         **slot_kwargs)

            if tabShoulder:
                chip.add(RoundRect(struct().getPos((-tabShoulderLength, tabShoulderWidth / 2)), tabShoulderLength, tabShoulderWidth / 2, min(tabShoulderRadius, (tabShoulderWidth / 2) / 2),
                                   roundCorners=[0, 0, 0, 1], rotation=struct().direction, bgcolor=bgcolor, **kwargs))
                chip.add(RoundRect(struct().getPos((-tabShoulderLength, -tabShoulderWidth / 2)), tabShoulderLength, tabShoulderWidth / 2, min(tabShoulderRadius, (tabShoulderWidth / 2) / 2),
                                   roundCorners=[1, 0, 0, 0], valign=const.TOP, rotation=struct().direction, bgcolor=bgcolor, **kwargs))
    else:
        print('You chose not flipped. This is not yet written, sorry.')

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
                
                # Transmon3D with shunt (scaled to fit 500 micron field)
                Transmon3DWithShunt(self, (cx + CHIPLET_SIZE_x/2 - TMON_PAD_WIDTH/2, cy + CHIPLET_SIZE_y/2 + TMON_SEPARATION/2), 
                                    padw=TMON_PAD_WIDTH, padh=TMON_PAD_HEIGHT, 
                                    leadw=TMON_LEAD_WIDTH, leadh=TMON_LEAD_HEIGHT, 
                                    padradius=TMON_PAD_RADIUS, 
                                    separation=TMON_SEPARATION, shunt=True, 
                                    shunt_width=TMON_SHUNT_WIDTH, shunt_dist=TMON_SHUNT_DIST, 
                                    shunt_length=TMON_SHUNT_LENGTH, shunt_side='left', flipped=True, 
                                    tab=True, tab_offset_x=CHIPLET_SIZE_x, tab_offset_y=CHIPLET_SIZE_y, tab_shift_x=TMON_PAD_WIDTH/2,
                                    layer=wafer.lyr(METAL_LAYER))
        
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
        # Temporarily override wafer chip dimensions so Chip sets the correct origin offset.
        original_frame = wafer.frame
        original_chip_x = wafer.chipX
        original_chip_y = wafer.chipY
        wafer.frame = False
        wafer.chipX = CORNER_CHIP_SIZE + wafer.sawWidth
        wafer.chipY = CORNER_CHIP_SIZE + wafer.sawWidth
        m.Chip.__init__(self, wafer, chipID, layer, defaults={'w':200, 'r_out':10, 'r_ins':0})
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
                Transmon3DWithShunt(self, (cx + CORNER_CHIP_SIZE/2 - TMON_PAD_WIDTH/2, cy + CORNER_CHIP_SIZE/2 + TMON_SEPARATION/2), 
                                    padw=TMON_PAD_WIDTH, padh=TMON_PAD_HEIGHT, 
                                    leadw=TMON_LEAD_WIDTH, leadh=TMON_LEAD_HEIGHT, 
                                    padradius=TMON_PAD_RADIUS, 
                                    separation=TMON_SEPARATION, shunt=True, 
                                    shunt_width=TMON_SHUNT_WIDTH, shunt_dist=TMON_SHUNT_DIST, 
                                    shunt_length=TMON_SHUNT_LENGTH, shunt_side='left', flipped=True, 
                                    tab=True, tab_offset_x=CORNER_CHIP_SIZE, tab_offset_y=CORNER_CHIP_SIZE, tab_shift_x=TMON_PAD_WIDTH/2,
                                    layer=wafer.lyr(METAL_LAYER))
        

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