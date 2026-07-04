#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 13:52:03 2021

@author: sasha
Edited by Agrim, 2026 (junction stack, centered tabs, corner-based
coordinates, 3D parameter sweep via maskLib.arrayLib)
"""

import time
SCRIPT_START = time.time()

import numpy as np

import maskLib.MaskLib as m
from maskLib.arrayLib import Sweep3D, dose_layer, export_ldt
from maskLib.junctionLib import setupJunctionLayers, JcalcTabDims, JContact_slot, Transmon3DWithShunt
from maskLib.fluxoniumLib import smallJJ, leads_for_tmon_dosearray_custom
from maskLib.Entities import SolidPline, RoundRect
from maskLib.markerLib import MarkerSquare, MarkerCross
from maskLib.utilities import doMirrored, cornerRound
from dxfwrite import DXFEngine as dxf
from dxfwrite import const

# ===============================================================================
# 1) Wafer & field-array geometry
# ===============================================================================
FIELD_PADDING = 500     # Padding inside each chiplet (keeps fields away from the edge)
FIELD_SIZE = 500        # Field size in microns
FIELD_SAW = 0           # Spacing between fields (saw width) in microns
CHIPLET_SIZE_x = 21000 + FIELD_SAW      # Large chiplet dimensions in microns
CHIPLET_SIZE_y = 21000 + FIELD_SAW      # Large chiplet dimensions in microns
WAFER_DIAMETER_NAME = '4in'             # key into m.waferDiameters (also used in the optical-only file name)

# field grid dimensions (derived): 40 x 40 for the default sizes above
FIELD_STEP = FIELD_SIZE + FIELD_SAW
GRID_NX = int((CHIPLET_SIZE_x - 2 * FIELD_PADDING) / FIELD_STEP)
GRID_NY = int((CHIPLET_SIZE_y - 2 * FIELD_PADDING) / FIELD_STEP)

# corner chips (smaller chips placed at the wafer corners)
CORNER_CHIP_SIZE = int((CHIPLET_SIZE_x-FIELD_SAW-2*FIELD_PADDING)//2 + 2*FIELD_PADDING)
CORNER_GRID_N = int((CORNER_CHIP_SIZE - 2 * FIELD_PADDING) / FIELD_STEP)  # fields per side (20)

FRAME_LAYER = '703/0'      # Layer for chip frame boundary
METAL_LAYER = 'BASEMETAL'  # Layer for metal structures

# ------------------------------------------------------------------------------
# Workflow toggles -- the two bools below pick what this run generates
# ------------------------------------------------------------------------------
# OPTICAL_ONLY = True : draw the FULL WAFER (all chiplets + corner chips
#                       placed) with optical layers only -- no ebeam
#                       structures, no standalone chip DXFs, no xlsx.
#                       For wafer-level layout checks. Output name gets an
#                       '_opticalOnly' suffix.
# OPTICAL_ONLY = False: no full wafer; generate ONE standalone chip DXF with
#                       all ebeam structures plus its sweep xlsx. Which chip
#                       is picked by the second toggle:
#     GENERATE_CORNER_CHIP = False: the main chiplet (40x40 fields)
#     GENERATE_CORNER_CHIP = True : the corner chip (20x20 fields) with its
#                                   own 20x20 xlsx
OPTICAL_ONLY = False
GENERATE_CORNER_CHIP = False

REUSE_IDENTICAL_CHIPS = True       # Reuse one generated block for repeated chiplets/corner chips

# derived output flags -- set the two workflow toggles above, not these
RENDER_FULL_WAFER = OPTICAL_ONLY
DRAW_EBEAM_LAYERS = not OPTICAL_ONLY
EXPORT_SAMPLE_CHIPLET_DXF = (not OPTICAL_ONLY) and (not GENERATE_CORNER_CHIP)
EXPORT_CORNER_CHIP_DXF = (not OPTICAL_ONLY) and GENERATE_CORNER_CHIP
EXPORT_SWEEP_MAP = not OPTICAL_ONLY

# ===============================================================================
# 2) Transmon & junction parameters (defaults; the sweep overrides per field)
# ===============================================================================
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
# Each lead starts with a 'home plate' contact pad sitting on the transmon
# pad edge, covering the contact tab slot cut into the pad; a wedge tapers
# it down to the thin lead that runs to the junction.
JJ_LEAD_WIDTH = 1            # Width of thin leads in microns
JJ_CONTACT_WIDTH = 20        # Width of contact pads (home plates) in microns
JJ_CONTACT_LENGTH = 10       # Length of contact pads in microns
JJ_WEDGE_LENGTH = 10         # Length of wedge taper from contact pad to thin lead in microns
JJ_PLATE_MARGIN = 1.5        # How far the home plate reaches past the tab slot into the pad in microns
JJ_FINGER_LENGTH = 1.5       # Length of small/big fingers in microns
JJ_SMALLFINGER_WIDTH = 0.140 # Width of small finger in microns
JJ_BIGFINGER_WIDTH = 0.340   # Width of big finger in microns (small finger + 0.2)
JJ_BRIDGE_WIDTH = 0.840      # Width of bridge in microns (small finger + 0.7)
JJ_BRIDGE_LENGTH = 0.250     # Length of bridge in microns
JJ_UNDERCUT = 0.2            # Undercut width in microns

# ===============================================================================
# 3) Parameter sweep configuration (3D: tile columns x tile rows x tiles)
# ===============================================================================
# The field grid is tiled with identical TILE_NX x TILE_NY blocks ("tiles").
# Each sweep axis is a dict {parameter_name: (start, step)} -- the number of
# steps comes from the axis (TILE_NX, TILE_NY, or the number of tiles) and
# the final value is calculated automatically. Several parameters can share
# one axis (they sweep in lockstep). Use {} to disable an axis.
#
# Geometry parameters change the drawn shapes; dose parameters redirect
# shapes onto auto-generated per-dose layers (e.g. BRIDGE_400) -- never add
# those to SetupLayers by hand.
#
# Every field gets a corner label like 'A0103' (tile A, column 01, row 03).
# The label -> parameter mapping is drawn as a legend in the chip margin and
# exported to <wafer>_sweep_map.xlsx (parameter table, label minimap, and a
# gradient-colored value map per swept parameter).

TILE_NX = 10   # tile size in fields; must evenly divide the chiplet field grid
TILE_NY = 20

SWEEP_COL = {                                # 10 steps along x in each tile
    'smallfinger_width': (0.100, 0.010),     # 0.10 -> 0.19 um
    'bigfinger_width':   (0.300, 0.010),     # locked to smallfinger + 0.2
}
SWEEP_ROW = {                                # 20 steps along y in each tile
    'bridge_dose': (400, 40),                # 400 -> 1160
}
SWEEP_TILE = {                               # one step per tile (8 tiles A-H)
    'smallfinger_dose': (800, 100),          # 800 -> 1500
}

# which JJ geometry kwarg each geometry parameter controls
GEOMETRY_PARAMS = {
    'smallfinger_width': 'smallfingerW',
    'bigfinger_width':   'bigfingerW',
    'bridge_width':      'bridgeW',
    'bridge_length':     'bridgeL',
    'finger_length':     'fingerL',
}
# which layer family each dose parameter redirects (layer = FAMILY_<dose>)
DOSE_PARAMS = {
    'pads_dose':        'BASEMETAL',
    'leads_dose':       'LEADS',
    'smallfinger_dose': 'SMALLFINGER',
    'bigfinger_dose':   'BIGFINGER',
    'bridge_dose':      'BRIDGE',
    'undercut_dose':    'UNDERCUT',
    'shift_dose':       'SHIFT',
}

# doses for ebeam layers whose dose is NOT being swept -- written to the
# Elionix .ldt dose table alongside the swept per-dose layers. Families that
# ARE swept (e.g. BRIDGE while bridge_dose is a sweep axis) are ignored here.
EBEAM_BASE_DOSES = {
    'LEADS':       1200,
    'SMALLFINGER': 1000,
    'BIGFINGER':   1000,
    'BRIDGE':      500,
    'UNDERCUT':    200,
    'SHIFT':       400,
    'LABELS':      600,   # field labels + legend (large features, lower dose)
}

sweep = Sweep3D(GRID_NX, GRID_NY, TILE_NX, TILE_NY,
                col=SWEEP_COL, row=SWEEP_ROW, tile=SWEEP_TILE,
                geometry_params=GEOMETRY_PARAMS, dose_params=DOSE_PARAMS)
sweep.print_summary()

# Output file names are built from the design parameters, so the DXF and
# xlsx names always describe what was generated. Only DESIGN_NAME is
# hardcoded (change it freely); everything after it is auto-populated.
DESIGN_NAME = 'ShuntedJJArray'
if OPTICAL_ONLY:
    # Full-wafer optical layout check: the per-chiplet sweep detail (grid,
    # tile counts, sweep type) is irrelevant here, so name it by the wafer
    # diameter only, e.g. 'ShuntedJJArray_4in_opticalOnly'
    WAFER_NAME = '%s_%s_opticalOnly' % (DESIGN_NAME, WAFER_DIAMETER_NAME)
else:
    # Standalone chiplet: describe the sweep as steps-per-axis in
    # col x row x tile order, matching the sweep_type words,
    # e.g. 'ShuntedJJArray_21mm_10x20x8_SizeDoseDose'
    WAFER_NAME = '%s_%gmm_%dx%dx%d_%s' % (
        DESIGN_NAME, CHIPLET_SIZE_x / 1000,
        sweep.tile_nx, sweep.tile_ny, sweep.n_tiles, sweep.sweep_type())
print('Output name:', WAFER_NAME)

# ===============================================================================
# 4) Drawing functions
# ===============================================================================
def JunctionWithLeads(chip, pos, params=None):
    '''
    Draw the junction stack in the gap between the transmon pads: from each
    pad a 'home plate' contact pad covering the pad's contact tab slot,
    tapering through a wedge into a thin lead that meets a Dolan junction
    (small finger, bridge, big finger, with undercut and shift layers)
    centered in the gap.

    pos is the same point passed to Transmon3DWithShunt (bottom-left corner
    of the top pad). Dimensions default to the JJ_* design parameters;
    params (from sweep.field) overrides geometry values and redirects
    shapes onto per-dose layers (e.g. BRIDGE_400).
    '''
    if params is None:
        params = {}

    # geometry: JJ_* defaults, overridden by any swept geometry parameters
    geo = {'smallfingerW': JJ_SMALLFINGER_WIDTH, 'bigfingerW': JJ_BIGFINGER_WIDTH,
           'bridgeW': JJ_BRIDGE_WIDTH, 'bridgeL': JJ_BRIDGE_LENGTH,
           'fingerL': JJ_FINGER_LENGTH}
    for pname, key in GEOMETRY_PARAMS.items():
        if pname in params:
            geo[key] = params[pname]

    total_JJ_length = 2*geo['fingerL'] + geo['bridgeL']

    # The home plate sits on the pad edge, covering the contact tab slot cut
    # into the pad (dims from JcalcTabDims, same defaults Transmon3DWithShunt
    # uses for the tab): its flat end reaches JJ_PLATE_MARGIN past the slot
    # into the pad, the rest sticks into the gap and tapers via the wedge to
    # the thin lead running to the junction.
    slot_depth, slot_halfwidth = JcalcTabDims(chip, pos)
    plate_start = slot_depth + JJ_PLATE_MARGIN  # lead start, measured into the pad from its edge
    assert JJ_CONTACT_WIDTH/2 > slot_halfwidth, 'Home plate too narrow to cover the tab slot: increase JJ_CONTACT_WIDTH'
    assert JJ_CONTACT_LENGTH > plate_start, 'Home plate too short to cover the tab slot: increase JJ_CONTACT_LENGTH or reduce JJ_PLATE_MARGIN'
    wedge_to_JJ = TMON_SEPARATION/2 + plate_start - JJ_CONTACT_LENGTH - JJ_WEDGE_LENGTH - total_JJ_length/2
    assert wedge_to_JJ > 0, 'Leads do not fit: shrink contact/wedge lengths or increase TMON_SEPARATION'

    # chip.add() shifts objects with a .points attribute (like the transmon's
    # SolidPline pads) by chip.origin_offset, but not the plain
    # polylines/rectangles drawn here -- apply the offset manually so both
    # always land in the same place. With centerChip=False this is (0,0).
    x0 = pos[0] + chip.origin_offset[0]
    y0 = pos[1] + chip.origin_offset[1]

    xc = x0 + TMON_PAD_WIDTH/2  # center line of the pads

    for side, ystart in (('top', y0 + plate_start),
                         ('bottom', y0 - TMON_SEPARATION - plate_start)):
        leads_for_tmon_dosearray_custom(chip, m.Structure(chip, start=(xc, ystart), direction=0),
                                        toporbottom=side,
                                        layer=dose_layer('LEADS', params, 'leads_dose'),
                                        leadLpadtocontact=0, leadLcontacttoJJ=wedge_to_JJ,
                                        leadW=JJ_LEAD_WIDTH, contactW=JJ_CONTACT_WIDTH,
                                        contactL=JJ_CONTACT_LENGTH, wedgeL=JJ_WEDGE_LENGTH)

    # junction centered vertically in the gap (x nudged by leadW/2 to counter
    # the perpendicular lead offset smallJJ applies internally)
    smallJJ(chip, m.Structure(chip, start=(xc + JJ_LEAD_WIDTH/2, y0 - TMON_SEPARATION/2), direction=90),
            smallfingerlayer=dose_layer('SMALLFINGER', params, 'smallfinger_dose'),
            bigfingerlayer=dose_layer('BIGFINGER', params, 'bigfinger_dose'),
            bridgelayer=dose_layer('BRIDGE', params, 'bridge_dose'),
            Ulayer=dose_layer('UNDERCUT', params, 'undercut_dose'),
            Slayer=dose_layer('SHIFT', params, 'shift_dose'),
            gap=geo['bridgeW'], leadW=JJ_LEAD_WIDTH, fingerL=geo['fingerL'],
            bigfingerW=geo['bigfingerW'], smallfingerW=geo['smallfingerW'],
            bridgeW=geo['bridgeW'], bridgeL=geo['bridgeL'], undercut=JJ_UNDERCUT)


def draw_field(chip, wafer, cx, cy, params, flabel):
    '''one complete field at center (cx, cy): frame, transmon, junction, label'''
    # Field frame on 703/0 layer
    chip.add(dxf.rectangle(
        (cx - FIELD_SIZE/2, cy - FIELD_SIZE/2),
        FIELD_SIZE, FIELD_SIZE,
        layer=wafer.lyr(FRAME_LAYER)
    ))

    # Transmon3D with shunt (scaled to fit 500 micron field)
    tpos = (cx - TMON_PAD_WIDTH/2, cy + TMON_SEPARATION/2)
    Transmon3DWithShunt(chip, tpos,
                        padw=TMON_PAD_WIDTH, padh=TMON_PAD_HEIGHT,
                        leadw=TMON_LEAD_WIDTH, leadh=TMON_LEAD_HEIGHT,
                        padradius=TMON_PAD_RADIUS,
                        separation=TMON_SEPARATION, shunt=True,
                        shunt_width=TMON_SHUNT_WIDTH, shunt_dist=TMON_SHUNT_DIST,
                        shunt_length=TMON_SHUNT_LENGTH, shunt_side='left', flipped=True,
                        tab=True, tab_shift_x=TMON_PAD_WIDTH/2 - TMON_LEAD_WIDTH/2,
                        layer=wafer.lyr(dose_layer(METAL_LAYER, params, 'pads_dose')))

    # Leads and Dolan junction in the pad gap (ebeam structures)
    if DRAW_EBEAM_LAYERS:
        JunctionWithLeads(chip, tpos, params)

    # field label (tile letter + column + row, e.g. A0103) in the
    # bottom-left corner, clear of the pads and shunt. LABELS is an ebeam
    # layer, so skip it (like the junctions) on the optical-only wafer.
    if DRAW_EBEAM_LAYERS:
        chip.add_chip_label(flabel,
                            (cx - FIELD_SIZE/2 + 50, cy - FIELD_SIZE/2 + 42),
                            height=22, layer='LABELS')


# ===============================================================================
# wafer setup
# ===============================================================================

w = m.Wafer(WAFER_NAME,'DXF/',CHIPLET_SIZE_x,CHIPLET_SIZE_y,padding=1500,waferDiameter=m.waferDiameters[WAFER_DIAMETER_NAME],sawWidth=500,singleChipColumn=False, centerChip=False, frame=True, markers=False
)
    #set wafer properties
    # w.frame: draw frame layer?
    # w.solid: draw things solid?
    # w.multiLayer: draw in multiple layers?
    # w.singleChipColumn: only make one column of chips?

# ------------------------------------------------------------------------------
# Layers, grouped by fabrication role
# ------------------------------------------------------------------------------
# INACTIVE  -- DXF bookkeeping only, never fabricated: '0' and 'VIEWPORTS'
#              (managed by dxfwrite; layer-0 content inside blocks inherits
#              the chip's layer at insert time)
# GUIDE     -- drawn for reference, written to no mask: field/chip frames on
#              703/0, plus WAFER_FRAME (wafer outline circles, added
#              automatically by w.init() because frame=True)
# OPTICAL 1 -- first photolithography exposure
# OPTICAL 2 -- second photolithography exposure
#              (Opt_Mark alignment marks are shared by both optical steps)
# EBEAM     -- electron-beam exposure: the junction structures; per-dose
#              copies (e.g. BRIDGE_400) are auto-generated from the sweep
#
# NOTE: the first layer passed to SetupLayers becomes the wafer default
# layer, so LAYERS_OPTICAL_1 (BASEMETAL) must stay first.

LAYERS_GUIDE = [
    ['703/0', 9],          # field + chip frame boundaries (light gray)
]
LAYERS_OPTICAL_1 = [
    ['BASEMETAL', 5],      # transmon pads + shunts
    ['DICEBORDER', 4],     # dicing saw guides
    ['Opt_Mark', 3],       # optical alignment marks (also used in step 2)
]
LAYERS_OPTICAL_2 = [
    ['TiW_Mark', 6],       # TiW marker crosses
]
LAYERS_EBEAM = [
    ['LEADS', 1],          # base (undosed) ebeam layers; swept structures
    ['SMALLFINGER', 30],   # move to the per-dose copies below
    ['BIGFINGER', 32],
    ['BRIDGE', 34],
    ['UNDERCUT', 36],
    ['SHIFT', 38],
    ['LABELS', 2],         # field labels + sweep legend (ebeam: chiplet DXFs only, not the optical-only wafer)
]
# per-dose ebeam layers from the sweep (e.g. BRIDGE_400 ... SMALLFINGER_1500)
LAYERS_EBEAM_DOSES = [[name, 50 + k % 200] for k, name in enumerate(sweep.dose_layers())]

if not DRAW_EBEAM_LAYERS:
    LAYERS_EBEAM = []
    LAYERS_EBEAM_DOSES = []

w.SetupLayers(LAYERS_OPTICAL_1 + LAYERS_OPTICAL_2 + LAYERS_GUIDE
              + LAYERS_EBEAM + LAYERS_EBEAM_DOSES)

# (index l[0]: SetupLayers appends a layer number to each entry in place)
print('Layer manifest by fabrication role:')
print('  inactive (never written): 0, VIEWPORTS')
print('  guide    (not written):   WAFER_FRAME, ' + ', '.join(l[0] for l in LAYERS_GUIDE))
print('  optical 1: ' + ', '.join(l[0] for l in LAYERS_OPTICAL_1))
print('  optical 2: ' + ', '.join(l[0] for l in LAYERS_OPTICAL_2) + '  (+ Opt_Mark alignment)')
if DRAW_EBEAM_LAYERS:
    print('  ebeam:     ' + ', '.join(l[0] for l in LAYERS_EBEAM)
          + ' + %d per-dose layers' % len(LAYERS_EBEAM_DOSES))
else:
    print('  ebeam:     DISABLED (DRAW_EBEAM_LAYERS = False)')

#initialize the wafer (remember to finalize any wafer properties like layers before initializing!)
w.init()


#do dicing border (by default located on layer 'MARKERS', so let's put it on layer 'DICEBORDER' instead)
w.DicingBorder(layer='DICEBORDER')

#do optical markers
#(note: mirrorX and mirrorY are true by default, but I've exposed them here to demonstrate how they work)
doMirrored(MarkerCross, w, (0,45000),(500,500), 20,layer='Opt_Mark',mirrorX=True,mirrorY=True)
doMirrored(MarkerCross, w, (45000,0),(500,500), 20,layer='Opt_Mark',mirrorX=True,mirrorY=True)

# ===============================================================================
# chiplet class definition (named for its size: 21000 um = 21 mm square)
# ===============================================================================
class Chiplet21000um(m.Chip):
    def __init__(self, wafer, chipID, layer, defaults=None, **kwargs):
        # centerChip=False: chip.add() applies no origin shift, so all entity
        # types share one corner-based coordinate system (no +chipsize/2 hacks)
        m.Chip.__init__(self, wafer, chipID, layer, defaults={'w':200, 'r_out':10, 'r_ins':0}, centerChip=False)
        if defaults is not None:
            for d in defaults:
                self.defaults[d] = defaults[d]

        print(f"Fields per chiplet: {GRID_NX} x {GRID_NY} ({GRID_NX*GRID_NY} total)")

        # Center the field grid in the large chiplet
        offset_x = (CHIPLET_SIZE_x - GRID_NX * FIELD_STEP) / 2 + FIELD_STEP / 2
        offset_y = (CHIPLET_SIZE_y - GRID_NY * FIELD_STEP) / 2 + FIELD_STEP / 2

        for ix in range(GRID_NX):
            for iy in range(GRID_NY):
                params, flabel = sweep.field(ix, iy)
                draw_field(self, wafer,
                           offset_x + ix * FIELD_STEP, offset_y + iy * FIELD_STEP,
                           params, flabel)

        # sweep legend in the bottom margin, below the field grid
        # (LABELS is an ebeam layer: skipped on the optical-only wafer)
        if DRAW_EBEAM_LAYERS:
            for k, line in enumerate(sweep.legend_lines()):
                self.add_chip_label(line, (CHIPLET_SIZE_x/2, 400 - 110*k),
                                    height=70, layer='LABELS')

        # sweep workbook for the lab notebook (parameter table, label
        # minimap, and a gradient-colored value map per swept parameter)
        if EXPORT_SWEEP_MAP:
            xlsx_path = wafer.path + wafer.fileName + '_sweep_map.xlsx'
            sweep.export_workbook(xlsx_path)
            print('Sweep map saved to', xlsx_path)

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
# Corner chip class definition (11000 um square chips for the wafer corners)
# ===============================================================================
class CornerChip11000um(m.Chip):
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

        mx = my = CORNER_GRID_N
        print(f"Corner chip size: {CORNER_CHIP_SIZE} microns")
        print(f"Fields per corner chip: mx={mx}, my={my} ({mx*my} total)")

        # Center the field grid in the corner chip with step offset (same as SimpleChiplet)
        offset_x = (CORNER_CHIP_SIZE - mx * FIELD_STEP) / 2 + FIELD_STEP / 2
        offset_y = (CORNER_CHIP_SIZE - my * FIELD_STEP) / 2 + FIELD_STEP / 2

        for ix in range(mx):
            for iy in range(my):
                # strict=False: the corner chip grid may hold fewer tiles
                # than the main chiplet (tile parameter values wrap around)
                params, flabel = sweep.field(ix, iy, grid_nx=mx, grid_ny=my, strict=False)
                draw_field(self, wafer,
                           offset_x + ix * FIELD_STEP, offset_y + iy * FIELD_STEP,
                           params, flabel)

        # sweep legend in the bottom margin
        # (LABELS is an ebeam layer: skipped on the optical-only wafer)
        if DRAW_EBEAM_LAYERS:
            for k, line in enumerate(sweep.legend_lines()):
                self.add_chip_label(line, (CORNER_CHIP_SIZE/2, 400 - 110*k),
                                    height=70, layer='LABELS')

        # corner-chip sweep workbook, covering only its own grid (20x20;
        # tile parameter values wrap around, matching what is drawn)
        if EXPORT_SWEEP_MAP and GENERATE_CORNER_CHIP:
            xlsx_path = wafer.path + wafer.fileName + '_CORNER_%dx%d_sweep_map.xlsx' % (mx, my)
            sweep.export_workbook(xlsx_path, grid_nx=mx, grid_ny=my, strict=False)
            print('Corner sweep map saved to', xlsx_path)


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
# generate chips (what gets built is set by the workflow toggles up top)
# ===============================================================================
# chip IDs carry each chip's own size and field grid (the wafer-name prefix
# describes the main chiplet design; this suffix describes the chip itself)

def export_ldt_for(fileName):
    '''Elionix layer dose table for one chip's ebeam job, named to match its
    DXF. Layer numbers come from the wafer layer table PLUS ONE: the DXF->GDS
    converter numbers layers 1-based (verified against a converted GDS:
    BASEMETAL table index 1 -> GDS layer 2, LEADS 7 -> 8, SHIFT 12 -> 13).
    Doses come from the sweep values plus EBEAM_BASE_DOSES for unswept
    families.'''
    path = w.path + fileName + '.ldt'
    export_ldt(path, sweep.ldt_entries(lambda name: w.layerNums[name] + 1, EBEAM_BASE_DOSES))
    print('Elionix dose table saved to', path)

# --- main chiplet: needed for the full wafer, or when it is the requested chip
if OPTICAL_ONLY or not GENERATE_CORNER_CHIP:
    default_chiplet = Chiplet21000um(w, f'CHIPLET_{CHIPLET_SIZE_x}um_{GRID_NX}x{GRID_NY}', w.defaultLayer)
    w.setDefaultChip(default_chiplet)

    # Populate remaining chiplets in buffer
    if REUSE_IDENTICAL_CHIPS:
        for i in range(1, len(w.chips)):
            w.setChipBuffer(default_chiplet, i)
    else:
        for i in range(1, len(w.chips)):
            w.setChipBuffer(Chiplet21000um(w, f'CHIPLET_{CHIPLET_SIZE_x}um_{GRID_NX}x{GRID_NY}_{i}', w.defaultLayer).save(w), i)

    # Save a standalone chiplet DXF + its ebeam dose table if requested
    # (named exactly WAFER_NAME: only one chip DXF is generated per run, and
    # the wafer name already fully describes it)
    if EXPORT_SAMPLE_CHIPLET_DXF and len(w.chips) > 1:
        w.chips[1].save(w, drawCopyDXF=True, dicingBorder=False, fileName=w.fileName)
        export_ldt_for(w.fileName)

# --- corner chips: needed for the full wafer, or when they are the requested chip
if OPTICAL_ONLY or GENERATE_CORNER_CHIP:
    # Calculate corner positions (at ~55% of wafer radius, at 45 degrees)
    wafer_radius = m.waferDiameters[WAFER_DIAMETER_NAME] / 2
    corner_distance = wafer_radius * 0.55

    corner_positions = [
        (-corner_distance, -corner_distance),  # Bottom-left
        (corner_distance, -corner_distance),   # Bottom-right
        (-corner_distance, corner_distance),   # Top-left
        (corner_distance, corner_distance),    # Top-right
    ]

    # corner-chip DXF/ldt name: wafer name + corner-chip size (the wafer
    # name's own size describes the main chiplet)
    corner_name = w.fileName + '_CORNER_%gmm' % (CORNER_CHIP_SIZE / 1000)

    if REUSE_IDENTICAL_CHIPS:
        corner_chip = CornerChip11000um(w, f'CORNER_{CORNER_CHIP_SIZE}um_{CORNER_GRID_N}x{CORNER_GRID_N}', w.defaultLayer)
        corner_chip.save(w, drawCopyDXF=EXPORT_CORNER_CHIP_DXF, dicingBorder=False, fileName=corner_name)
        if EXPORT_CORNER_CHIP_DXF:
            export_ldt_for(corner_name)

        if RENDER_FULL_WAFER:
            for cx, cy in corner_positions:
                # Adjust insertion point to compensate for CORNER_CHIP_SIZE/2 shift
                adj_x = cx - CORNER_CHIP_SIZE / 2 - FIELD_SIZE
                adj_y = cy - CORNER_CHIP_SIZE / 2 - FIELD_SIZE
                insert_pt = w.chipSpace((adj_x, adj_y))
                w.drawing.add(dxf.insert(corner_chip.ID, insert=insert_pt, layer=w.lyr(corner_chip.layer)))
    else:
        for idx, (cx, cy) in enumerate(corner_positions):
            corner_chip = CornerChip11000um(w, f'CORNER_{CORNER_CHIP_SIZE}um_{CORNER_GRID_N}x{CORNER_GRID_N}_{idx}', w.defaultLayer)
            corner_chip.save(w, drawCopyDXF=EXPORT_CORNER_CHIP_DXF and idx == 0, dicingBorder=False, fileName=corner_name)
            if EXPORT_CORNER_CHIP_DXF and idx == 0:
                export_ldt_for(corner_name)
            # Adjust insertion point to compensate for CORNER_CHIP_SIZE/2 shift
            adj_x = cx - CORNER_CHIP_SIZE / 2 - FIELD_SIZE
            adj_y = cy - CORNER_CHIP_SIZE / 2 - FIELD_SIZE
            insert_pt = w.chipSpace((adj_x, adj_y))
            if RENDER_FULL_WAFER:
                w.drawing.add(dxf.insert(corner_chip.ID, insert=insert_pt, layer=w.lyr(corner_chip.layer)))

# --- full wafer DXF (OPTICAL_ONLY mode)
if RENDER_FULL_WAFER:
    w.populate()
    w.save()

print('Total runtime: %.1f s' % (time.time() - SCRIPT_START))
