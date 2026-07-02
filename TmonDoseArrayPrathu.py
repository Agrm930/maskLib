#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 13:52:03 2021

Author: Tom
"""

# relative import (i.e. if you're running this script from the maskLib\Localonly directory or similar)
###relative imports to work on S: drive
import os
import sys
import numpy as np

# Print the current working directory
print('Current directory before change:', os.getcwd())

# Change the directory
os.chdir(r'c:\Users\Agrim\OneDrive\Documents\GitHub\maskLib-master')

# Print the current working directory again to see the change
print('Current directory after change:', os.getcwd())

# Add the parent directory to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Add the grandparent directory to sys.path to ensure maskLib can be found
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# Absolute masklib import (i.e. if you've installed it with pip)
import maskLib.MaskLib as m
from maskLib.microwaveLib import Strip_straight
from maskLib.Entities import SolidPline
from maskLib.utilities import cornerRound
from dxfwrite.entities import Polyline

from maskLib.junctionLib import FlagPads, Transmon3DWithShunt, DolanJunction, JContact_slot, JContact_tab, JcalcTabDims
from maskLib.fluxoniumLib import junction_chain, smallJ, smallJJ, no_loop_leads, leads_for_tmon_dosearray_custom
from maskLib.qubitLib import Transmon3D, qubit_defaults

from maskLib.markerLib import MarkerSquare, MarkerCross
from maskLib.utilities import doMirrored

import ezdxf
import dxfwrite
from dxfwrite import DXFEngine as engine
from ezdxf.addons import text2path
#from ezdxf.fonts import FontFace
import time
import sys


# ===============================================================================
# wafer setup
# ===============================================================================

def largest_square_in_circle(radius):
    # The side length of the largest square that can be inscribed in a circle
    # is equal to the diameter of the circle divided by the square root of 2.
    side_length = radius * (2 ** 0.5)
    return side_length

chip_square_side_length = largest_square_in_circle(22000)
# chip_square_side_length = 15000


w = m.Wafer('TmonDA','DXF/',chip_square_side_length,chip_square_side_length,padding=2500,waferDiameter=m.waferDiameters['2in'],sawWidth=300,
                frame=True,solid=False,multiLayer=True,singleChipRow=True,singleChipColumn=True)

w.SetupLayers([
    # ['BASEMETAL',3],
    # ['FT', 1],
    ['PADS',3],
    ['DICEBORDER',1],
    ['MARKERS',2],
    # ['MARKERS2',3],
    # ['MARKERS3',4],
    ['DOSEARRAY',4],
    
    ])

print("dims of chip = ", w.chipY, 'x', w.chipX, 'microns')

#initialize the wafer (remember to finalize any wafer properties like layers before initializing!)
w.init()

#do dicing border (by default located on layer 'MARKERS', so let's put it on layer 'DICEBORDER' instead)
# w.DicingBorder(layer='DICEBORDER')
#wafer label
# w.add_text("Wafer Label", (1000, 1000), 250, "WAFERLABEL", layer='TEXT')
# #do optical markers
# doMirrored(MarkerCross, w, (30000,20000),(200,200), 5,layer='MARKERS1',mirrorX=True,mirrorY=True)


#do ebeam markers--done below instead
# markerpts = [(15000,15000),(14000,14000),(13000,13000),(12000,12000)]
radius = 22000
markerpts = [(0,radius), (radius,0), (radius*np.cos(np.pi/4),radius*np.cos(np.pi/4)), (0,0)]
for pt in markerpts:
    #(note: mirrorX and mirrorY are true by default)
    doMirrored(MarkerSquare, w, pt, 20,layer='MARKERS')

arraychip = w.chips[0]
arraychip.chipID = 'TmonDA_mainchip'
arraychip.ID = 'TmonDA_mainchip'


class TmonDoseArray(m.Chip):
    def __init__(self, wafer, chipID, layer, chip_id_loc=(6100, 0), arraydims=(10,10), arrayspacing_x=2050, arrayspacing_y=2000, transmon_number_label=False, 
                 defaults=None, JJparams=None, JJparams_label=False, Doses=None, Doselabels=False, **kwargs):
        m.Chip.__init__(self, wafer, chipID, layer, defaults={'w': 200, 'r_out': 10, 'r_ins': 0})
        startpoint = (1500, 2350)
        if defaults is not None:
            for d in defaults:
                self.defaults[d] = defaults[d]
        self.add_JJ_dose_array(arraydims=arraydims, arrayspacing_x=arrayspacing_x, arrayspacing_y=arrayspacing_y, basedose=1000,printdose=False,layer='DOSEARRAY', FT=False, 
                               transmon_number_label=transmon_number_label, JJparams=JJparams, JJparams_label=JJparams_label, Doses=Doses, Doselabels=Doselabels)
        
    def add_JJ_dose_array(self, arraydims=(5,5), arrayspacing_x=500, arrayspacing_y=1000, basedose=1000,printdose=False,layer='DOSEARRAY', FT=False,
                          transmon_number_label=False, JJparams=None, JJparams_label=False, Doses=None, Doselabels=False):
        
        def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█', print_end="\r"):
            percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
            filled_length = int(length * iteration // total)
            bar = fill * filled_length + '-' * (length - filled_length)
            print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=print_end)
            if iteration == total:
                print()

        total_iterations = arraydims[0] * arraydims[1]
        iteration = 0

        # Add a dose array to the chip
        # if doses is None:
        #     doses = np.ones(arraydims)*basedose
        for i in range(arraydims[0]):
            for j in range(arraydims[1]): 

                iteration += 1
                print_progress_bar(iteration, total_iterations, prefix='Progress:', suffix='Complete', length=50)

                params = {
                    'startpoint': (500,1000),
                    'large_rect_length': 0,
                    'large_rect_width': 0,
                    'small_rect_length': 500,
                    'small_rect_width': 100,
                    'conductor_width': 10,
                    'Y-offset': 0,  # smallrect offset from the +Y side of the largerect
                    'X-offset': 0,  # offset of the +Y largerect line from the +X end of the small rectangle
                    'outer_radius': 20,
                    'inner_radius': 10
                }
                if FT: # params for "FT" and transmon 3D pads
                    FT_start = (params['startpoint'][0] + 75, params['startpoint'][1] - 0.5*params['small_rect_width'])
                    # Generate dummy FT coordinates using params variables. The dummy FT is just a filleted rectangle.
                    dummy_FT_pts_outer = [
                        FT_start,
                        (FT_start[0] + 500, FT_start[1]),
                        (FT_start[0] + 500, FT_start[1] - params['small_rect_width']),
                        (FT_start[0], FT_start[1] - params['small_rect_width']),
                        FT_start
                    ]
                    dummy_FT_quadrants_outer = [2,1,4,3]
                    dummy_FT_clockwises_outer = [True, True, True, True]
                    filleted_points_outer = []
                    for point, quadrant, clockwise in zip(dummy_FT_pts_outer, dummy_FT_quadrants_outer, dummy_FT_clockwises_outer):
                        radius = params['inner_radius'] if not clockwise else params['outer_radius']
                        filleted_points_outer.extend(cornerRound(point, quadrant, radius, clockwise=clockwise))
                    self.add(SolidPline((i*arrayspacing_x,j*arrayspacing_y), points=filleted_points_outer, layer=layer))
                    

                    dummy_FT_pts_inner = [
                        (FT_start[0] + 10, FT_start[1] - 10),
                        (FT_start[0] + 490, FT_start[1] - 10),
                        (FT_start[0] + 490, FT_start[1] - params['small_rect_width'] + 10),
                        (FT_start[0] + 10, FT_start[1] - params['small_rect_width'] + 10),
                        (FT_start[0] + 10, FT_start[1] - 10)
                    ]
                    dummy_FT_quadrants_inner = dummy_FT_quadrants_outer
                    dummy_FT_clockwises_inner = dummy_FT_clockwises_outer
                    filleted_points_inner = []
                    for point, quadrant, clockwise in zip(dummy_FT_pts_inner, dummy_FT_quadrants_inner, dummy_FT_clockwises_inner):
                        radius = params['outer_radius'] if not clockwise else params['inner_radius']
                        filleted_points_inner.extend(cornerRound(point, quadrant, radius, clockwise=clockwise))
                    self.add(SolidPline((i*arrayspacing_x,j*arrayspacing_y), points=filleted_points_inner, layer=layer))



                separation = 400
                padh = 500
                padw = 500
                padradius=20
                shunt_length = separation+2*padh

                # load doses
                pads_shunt_dose = Doses['pads_shunt_dose'][i,j]
                leads_contactpads_dose = Doses['leads_contactpads_dose'][i,j]
                bridge_dose = Doses['bridge_dose'][i,j]
                bigfinger_dose = Doses['bigfinger_dose'][i,j]
                smallfinger_dose = Doses['smallfinger_dose'][i,j]
                undercut_dose = Doses['undercut_dose'][i,j]
                shift_dose = Doses['shift_dose'][i,j]
                label_dose = Doses['label_dose'][i,j]
                

                # JJ params
                # load JJ params
                finger_width = JJparams['finger_width'][i,j]
                finger_length = JJparams['finger_length'][i,j]
                bridge_length = JJparams['bridge_length'][i,j]
                bigfinger_length = JJparams['bigfinger_length'][i,j]
                # bridge_offset = JJparams['bridge_offset'][i,j]
                # bigfinger_offset = JJparams['bigfinger_offset'][i,j]
                bridge_width = JJparams['bridge_width'][i,j]
                bigfinger_width = JJparams['bigfinger_width'][i,j]
                
                gap = bridge_width
                fingerL = finger_length
                bigfingerW = bigfinger_width
                smallfingerW = finger_width
                bridgeW = bridge_width
                bridgeL = bridge_length
                undercut = 0.2
                total_JJ_length = 2*fingerL + bridgeL
                # print('total_JJ_length = ', total_JJ_length)

                # self.wafer.defaultLayer='PADS_'+str(pads_shunt_dose)
                # print(self.wafer.defaultLayer)
                # # Draw a transmon 3D with shunt
                Transmon3DWithShunt(self, (params['startpoint'][0]+i*arrayspacing_x, params['startpoint'][1]+j*arrayspacing_y+100) , 
                                    padw=padw, padh=padh, leadw=100, leadh=2000, padradius=padradius, separation=separation, shunt=True, 
                                    shunt_width=10, shunt_dist=300, shunt_length=shunt_length, shunt_side='left', flipped=True, layer='PADS_'+str(pads_shunt_dose))
    
                # draw custom leads for this design.
                # lead and contacat pad params
                wedgetowedgeL = 2*81.5 + 17
                wedgetoJJL = wedgetowedgeL/2 - total_JJ_length/2
                leadW = 1
                contactL = 10
                wedgeL = 10
                contactW = 20
                overlap_between_leads_and_bigpads = 20
                leadLpadtocontact = separation/2 - wedgetowedgeL/2 - wedgeL - contactL + overlap_between_leads_and_bigpads
                # print('leadLpadtocontact = ', leadLpadtocontact, ' wedgetoJJL = ', wedgetoJJL)
                
                
                # top one. 
                leads_for_tmon_dosearray_custom(self, m.Structure(self, start=(params['startpoint'][0]+i*arrayspacing_x+padw/2, params['startpoint'][1]+j*arrayspacing_y+100+overlap_between_leads_and_bigpads),direction=0),
                              start=(0,0), 
                              contactpads=False, layer='LEAD_'+str(leads_contactpads_dose), toporbottom='top',
                              leadLpadtocontact=leadLpadtocontact, leadLcontacttoJJ=wedgetoJJL, leadW=leadW, startshifty=1.5, contactW=contactW, contactL=contactL, wedgeL=wedgeL
                )
                
                # bottom one
                leads_for_tmon_dosearray_custom(self, m.Structure(self, start=(params['startpoint'][0]+i*arrayspacing_x+padw/2, params['startpoint'][1]+j*arrayspacing_y-(separation-100+overlap_between_leads_and_bigpads)),direction=0),
                              start=(0,0),  
                              contactpads=False, layer='LEAD_'+str(leads_contactpads_dose), toporbottom='bottom',
                              leadLpadtocontact=leadLpadtocontact, leadLcontacttoJJ=wedgetoJJL, leadW=leadW, startshifty=1.5, contactW=contactW, contactL=contactL, wedgeL=wedgeL              
                )


             
                # Add a small JJ
                smallJJ(self, m.Structure(self, start=(params['startpoint'][0]+i*arrayspacing_x+padw/2+leadW/2, params['startpoint'][1]+j*arrayspacing_y-wedgeL-wedgetoJJL-total_JJ_length/2),direction=90), 
                        smallfingerlayer='SMALLFINGER_'+str(smallfinger_dose), bigfingerlayer='BIGFINGER_'+str(bigfinger_dose), bridgelayer='BRIDGE_'+str(bridge_dose), 
                        Undercutlayer='UNDERCUT_'+str(undercut_dose), shiftlayer = 'SHIFT_'+str(shift_dose),
                        gap=gap, leadW = leadW, fingerL=fingerL, bigfingerW=bigfingerW, smallfingerW=smallfingerW, bridgeW=bridgeW, bridgeL=bridgeL, undercut=undercut)
                
                # labels
                fontsize=30
                xpos = params['startpoint'][0]+i*arrayspacing_x+1.5*padw
                ypos = params['startpoint'][1]+j*arrayspacing_y+padh

                

                if transmon_number_label == True:
                    self.add_chip_label('JJ ('+str(i)+','+str(j)+')', layer='LABEL', 
                                        position=(xpos, params['startpoint'][1]+j*arrayspacing_y+padh), height=fontsize
                                        )
                
                if JJparams_label == True:
                    ypos = ypos-fontsize-15
                    label_list =               [finger_length,  finger_width,   bigfinger_length, bigfinger_width, bridge_length , bridge_width]
                    for ii,label in enumerate(['SmallFingerL', 'SmallFingerW', 'BigfingerL',     'BigFingerW',    'BridgeL',      'BridgeW' ]):
                        self.add_chip_label(label+ ' = ' + str(round(label_list[ii], 3)), 
                                        layer='LABEL', 
                                        position=(xpos, ypos), height=fontsize
                                        )
                        ypos-=fontsize+10
                if Doselabels == True:
                    ypos = params['startpoint'][1]+j*arrayspacing_y-separation
                    label_list = [pads_shunt_dose, leads_contactpads_dose, bridge_dose, bigfinger_dose, smallfinger_dose, undercut_dose, shift_dose, label_dose]
                    for ii, label in enumerate(['pads_shunt_dose', 'leads_contactpads_dose', 'bridge_dose', 'bigfinger_dose', 'smallfinger_dose', 'undercut_dose', 
                                                'shift_dose', 'label_dose']):
                        if label == 'leads_contactpads_dose':
                            xpos+=70
                        self.add_chip_label(label+ ' = ' + str(round(label_list[ii], 3)), 
                                        layer='LABEL', 
                                        position=(xpos, ypos), height=fontsize
                                        )
                        ypos-=fontsize+10
                
                
                

# small JJ from fluxonium.py

# TmonDA = TmonDoseArray(w, arraychip.chipID, w.defaultLayer, jfingerw=1.08)
# TmonDA.add_chip_label(arraychip.chipID, (0, 1000), height=500)


# This will set the default chip for the wafer, filling the chip buffer with this chip


# for 15x15 = 225 dose array:
arraydims = (15,15)
arrayspacing_x = 2050
arrayspacing_y = 2000


base_array = np.ones(arraydims)

finger_width = np.ones(arraydims)
finger_width[:,0:5] = 0.100
finger_width[:,5:10] = 0.140
finger_width[:,10:15] = 0.180
finger_length = 1.5 * np.ones(arraydims)

bigfinger_offset = 0.2 * np.ones(arraydims)
bigfinger_width = finger_width + bigfinger_offset
bigfinger_length = 1.5 * np.ones(arraydims)

bridge_offset = 0.7 * np.ones(arraydims)
bridge_width = finger_width + bridge_offset
bridge_length = 0.250 * np.ones(arraydims)



#things we want to dose separately, for a given element of the transmon array
pads_shunt_dose = np.ones(arraydims)
pads_shunt_dose[0:5] = 600
pads_shunt_dose[5:10] = 650
pads_shunt_dose[10:15] = 700
leads_contactpads_dose = np.ones(arraydims) * 1200

bridge_dose = np.ones(arraydims)
bridge_dose[:,0::5] = 400
bridge_dose[:,1::5] = 450
bridge_dose[:,2::5] = 500
bridge_dose[:,3::5] = 550
bridge_dose[:,4::5] = 600

bigfinger_dose = np.ones(arraydims)
bigfinger_dose[0::5, :] = 800
bigfinger_dose[1::5, :] = 900
bigfinger_dose[2::5, :] = 1000
bigfinger_dose[3::5, :] = 1100
bigfinger_dose[4::5, :] = 1200
smallfinger_dose = np.ones(arraydims)
smallfinger_dose[0::5, :] = 800
smallfinger_dose[1::5, :] = 900
smallfinger_dose[2::5, :] = 1000
smallfinger_dose[3::5, :] = 1100
smallfinger_dose[4::5, :] = 1200

undercut_dose = np.ones(arraydims) * 200
shift_dose = np.ones(arraydims) * 400

label_dose = np.ones(arraydims) * 700



JJparams = {
    'finger_width': finger_width,
    'finger_length': finger_length,
    'bridge_length': bridge_length,
    'bigfinger_length': bigfinger_length,
    'bridge_offset': bridge_offset,
    'bigfinger_offset': bigfinger_offset,
    'bridge_width': bridge_width,
    'bigfinger_width': bigfinger_width
}


Doses = {
    'pads_shunt_dose': pads_shunt_dose,
    'leads_contactpads_dose': leads_contactpads_dose,
    'bridge_dose': bridge_dose,
    'bigfinger_dose': bigfinger_dose,
    'smallfinger_dose': smallfinger_dose,
    'undercut_dose': undercut_dose,
    'shift_dose': shift_dose,
    'label_dose': label_dose
}


w.setDefaultChip(TmonDoseArray(w, 'TmonDA_mainchip', w.defaultLayer, 
                               arraydims=arraydims, arrayspacing_x=arrayspacing_x,arrayspacing_y=arrayspacing_y,transmon_number_label=True, 
                               JJparams = JJparams, JJparams_label = True, Doses = Doses, Doselabels=True,
                               defaults={'r_out': 0, 'r_ins': 0}))

# This goes through the chip buffer and sets each entry to a new chip we define.
# for i in range(0, len(w.chips)):
#     chip = TmonDoseArray(w, 'tmonDA_wafer' + str(i), w.defaultLayer, jfingerw=junc_ws[i])
#     chip.add_chip_label(f'Chip {i}', (40000, 1000),height=500)  # Add label to the chip's block
#     w.setChipBuffer(chip.save(w), i) # Save the chip, which adds the block to the wafer's block list.

# print(w.layernames)

# Now that all chips are saved in the blocks section, write instances of the chips at the right spots on the wafer
w.populate()
w.save()