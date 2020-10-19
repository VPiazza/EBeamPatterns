###################################
####### RUN ON A POWERSHELL ########
###################################
## Path: shell: 'C:\Users\valee\VPScripts\Work\E-BeamPatterns\100 Large Area Patterns' ##
## To run it: python 1cmx1cm.py

## 1. ADJUST POSITION SMALL FIELDS
## 2. DEFINE SLITS IN SMALL FIELDS


# -*- coding: utf-8 -*-
"""
Created on Fri Sept 11 08:51:31 2020

@author: Valerio Piazza

1.4.1 - Added cleaved cross section wires to quickly check the cross section of many wires by cleaving + SEM
1.5 - Added thicker slits in cleaving X-section region and theory cell: see if we can grow single NWs using wider slits
1.5.1 - Changed slit widths to (only) 20/35/50nm and got rid of other slits sizes to simplify
1.5.1 - Got rid of 2um pitch field, as growth was not good here
1.5.1 - Removed branched theory fields to reduce writing time as these were not being used
1.5.1 - Changed all slits to path objects (instead of rectangles) 

2.0.1 - Implementation of gdsCAD for device-specific nanopatterning. Based on the work of Martin Friedl.
2.1.1 - Changing the structure of the chip: from 4-field based to concentric frame.
PROBLEMS:
12.10.2020:     overlapping PAMMs/contacts
                definition small fingers
"""

import sys
module_path = "C:\\Users\\valee\\VPScripts\\Work\\E-BeamPatterns" # Each time you change the path of .py file, you MUST change this path accordingly
sys.path.append(module_path)
print(sys.path)

import itertools
import json
import os.path
from datetime import date
from random import choice as random_choice

import numpy as np

from Patterns.GrowthTheoryCell import make_theory_cell
from Patterns.QuantumPlayground_100_v1 import make_qp
from gdsCAD_py3.core import Cell, Boundary, CellArray, Layout, Path
from gdsCAD_py3.shapes import Box, Rectangle, Label
from gdsCAD_py3.templates100 import Wafer_GridStyle, dashed_line

WAFER_ID = '000058074228'  # CHANGE THIS FOR EACH DIFFERENT WAFER
PATTERN = 'SQ2.0.1'
putOnWafer = True  # Output full wafer or just a single pattern?
HighDensity = False  # High density of triangles?
glbAlignmentMarks = False
tDicingMarks = 10.  # Dicing mark line thickness (um)
pad_size = 200.     # Square contact pad side lenght
finger_width = 50.  # Connecting metal bars
finger_length = 125.
rotAngle = 0.  # Rotation angle of the membranes
wafer_r = 25e3
# waferVer = "100 Membranes SQ1.5.1".format(int(wafer_r / 1000))            #Error debugging(Valerio)
# waferVer = "100 Membranes SQ1.5.1".format()                       >> VP_mod: disabled and chagend (below) <<
waferVer = "100_Membranes".format()



waferLabel = waferVer + '\n' + date.today().strftime("%d%m%Y")
# Layers
l_smBeam = 0
l_lgBeam = 1
l_LGBeam = 10
l_drawing = 100         


# %% Wafer template for MBE growth
class MBE100Wafer(Wafer_GridStyle):
    """
    A 2" wafer divided into square cells
    """

    def __init__(self, name, cells=None):
        Wafer_GridStyle.__init__(self, name=name, cells=cells, block_gap=1200.)

        # The placement of the wafer alignment markers
        am_x = 1.5e4
        am_y = 1.5e4
        self.align_pts = np.array([am_x, am_y])
        self.align_pts = np.vstack((self.align_pts, self.align_pts *
                                    (-1, 1)))  # Reflect about y-axis
        self.align_pts = np.vstack((self.align_pts, self.align_pts *
                                    (1, -1)))  # Reflect about x-axis

        # Size of the chip and placement:
        self.wafer_r = 25e3                                 # 2" wafer
        self.block_size = np.array([10e3, 10e3])            # 1 cm chip                          
        self._place_blocks(radius=self.wafer_r + 5e3)       # top-right corner of the chip array
        self.add_blocks()

        self.add_wafer_outline(layers=l_drawing)                # draws wafer 
        self.add_dashed_dicing_marks(layers=[l_lgBeam])         # Add dicing lines
        self.add_subdicing_marks(200, 5, layers=[l_lgBeam])     # Perpendicular dicing lines
        
        self.add_block_labels(l_lgBeam, unique_ids=True, load_ids=True)     # Chip ID
        self.add_prealignment_markers(layers=[l_lgBeam])                    # Pre-align. Marks
        self.add_contacts(pad_size,finger_width, finger_length, l_LGBeam)
        #self.add_tem_membranes([0.02, 0.035, 0.050], 1000, 1, l_smBeam)     # TEM Sections 
        #self.add_theory_cells()                                             # Growth Section
        #self.add_cleave_xsection_nws()                                      # Cross-section Section                           
        self.add_chip_labels()                                              # Wafer and Pattern name on each chip

        bottom = np.array([0, -self.wafer_r * 0.9])
        self.add_waferLabel(waferLabel, l_drawing, pos=bottom)

    def add_dashed_dicing_marks(self, layers):
        if type(layers) is not list:
            layers = [layers]
        width = 10. / 2
        dashlength = 2000
        r = self.wafer_r
        rng = np.floor(self.wafer_r / self.block_size).astype(int)
        dmarks = Cell('DIC_MRKS')
        for l in layers:
            for x in np.arange(-rng[0], rng[0] + 1) * self.block_size[0]:
                y = np.sqrt(r ** 2 - x ** 2)
                vm = dashed_line([x, y], [x, -y], dashlength, width, layer=l)
                dmarks.add(vm)

            for y in np.arange(-rng[1], rng[1] + 1) * self.block_size[1]:
                x = np.sqrt(r ** 2 - y ** 2)
                hm = dashed_line([x, y], [-x, y], dashlength, width, layer=l)
                dmarks.add(hm)
        self.add(dmarks)

    def add_subdicing_marks(self, length, width, layers):
        if type(layers) is not list:
            layers = [layers]

        for l in layers:
            mark_cell = Cell("SubdicingMark")
            line = Path([[0, 0], [0, length]], width=width, layer=l)
            mark_cell.add(line)

            for block in self.blocks:
                block.add(mark_cell, origin=(self.block_size[0] / 2., 0), rotation=0)
                block.add(mark_cell, origin=(0, self.block_size[1] / 2.), rotation=-90)
                block.add(mark_cell, origin=(self.block_size[0], self.block_size[1] / 2.), rotation=90)
                block.add(mark_cell, origin=(self.block_size[0] / 2., self.block_size[1]), rotation=180)

    def add_block_labels(self, layers, unique_ids=False, save_ids=True, load_ids=True):
        if type(layers) is not list:
            layers = [layers]

        txtSize = 700
        blockids = []

        if not unique_ids:
            for (i, pt) in enumerate(self.block_pts):
                blockids.append(self.blockcols[pt[0]] + self.blockrows[pt[1]])
        else:
            existing_ids = {}
            existing_id_set = set()

            # Load the previously-used IDs from a JSON file
            if load_ids:
                master_db = '../../../ChipIDs_DB.json'
                if os.path.isfile(master_db):
                    with open(master_db, 'r') as f:
                        try:
                            existing_ids = json.load(f)
                            existing_id_set = set([item for sublist in list(existing_ids.values()) for item in sublist])

                            # Check if wafer is in the loaded database
                            if load_ids and WAFER_ID in existing_ids:
                                blockids = existing_ids[WAFER_ID]

                        # If there is a reading error then proceed with a warning
                        except json.decoder.JSONDecodeError:
                            print("Json Error: Couldn't load chip IDs from database!")
                            existing_id_set = set()

            # If the IDs haven't already been set by loading them from the database
            if not blockids:
                # Generate some new IDs, but only use the ones that haven't previously been used
                unique_label_string = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890'
                possible_label_set = set(["".join(x) for x in itertools.product(unique_label_string, repeat=2)])
                possible_label_set = possible_label_set - existing_id_set  # Remove chip lbls from the set of possible lbls

                blockids_set = set()
                while len(blockids_set) < len(self.blocks):
                    blockids_set.add(random_choice(list(possible_label_set)))

                blockids = list(blockids_set)

        # Save the labels to a file
        if save_ids:
            existing_ids.update({WAFER_ID: blockids})
            json_string = json.dumps(existing_ids)
            json_string = json_string.replace("], ", "],\n")  # Make the file a bit easier to read in notepad
            with open(master_db, 'w') as f:
                f.write(json_string)

        # Write the labels to the cells
        for i, block in enumerate(self.blocks):
            blocklabel = Cell('LBL_B_' + blockids[i])
            for l in layers:
                txt = Label(blockids[i], txtSize, layer=l)
                bbox = txt.bounding_box
                offset = (0, 0)
                txt.translate(-np.mean(bbox, 0))  # Center text around origin
                txt.translate(offset)  # Translate it to bottom of wafer
                blocklabel.add(txt)
                block.add(blocklabel, origin=(self.block_size[0] / 2., self.block_size[1] / 2. + 000)) #Shifted by 0.0 mm in y. ???

    def add_prealignment_markers(self, layers, mrkr_size=7):
        if mrkr_size % 2 == 0:  # Number is even, but we need odd numbers
            mrkr_size += 1
        if type(layers) is not list:
            layers = [layers]

        for l in layers:
            rect_size = 10.  # 10 um large PAMM rectangles
            marker_rect = Rectangle([-rect_size / 2., -rect_size / 2.], [rect_size / 2., rect_size / 2.], layer=l) # Coordinates
            marker = Cell('10umMarker')
            marker.add(marker_rect)

            # Make one arm of the PAMM array
            marker_arm = Cell('PAMM_Arm')
            # Define the positions of the markers, they increase in spacing by 1 um each time:
            mrkr_positions = [75 * n + (n - 1) * n // 2 for n in range(1, (mrkr_size - 1) // 2 + 1)]
            for pos in mrkr_positions:
                marker_arm.add(marker, origin=[pos, 0])

            # Build the final PAMM Marker
            pamm_cell = Cell('PAMM_Marker')
            pamm_cell.add(marker)  # Center marker
            pamm_cell.add(marker_arm)  # Right arm
            pamm_cell.add(marker_arm, rotation=180)  # Left arm
            pamm_cell.add(marker_arm, rotation=90)  # Top arm
            pamm_cell.add(marker_arm, rotation=-90)  # Bottom arm
            for pos in mrkr_positions:
                pamm_cell.add(marker_arm, origin=[pos, 0], rotation=90)  # Top arms
                pamm_cell.add(marker_arm, origin=[-pos, 0], rotation=90)
                pamm_cell.add(marker_arm, origin=[pos, 0], rotation=-90)  # Bottom arms
                pamm_cell.add(marker_arm, origin=[-pos, 0], rotation=-90)

            # Make the 4 tick marks that mark the center of the array
            h = 30.
            w = 70.
            tick_mrk = Rectangle([-w / 2., -h / 2.], [w / 2, h / 2.], layer=l)
            tick_mrk_cell = Cell("TickMark")
            tick_mrk_cell.add(tick_mrk)
            pos = mrkr_positions[-1] + 50 + w / 2.
            pamm_cell.add(tick_mrk_cell, origin=[pos, 0])
            pamm_cell.add(tick_mrk_cell, origin=[-pos, 0])
            pamm_cell.add(tick_mrk_cell, origin=[0, pos], rotation=90)
            pamm_cell.add(tick_mrk_cell, origin=[0, -pos], rotation=90)

        center_x, center_y = (5000, 5000)                                           # Positioning onto the chip
        for block in self.blocks:
            block.add(pamm_cell, origin=(center_x + 2000, center_y + 2000))
            block.add(pamm_cell, origin=(center_x - 2000, center_y - 2000))
        
## First attempt to design the contact pads into the outer 1.5 mm frame
    def add_contacts(self, pad_size, finger_width, finger_length, layers):
        contact_frames = [2,4,6,8]
        spacing = pad_size/2
        for frame in contact_frames:
            frame_length = (10-frame)*1000
            n_cont = int((frame_length)/(pad_size+spacing))-1
            

            corner_pos = pad_size/2

            contact_pads = Cell('Contact_Pads')
            pad =  Rectangle((-corner_pos,-corner_pos), (corner_pos,corner_pos), layer=layers)
            pad_cell = Cell('Pad_Cell')
            pad_cell.add(pad)
            finger = Rectangle((-finger_width/2,-finger_length/2), (finger_width/2,finger_length/2), layer=layers)
            finger_cell = Cell('Finger Cell')
            finger_cell.add(finger)

            curr_x = (10000 - ((n_cont-1)*(pad_size + spacing)))/2
            curr_y = (10000 - frame_length)/2
            pad_array = CellArray(pad_cell, n_cont, 1, (pad_size + spacing, pad_size + spacing), origin = (curr_x, curr_y))
            finger_array_tr = CellArray(finger_cell, n_cont, 1, (pad_size + spacing, pad_size + spacing), origin = (curr_x - pad_size/2 + finger_width, curr_y + corner_pos + finger_length/2))
            finger_array_tl = CellArray(finger_cell, n_cont, 1, (pad_size + spacing, pad_size + spacing), origin = (curr_x + pad_size/2 - finger_width, curr_y + corner_pos + finger_length/2))
            finger_array_br = CellArray(finger_cell, n_cont, 1, (pad_size + spacing, pad_size + spacing), origin = (curr_x - pad_size/2 + finger_width, curr_y - corner_pos - finger_length/2))
            finger_array_bl = CellArray(finger_cell, n_cont, 1, (pad_size + spacing, pad_size + spacing), origin = (curr_x + pad_size/2 - finger_width, curr_y - corner_pos - finger_length/2))

            
            contact_pads.add(pad_array)
            if frame < 8:
                contact_pads.add(finger_array_tr)
                contact_pads.add(finger_array_tl)
            if frame > 2:
                contact_pads.add(finger_array_br)
                contact_pads.add(finger_array_bl)


            for block in self.blocks:
                block.add(contact_pads)
                block.add(contact_pads, origin = (10000,0), rotation = 90)
                block.add(contact_pads, origin = (10000,10000), rotation = 180)
                block.add(contact_pads, origin = (0, 10000), rotation = 270)


    def add_tem_membranes(self, widths, length, pitch, layer):
        tem_membranes = Cell('TEM_Membranes')
        n = 4  # Number of NMs in each width group (side-by-side)
        curr_y = 0
        for width in widths:
            membrane = Path([(-length / 2., 0), (length / 2., 0)], width=width, layer=layer)
            membrane_cell = Cell('Membrane_w{:.0f}'.format(width * 1000))
            membrane_cell.add(membrane)
            membrane_array = CellArray(membrane_cell, 1, n, (0, pitch))
            membrane_array_cell = Cell('MembraneArray_w{:.0f}'.format(width * 1000))
            membrane_array_cell.add(membrane_array)
            tem_membranes.add(membrane_array_cell, origin=(0, curr_y))
            curr_y += n * pitch

        n2 = 4  # Number of times to repeat each width group in the overall array             
        tem_membranes2 = Cell('Many_TEM_Membranes')
        tem_membranes2.add(CellArray(tem_membranes, 1, n2, (0, n * len(widths) * pitch)))

        center_x, center_y = (5000, 5000)
        for block in self.blocks:
            block.add(tem_membranes2, origin=(center_x, center_y + 2000))
            block.add(tem_membranes2, origin=(center_x, center_y + 1500), rotation=45)

    def add_theory_cells(self):
        theory_cells = Cell('TheoryCells')
        theory_cells.add(make_theory_cell(wafer_orient='100'), origin=(70, 0))
        # theory_cells.add(make_theory_cell_3br(), origin=(0, 0))
        # theory_cells.add(make_theory_cell_4br(), origin=(400, 0))

        theory_cells.add(make_theory_cell(wafer_orient='100'), origin=(20, -400), rotation=45)
        # theory_cells.add(make_theory_cell_3br(), origin=(-50, -400), rotation=45)
        # theory_cells.add(make_theory_cell_4br(), origin=(370, -400), rotation=45)

        center_x, center_y = (5000, 5000)
        for block in self.blocks:
            block.add(theory_cells, origin=(center_x, center_y - 1700))

    def add_chip_labels(self):
        wafer_lbl = PATTERN + '\n' + WAFER_ID
        text = Label(wafer_lbl, 20., layer=l_lgBeam)
        text.translate(tuple(np.array(-text.bounding_box.mean(0))))  # Center justify label
        chip_lbl_cell = Cell('chip_label')
        chip_lbl_cell.add(text)

        center_x, center_y = (5000, 5000)
        for block in self.blocks:
            block.add(chip_lbl_cell, origin=(center_x, center_y - 500))

    def add_cleave_xsection_nws(self):
        pitches = [0.5, 1., 2., 4.]
        widths = [10., 20., 40., 60., 100., 160., 240.]
        n_membranes = 10
        length = 50
        spacing = 10

        cleave_xsection_cell = Cell("CleaveCrossSection")

        y_offset = 0
        for pitch in pitches:
            for width in widths:
                nm_cell = Cell("P{:.0f}W{:.0f}".format(pitch, width))
                slit = Path([(-length / 2., 0), (length / 2., 0)], width=width / 1000., layer=l_smBeam)
                nm_cell.add(slit)
                nm_cell_array = Cell("P{:.0f}W{:.0f}_Array".format(pitch, width))
                tmp = CellArray(nm_cell, 1.0, n_membranes, [0, pitch])
                nm_cell_array.add(tmp)
                cleave_xsection_cell.add(nm_cell_array, origin=(0, y_offset + pitch * n_membranes))
                y_offset += pitch * n_membranes + spacing

                text = Label("P{:.1f}W{:.0f}".format(pitch, width), 1.0, layer=l_smBeam)
                text.translate(tuple(np.array(-text.bounding_box.mean(0))))  # Center justify label
                txt_cell = Cell("lbl_P{:.1f}W{:.0f}".format(pitch, width))
                txt_cell.add(text)
                cleave_xsection_cell.add(txt_cell, origin=(length * 0.75, y_offset - 8.0))
                cleave_xsection_cell.add(txt_cell, origin=(-length * 0.75, y_offset - 8.0))

            y_offset += spacing * 3

        center_x, center_y = (5000, 5000)
        for block in self.blocks:
            block.add(cleave_xsection_cell, origin=(center_x + 1150, center_y - 463))
            # block.add(cleave_xsection_cell, origin=(center_x - 350, center_y + 350), rotation=45.)    # >> VP_mod: disabled <<
            block.add(cleave_xsection_cell, origin=(center_x + 463, center_y - 1150), rotation=90.)   # >> VP_mod: disabled<<


class Frame(Cell):
    """
    Make a frame for writing to with ebeam lithography
    Params:
    -name of the frame, just like when naming a cell
    -size: the size of the frame as an array [xsize,ysize]
    """

    def __init__(self, name, size, border_layers):
        if not (type(border_layers) == list):
            border_layers = [border_layers]
        Cell.__init__(self, name)
        self.size_x, self.size_y = size
        # Create the border of the cell
        for l in border_layers:
            self.border = Box(
                (-self.size_x / 2., -self.size_y / 2.),
                (self.size_x / 2., self.size_y / 2.),
                1,
                layer=l)
            self.add(self.border)  # Add border to the frame

        self.align_markers = None

    def make_align_markers(self, t, w, position, layers, joy_markers=False, camps_markers=False):
        if not (type(layers) == list):
            layers = [layers]
        top_mk_cell = Cell('AlignmentMark')
        for l in layers:
            if not joy_markers:
                am0 = Rectangle((-w / 2., -w / 2.), (w / 2., w / 2.), layer=l)
                rect_mk_cell = Cell("RectMarker")
                rect_mk_cell.add(am0)
                top_mk_cell.add(rect_mk_cell)
            elif joy_markers:
                crosspts = [(0, 0), (w / 2., 0), (w / 2., t), (t, t), (t, w / 2), (0, w / 2), (0, 0)]
                crosspts.extend(tuple(map(tuple, (-np.array(crosspts)).tolist())))
                am0 = Boundary(crosspts, layer=l)  # Create gdsCAD shape
                joy_mk_cell = Cell("JOYMarker")
                joy_mk_cell.add(am0)
                top_mk_cell.add(joy_mk_cell)

            if camps_markers:
                emw = 20.  # 20 um e-beam marker width
                camps_mk = Rectangle((-emw / 2., -emw / 2.), (emw / 2., emw / 2.), layer=l)
                camps_mk_cell = Cell("CAMPSMarker")
                camps_mk_cell.add(camps_mk)
                top_mk_cell.add(camps_mk_cell, origin=[100., 100.])
                top_mk_cell.add(camps_mk_cell, origin=[100., -100.])
                top_mk_cell.add(camps_mk_cell, origin=[-100., 100.])
                top_mk_cell.add(camps_mk_cell, origin=[-100., -100.])

            self.align_markers = Cell("AlignMarkers")
            self.align_markers.add(top_mk_cell, origin=np.array(position) * np.array([1, -1]))
            self.align_markers.add(top_mk_cell, origin=np.array(position) * np.array([-1, -1]))
            self.align_markers.add(top_mk_cell, origin=np.array(position) * np.array([1, 1]))
            self.align_markers.add(top_mk_cell, origin=np.array(position) * np.array([-1, 1]))
            self.add(self.align_markers)

### Define slits for device analysis ##
                        
    def make_slit_patterns(self, sflabels, _pitches, spacing, _widths, _lengths, rot_angle,
                        array_height, array_width, array_spacing, layers):
        if not (type(layers) == list):
            layers = [layers]
        if not (type(_pitches) == list):
            _pitches = [_pitches]
        if not (type(_lengths) == list):
            _lengths = [_lengths]
        if not (type(_widths) == list):
            _widths = [_widths]
        manyslits = i = j = None
        for l in layers:
            i = -1
            j = -1
            manyslits = Cell("SlitArray")
            pitch = _pitches[0]
            for length in _lengths:
                j += 1
                i = -1

                for width in _widths:
                    # for pitch in pitches:
                    i += 1
                    if i % 3 == 0:
                        j += 1  # Move to array to next line
                        i = 0  # Restart at left

                    nx = int(array_width / (length + spacing))
                    ny = int(array_height / pitch)
                    # Define the slits
                    slit = Cell("Slit_w{:.0f}".format(width * 1000))
                    slit_path = Path([(-length / 2., 0), (length / 2., 0)], width=width, layer=l)
                    slit.add(slit_path)
                    slits = CellArray(slit, nx, ny, (length + spacing, pitch))
                    slits.translate((-(nx - 1) * (length + spacing) / 2., -(ny - 1) * pitch / 2.))
                    slit_array = Cell("SlitArray_w{:.0f}".format(width * 1000))
                    slit_array.add(slits)
                    text = Label('w/p/l\n%i/%i/%i' % (width * 1000, pitch, length), 5, layer=l)
                    lbl_vertical_offset = 1.35
                    if j % 2 == 0:
                        text.translate(
                            tuple(np.array(-text.bounding_box.mean(0)) + np.array((
                                0, -array_height / lbl_vertical_offset))))  # Center justify label
                    else:
                        text.translate(
                            tuple(np.array(-text.bounding_box.mean(0)) + np.array((
                                0, array_height / lbl_vertical_offset))))  # Center justify label
                    slit_array.add(text)
                    manyslits.add(slit_array,
                                  origin=((array_width + array_spacing) * i, (
                                          array_height + 2. * array_spacing) * j - array_spacing / 2.))
                    
        specific_label = Label(sflabels, 20, layer=l)
        specific_label.translate((-lbl_vertical_offset*smMarkerPosition, -lbl_vertical_offset*smMarkerPosition))  # Center Small Field
        slit_array.add(specific_label)

        # This is an ugly hack to center rotated slits, should fix this properly...
        if rot_angle == 45:  # TODO: fix this ugly thing
            hacky_offset_x = 200
            hacky_offset_y = -25
        elif rot_angle == 90:
            hacky_offset_x =  356
            hacky_offset_y =  96.5
        elif rot_angle == 180:              
            hacky_offset_x = 260
            hacky_offset_y = 452
        elif rot_angle == 270 or rot_angle == -90:             
            hacky_offset_x = -96.5
            hacky_offset_y = 356
        else:
            hacky_offset_x = 0
            hacky_offset_y = 0

        self.add(manyslits, origin=(-(i * (array_width + array_spacing)) / 2 + hacky_offset_x,
                                    -(j + 1.5) * (array_height + array_spacing) / 2 + hacky_offset_y), 
                                    rotation=rot_angle)






# %%

# Define parameters that we will use for the features inside the small fields
widths = [0.020, 0.035, 0.050, 0.020, 0.035, 0.050]
pitch = 3.0
length = 15.

smFrameSize = 100
smFrameSpacing = 100  # Spacing between the small fields

smMarkerPosition = 50.

# Create the Cell, Frames and Small Fields
topCell = Cell("TopCell")

# Center of the Chip
cent_x = -5000
cent_y = -5000

sm_frames = [3,5,7]                                     # Defines the position of the frames

for frame in sm_frames: 
    # Create the pattern that we want to write
    lgField_width = frame*1000
    lgField = Frame("LargeField", (lgField_width, lgField_width), [])  # Create the large write field
    lgField.make_align_markers(10., 200., (lgField_width/2, lgField_width/2), l_lgBeam, joy_markers=True, camps_markers=True)


    topCell.add(lgField, origin=(cent_x, cent_y))

    # Add everything together to a top cell
    
    
    sf_num = int((lgField_width-1000) / (smFrameSize+smFrameSpacing))+1
    
    if sf_num % 2 != 0:
        sf_num += 1

    sm_origin = ((sf_num*smFrameSize)+((sf_num-3)*smFrameSpacing))/2 


    for n in range(0, sf_num):
        smField = Frame("SmallField1", (smFrameSize, smFrameSize), [])
        smField.make_align_markers(2., 20., (smMarkerPosition, smMarkerPosition), l_lgBeam, joy_markers=True)
        #smField.make_slit_patterns("Single Slits", pitch, slitColumnSpacing, widths, length, 0, 100, 100, 30, l_smBeam)

        
        sf_pos = sm_origin - n*(smFrameSize + smFrameSpacing)+0.5*smFrameSize
        distance_from_center = lgField_width/2
        topCell.add(smField,  origin= (cent_x - sf_pos, cent_y - distance_from_center), rotation = 0)
        topCell.add(smField,  origin= (cent_x + distance_from_center, cent_y - sf_pos), rotation = 90)
        topCell.add(smField,  origin= (cent_x + sf_pos, cent_y + distance_from_center), rotation = 180)
        topCell.add(smField,  origin= (cent_x  - distance_from_center, cent_y + sf_pos), rotation = 270)


    # topCell.spacing = np.array([4000., 4000.])        #Create an array 2x2

    

# %%Create the layout and output GDS file
filestring = str(waferVer) + '_' + WAFER_ID + '_' + date.today().strftime("%d%m%Y")

layout = Layout('LIBRARY')
if putOnWafer:  # Fit as many patterns on a 2inch wafer as possible
    wafer = MBE100Wafer('MembranesWafer', cells=[topCell])
    layout.add(wafer)
# layout.show()
else:  # Only output a single copy of the pattern (not on a wafer)
    layout.add(topCell)
    layout.show()
filename = filestring.replace(' ', '_') + '_wafer' + '.gds'
layout.save(filename)

# Output single chip for doing aligned jobs
cell_layout = Layout('LIBRARY')
cell_layout.add(wafer.blocks[0])
cell_layout.save(filestring.replace(' ', '_') + '_chip' + '.gds')

# Output topCell for doing aligned jobs
layout_field = Layout('LIBRARY')
layout_field.add(topCell)
layout_field.save(filestring.replace(' ', '_') + '_smallFields.gds')
