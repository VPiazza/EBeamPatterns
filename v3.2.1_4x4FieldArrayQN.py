###################################
####### RUN ON A POWERSHELL ########
###################################
## Path: shell: 'C:\Users\valee\VPScripts\Work\E-BeamPatterns\100 Large Area Patterns' ##

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

3.0.1 - 4x4 Fields containing 16 small fields each. Based on the work of Martin Friedl.
3.1.1 - 4x4 Fields containing 16 small fields each. Based on the work of Martin Friedl.
        The script is optimized for photodetection. Written in collaboration with Nicholas P. Morgan.
3.2.1 - 4x4 Fields containing 16 small fields each. Based on the work of Martin Friedl.
        The script is optimized for quantum network. Written in collaboration with Didem Dede.

        TO CHECK AT THE END: Modification of the large field and of the small field should modify consistently everything!!
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
from gdsCAD_py3.shapes import Box, Rectangle, Label, RegPolygon
from gdsCAD_py3.templates100 import Wafer_GridStyle, dashed_line

WAFER_ID = '000011111111'  # CHANGE THIS FOR EACH DIFFERENT WAFER
PATTERN = 'SQ3.0.1'
#putOnWafer = True  # Output full wafer or just a single pattern?
HighDensity = False  # High density of triangles?
glbAlignmentMarks = False
tDicingMarks = 10.          # Dicing mark line thickness (um)
lgField_size = 1500.        # Area for parameters variation
lgField_spacing = 2000.
lgField_num = int((10000-lgField_size)/lgField_spacing)     # ACTUAL: 4x4 Lg Fields
pad_size = 200.     # Square contact pad side lenght
smField_size = 100. # Area of each test 
sm_spacing = pad_size + smField_size                        # ACTUAL: 3x3 Sm Fields
smField_num = int((lgField_size-200)/sm_spacing-1)
rotAngle = 0.  # Rotation angle of the membranes
wafer_r = 25e3
waferVer = "100_Membranes".format()

waferLabel = waferVer + '\n' + date.today().strftime("%d%m%Y")
# Layers
l_smBeam = 0        # 2nd job small ebeam 
l_lgBeam = 1        # 1st job large ebeam 
l_PHBeam = 10       # 3rd job large ebeam 
l_drawing = 100         

### Checking geometrical parameters:
precheck= False
intro = "You are starting to write your layout."
wafer_type = "**\nAre writing on a 2\" wafer? [y/n]  "
block_type = "**\nDo you want to write 1cmx1cm chip? [y/n]  "
chip_type = "**\nDo you want to write " + str(lgField_num) + " lines of Large Fields? [y/n]  "
LgField_type = "**\nDo you want to write " + str(smField_num) + " lines of Small Fields in each Large field? [y/n]  "

positive_answer = ""
lets_go = "**\nLet's Go!\n**"
exit_sentence = "**\nThen something is wrong......"

if precheck:
    print(intro)
    if input(wafer_type) == positive_answer:
        if input(block_type) == positive_answer:
            if input(chip_type) == positive_answer:
                if input(LgField_type) == positive_answer:
                    print(lets_go)
                else:
                    print(exit_sentence)
                    quit()
            else:
                print(exit_sentence)
                quit()
        else:
            print(exit_sentence)
            quit()
    else:
        print(exit_sentence)
        quit()

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
        self._place_blocks(radius=self.wafer_r + 5e3)       # distance of the center of the chip with respect to the center of the wafer
        self.add_blocks()

        self.add_wafer_outline(layers=l_drawing)                # draws wafer 
        self.add_dashed_dicing_marks(layers=[l_lgBeam])         # Add dicing lines
        self.add_subdicing_marks(200, 5, layers=[l_lgBeam])     # Perpendicular dicing lines
        
        self.add_block_labels(l_lgBeam, unique_ids=True, load_ids=True)     # Chip ID
        self.add_prealignment_markers(layers=[l_lgBeam])                    # Pre-align. Marks                         
        self.add_chip_labels()                                              # Wafer and Pattern name on each chip
        self.add_contacts(layers= l_PHBeam)

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

        txtSize = 200
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

    def add_prealignment_markers(self, layers, mrkr_size=7):    #TO BE ADJUSTED 
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
            mrkr_positions = [50 * n + (n - 1) * n // 2 for n in range(1, (mrkr_size - 1) // 2 + 1)]
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
            h = 15.
            w = 50.
            tick_mrk = Rectangle([-w / 2., -h / 2.], [w / 2, h / 2.], layer=l)
            tick_mrk_cell = Cell("TickMark")
            tick_mrk_cell.add(tick_mrk)
            pos = mrkr_positions[-1] + 75 + w / 2.
            pamm_cell.add(tick_mrk_cell, origin=[pos, 0])
            pamm_cell.add(tick_mrk_cell, origin=[-pos, 0])
            pamm_cell.add(tick_mrk_cell, origin=[0, pos], rotation=90)
            pamm_cell.add(tick_mrk_cell, origin=[0, -pos], rotation=90)


        center_x, center_y = (5000, 5000)                                           # Positioning onto the chip
        for block in self.blocks:
            block.add(pamm_cell, origin=(center_x + 2000, center_y))
            block.add(pamm_cell, origin=(center_x - 2000, center_y))
            block.add(pamm_cell, origin=(center_x, center_y + 2000))
            # block.add(pamm_cell, origin=(center_x, center_y - 2000))

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
        wafer_lbl = PATTERN + "\n" + WAFER_ID
        text = Label(wafer_lbl, 20., layer=l_lgBeam)
        text.translate(tuple(np.array(-text.bounding_box.mean(0))))  # Center justify label
        chip_lbl_cell = Cell('chip_label')
        chip_lbl_cell.add(text)

        center_x, center_y = (5000, 5000)
        for block in self.blocks:
            block.add(chip_lbl_cell, origin=(center_x, center_y - 200))
            block.add(chip_lbl_cell, origin=(center_x, center_y - 4500))
            block.add(chip_lbl_cell, origin=(center_x + 4500, center_y), rotation= 90)
            block.add(chip_lbl_cell, origin=(center_x, center_y + 4500), rotation= 180)
            block.add(chip_lbl_cell, origin=(center_x - 4500, center_y), rotation= 270)

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

    ## Design of the contact pads into each Large Field
    def add_contacts(self, layers):
        corner_pos = pad_size/2
        finger_width = 20.  
        finger_length = 80.
        n_cont = smField_num + 1

        contact_pads = Cell('Contact_Pads')
        pad =  Rectangle((-corner_pos,-corner_pos), (corner_pos,corner_pos), layer=layers)
        pad_cell = Cell('Pad_Cell')
        pad_cell.add(pad)
        finger = Rectangle((-finger_width/2,-finger_length/2), (finger_width/2,finger_length/2), layer=layers)
        finger_cell = Cell('Finger Cell')
        finger_cell.add(finger)
        n_finger = n_cont - 1

        pad_array = CellArray(pad_cell, n_cont, n_cont, (sm_spacing, sm_spacing), origin = (0, 0))
        finger_array1 = CellArray(finger_cell, n_finger, n_finger, (sm_spacing, sm_spacing), origin=(corner_pos - finger_width, +corner_pos + finger_length/2))
        finger_array2 = CellArray(finger_cell, n_finger, n_finger, (sm_spacing, sm_spacing), origin=(sm_spacing -corner_pos + finger_width, sm_spacing -corner_pos - finger_length/2))
        finger_array3 = CellArray(finger_cell, n_finger, n_finger, (sm_spacing, sm_spacing), rotation = 90, origin=((n_cont-1)*sm_spacing - corner_pos - finger_length/2, corner_pos - finger_width))
        finger_array4 = CellArray(finger_cell, n_finger, n_finger, (sm_spacing, sm_spacing), rotation = 90, origin=((n_cont-2)*sm_spacing + corner_pos + finger_length/2, sm_spacing -corner_pos + finger_width))
       
        contact_pads.add(pad_array)
        contact_pads.add(finger_array1)
        contact_pads.add(finger_array2)
        contact_pads.add(finger_array3)
        contact_pads.add(finger_array4)

        center = -0.5*((n_cont-1)*smField_size + (n_cont-1)*pad_size)
        for block in self.blocks:
            for n in range (0, lgField_num):
                for i in range (0, lgField_num):
                    block.add(contact_pads, origin = (center + (n+1)*lgField_spacing, center+ (i+1)*lgField_spacing))

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
    def make_slits(self, length, width, nslit, pitch, rot_angle, layers):
        """
        Define a single slit or a slit array with a given length, width and pitch
        """
        slitField = Cell("slitField")

        slit = Cell("Single Slit")
        slit_path = Path([(-length / 2., 0), (length / 2., 0)], width = width, layer = layers)
        slit.add(slit_path)
        
        if nslit == 1:
            slitField.add(slit, origin=(0,0), rotation=rot_angle)
        elif nslit > 1:
            slits = CellArray(slit, 1, nslit, (0,pitch))
            slits.translate((0, -(nslit-1) * pitch / 2.))
            slit_array = Cell("Multiple Slit")
            slit_array.add(slits)
            slitField.add(slit_array, origin=(0,0), rotation=rot_angle)
        else:
            print("Error in the number of slits. Check the internal code"*50)
            quit()
        self.add(slitField)
 
    def make_finger_contacts(self, slit_length, length, nslit, pitch, rot_angle, layers):
        global margin

        margin = 2.5
        cont_to_cent = 60.
        global fing_width
        fing_width = 1.

        global rad_angle
        rad_angle = rot_angle/180*np.pi

        fing_ext_length = cont_to_cent - (slit_length/2-margin)*np.cos(rad_angle) 
        fing_ext_hook = 30. + nslit/2 * pitch + np.sin(rad_angle)*slit_length/2
        fing_int_length = cont_to_cent + nslit/2 * pitch + np.sin(rot_angle)*length + margin
        cont_conn_length = 2*margin

        global fake_slit_length
        fake_slit_length = cont_to_cent-fing_ext_length - length - margin - fing_width/2
        
        contact = Cell(" FingerContact")

        triangle = RegPolygon((0,0), cont_conn_length, 3, layer = layers)
        conn_triangle = Cell("ConnTriangle")
        conn_triangle.add(triangle)
        
        finger_rect = Rectangle((-fing_ext_length / 2., -fing_width/2), (fing_ext_length / 2., fing_width/2), layer = layers)
        finger_ext = Cell("Finger")
        finger_ext.add(finger_rect)
        
        hook = Rectangle(( -fing_ext_hook / 2., -fing_width/2), (fing_ext_hook / 2., fing_width/2), layer = layers)
        hook_finger = Cell("ContactFinger") 
        hook_finger.add(hook)

        finger_middle = Rectangle((-fing_int_length / 2., -fing_width/2), (fing_int_length / 2., fing_width/2), layer = layers)
        finger_int = Cell("Middle Finger")
        finger_int.add(finger_middle)

        # Relative coordinates for horizontal elements (origin = center of the triangle)
        h_finger_to_triangle_x = fing_ext_length/2
        h_hook_to_triangle_x = fing_ext_length - fing_width/2
        h_hook_to_triangle_y = -fing_ext_hook/2 + fing_width/2

        # Relative coordinates for vertical elements (origin = center of the triangle)
        v_finger_to_triangle_x = fing_int_length/2

        cont_ext = Cell("ExternalContact")
        cont_ext.add(conn_triangle)
        cont_ext.add(finger_ext, origin = (h_finger_to_triangle_x, 0))
        cont_ext.add(hook_finger, origin =(h_hook_to_triangle_x, h_hook_to_triangle_y), rotation = 90)

        cont_int = Cell("InternalContact")
        cont_int.add(conn_triangle)
        cont_int.add(finger_int, origin = (v_finger_to_triangle_x, 0))

        # Coordinates Horizontal and Vertical Contacts into the small field
        hor_x = -cont_to_cent
        hor_y = cont_to_cent/2-margin
        vert_x = (-length - fing_width/2)*np.cos(rad_angle) - (1-np.cos(rad_angle))*(fing_width/2)
        vert_y = cont_to_cent

        contact.add(cont_ext, origin = (hor_x, hor_y))
        contact.add(cont_int, origin = (vert_x, vert_y), rotation = -90)       

        self.add(contact)
        self.add(contact, rotation = 180)
        return 

    def make_slits_reservoir(self, nslit, pitch, contact_distance, layers): # 5 additional slits as material reservoir
        res_slit = 5

        gap = contact_distance + 2. + 2.  
        #res_length = (length - gap - 2.5*margin)/2 
        res_length = fake_slit_length
        res_width = width_std
        res_pitch = pitch

        resField = Cell("resField")

        reservoir = Cell("Single Reservoir")
        res_path = Path([(-res_length / 2., 0), (res_length / 2., 0)], width = res_width, layer = layers)
        reservoir.add(res_path)

        x_spac = (res_length + gap)/np.cos(rad_angle) 
        y_spac = res_pitch

        reservoirs= CellArray(reservoir, 2, res_slit, spacing = (x_spac, y_spac))
        x_transl = -(res_length + gap)/(2*np.cos(rad_angle)) + (margin)*np.sin(rad_angle)
        reservoirs.translate((x_transl,0))
        res_array = Cell("Multiple Slit")
        res_array.add(reservoirs)
        resField.add(res_array, origin=(0,0), rotation=rot_angle)

        if contact_distance > margin:
            add_slit = Cell("Additional Reservoir")
            add_res_path = Path([(-(contact_distance - margin) / 2., 0), ((contact_distance - margin) / 2., 0)], width = res_width, layer = layers)
            add_slit.add(add_res_path)

            add_reservoir = CellArray(add_slit, 1, res_slit, spacing = (0, res_pitch))
            add_reservoir.translate((0,0))
            add_res_array = Cell("Additional Multiple Slit")
            add_res_array.add(add_reservoir)
            resField.add(add_res_array, origin=(0,0), rotation=rot_angle)



        self.add(resField, origin= (0,(nslit+1) * pitch/2 )/np.cos(rad_angle))
        self.add(resField, origin= (0, -((nslit+1+(2*(res_slit-1))) * pitch/2 )/np.cos(rad_angle)))

# %%Create the pattern that we want to write

# Define parameters that we will use for the slits: 
length = 40.        # Length of the slit. TBN that the lenght is varied by the position of the contact.

# Standard values (used when these params are not varying)
width_std = 0.05
pitch_std = 1.
lenght_std = 10.


# Since the 3x3 arrays of small field, the array MUST contain 3 values.
widths_sm = [0.02, 0.05, 0.10]
widths_lg = [0.15, 0.20, 0.25]                  # 6 widths
lenghts_sm = [0.25, 0.5, 1.] 
lengths_lg = [2., 5., 10.]                      # 6 lenghts
pitches = [0.5, 1., 2.]                         # 3 pitches
num_slits_sm = [2, 5, 10] 
num_slits_lg = [15, 20, 25]                     # 6 multiple slits   25 slits is the maximum for the chosen geometry

topCell = Cell("TopCell")
sm_writer = False
lg_label = ""


# Crate large field array following the geometry set at the beginning. 
for lg_row in range(0, lgField_num):
    for lg_col in range(0, lgField_num):
        
        # Single slit, width/length dependence (orientation 0°).   lgField: 1x1,1x2,2x1,2x2 + 3x3,3x4,4x3,4x4
        if (lg_row+1 == 1 and lg_col+1 == 1) or (lg_row+1 == 2 and lg_col+1 == 1) or (lg_row+1 == 1 and lg_col+1 == 2) or(lg_row+1 == 2 and lg_col+1 == 2) or  (lg_row+1 == 3 and lg_col+1 == 3) or (lg_row+1 == 4 and lg_col+1 == 3) or (lg_row+1 == 3 and lg_col+1 == 4) or (lg_row+1 == 4 and lg_col+1 == 4):
            sm_writer = True
            lg_label = "Single Slit - 0deg  \nPitch = " + str(pitch_std) + "um"
            num_slit = [1,1,1]
            rot_angle = 0
            _pitch = [pitch_std,pitch_std,pitch_std]
            
            if lg_col % 2 == 0:
                _width = widths_sm
            elif lg_col % 2 == 1:
                 _width = widths_lg
            
            if lg_row % 2 == 0:
                _length = lenghts_sm
            elif lg_row % 2 == 1:
                _length = lengths_lg
       
        # Single slit, width/length dependence (orientation 45°):                   1x3,1x4,2x3,2x4 
        elif (lg_col+1 == 3 and lg_row+1 == 1) or (lg_col+1 == 3 and lg_row+1 == 2) or (lg_col+1 == 4 and lg_row+1 == 1) or (lg_col+1 == 4 and lg_row+1 == 2):
            sm_writer = True
            lg_label = "Single Slit - 45deg\nPitch = " + str(pitch_std) + "um"
            num_slit = [1,1,1]
            rot_angle = 45
            _pitch = [pitch_std,pitch_std,pitch_std]
            
            if lg_col+1 ==3:
                _width = widths_sm
            elif lg_col+1 ==4:
                 _width = widths_lg
            
            if lg_row+1 == 1:
                _length = lenghts_sm
            elif lg_row+1 == 2:
                _length = lengths_lg
                
        # Multiple slits, pitch dependence, width/length fixed (orientation 0°):    3x1,3x2,4x1,4x2
        elif (lg_col+1 == 1 and lg_row+1 == 3) or (lg_col+1 == 1 and lg_row+1 == 4) or (lg_col+1 == 2 and lg_row+1 == 3) or (lg_col+1 == 2 and lg_row+1 == 4):
            sm_writer = True
            _width = [width_std , width_std , width_std ]
            lg_label = "Multiple Slit - 0deg\nLength = " + str(lenght_std) + "um\nWidth = " + str(width_std) + "um"
            rot_angle = 0
            _length = [lenght_std, lenght_std, lenght_std] 

            _pitch = pitches
            
            if lg_col+1 == 1:
                num_slit = num_slits_sm
            elif lg_col+1 == 2:
                num_slit = num_slits_lg

        else: 
            lg_label = "ERROR"
            sm_writer = False

        #Coordinates of the large field and Marker positioning geometry
        lg_orig_x, lg_orig_y = -4*lgField_spacing + lg_col*lgField_spacing, -lgField_spacing - lg_row*lgField_spacing
        lgMark_margin = 50.
        lgMark_position = lgField_size/2 - lgMark_margin

        #Large Field 
        #lgField_name = "LargeField_"+str(lg_row+1)+"x"+str(lg_col+1)
        lgField = Frame("LargeField", (lgField_size/2, lgField_size/2), [])  # Create the large write field
        lgField.make_align_markers(10., 200., (lgMark_position, lgMark_position), l_lgBeam, joy_markers=True, camps_markers=True)
        lgField_label = Label(lg_label, 20., position = (lg_orig_x-150,lg_orig_y + lgField_size/2), layer=l_lgBeam)

        # Create the smaller write field aligned with the large field.
        if sm_writer:                               # the sm_writer will avoid errors if there are inconsistencies when changing the large and small field geometry
            for sm_row in range(0, smField_num):
                for sm_col in range(0, smField_num):
                    
                    #Coordinates of the small field and Marker positioning geometry
                    dx = int(smField_num/2)*sm_spacing
                    dy = -dx
                    sm_orig_x = lg_orig_x - (1-sm_col)*dx
                    sm_orig_y = lg_orig_y - (1-sm_row)*dy
                    smMark_margin = 10
                    smMarkerPosition = smField_size/2 - smMark_margin


                    #Small Field
                    #smField_name = "SmallField_"+str(sm_row+1)+"x"+str(sm_col+1)
                    smField = Frame("SmallField", (smField_size, smField_size), [])
                    smField.make_slits(length, _width[sm_col], num_slit[sm_col], _pitch[sm_row], rot_angle, layers=l_smBeam)
                    smField.make_align_markers(2., 20., (smMarkerPosition, smMarkerPosition), l_lgBeam, joy_markers=True)

                    # Finger Contacts (for the moment the distance between central fingers is set to)
                    smField.make_finger_contacts(length, _length[sm_row]/2, num_slit[sm_col], _pitch[sm_row], rot_angle, layers= l_PHBeam) #ADJUST LAYERS
                    
                    #Material reservoir
                    smField.make_slits_reservoir(num_slit[sm_col], _pitch[sm_row], _length[sm_row], layers = l_smBeam)
                    
                    # Outer Labels
                    if sm_row == 0:
                        if num_slit[sm_col] == 1:
                            sm_label_top = "w = " + str(_width[sm_col])
                        else:
                            sm_label_top = "n = " + str(num_slit[sm_col])
                        smField_label_top = Label(sm_label_top, 10., position = (sm_orig_x-30, sm_orig_y - dy), layer=l_lgBeam)
                        topCell.add(smField_label_top)

                    if sm_col == 0:
                        if num_slit[sm_col] == 1:
                            sm_label_lat = "l = " + str(_length[sm_row])
                        else:
                            sm_label_lat = "p = " + str(_pitch[sm_row])
                        smField_label_lat = Label(sm_label_lat, 10., position = (sm_orig_x-dx, sm_orig_y - 30), angle = 90, layer=l_lgBeam)
                        topCell.add(smField_label_lat)

                    topCell.add(smField, origin=(sm_orig_x, sm_orig_y))

        topCell.add(lgField, origin=(lg_orig_x, lg_orig_y)) 
        topCell.add(lgField_label)
                    

# %%Create the layout and output GDS file
wafer = MBE100Wafer('MembranesWafer', cells=[topCell])
fileID = "GenChar_"
filestring = fileID + str(waferVer) + '_' + WAFER_ID + '_' + date.today().strftime("%d%m%Y")


# Output the whole wafer 
layout = Layout('LIBRARY')
layout.add(wafer)
# layout.show()
layout.save(filestring.replace(' ', '_') + '_wafer' + '.gds')

# Output single chip 
cell_layout = Layout('LIBRARY')
cell_layout.add(wafer.blocks[0])
# cell_layout.show()
cell_layout.save(filestring.replace(' ', '_') + '_chip' + '.gds')

# Output topCell 
layout_field = Layout('LIBRARY')
layout_field.add(topCell)
# layout_field.show()
layout_field.save(filestring.replace(' ', '_') + '_smallFields.gds')
