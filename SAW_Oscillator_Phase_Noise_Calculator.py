####################################################################################################
#
# SAW_Oscillator_Phase_Noise_Calculator.py
#
# Tool for calculating the phase/frequency noise in a closed-loop surface acoustic wave (SAW) 
# oscillator where the SAW device can be a (two-port) delay line or a resonator.
#
# Phillip Durdaut, 2021-09-23
#
####################################################################################################
import os 
import PySimpleGUI as sg
import matplotlib
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{mathrsfs}'
matplotlib.rcParams['figure.figsize'] = (7, 6)
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, FigureCanvasAgg
from matplotlib.figure import Figure
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
from PIL import ImageTk

# global constants
SYSTEM_COLOR = '#f0f0f0'
LABEL_WIDTH = 30
INPUT_WIDTH = 13
RADIO_WIDTH = INPUT_WIDTH - 5
GROUP_ID_SAW_RADIO = 1
SAW_TYPE_DELAY_LINE = 0
SAW_TYPE_RESONATOR = 1

# global variables
f = 0
start_frequency_Hz = 0
stop_frequency_Hz = 0
number_of_points = 0

generic_saw_type = 0
generic_temperature_K = 0

saw_delay_line_synchronous_frequency_MHz = 0
saw_delay_line_insertion_loss_dB = 0
saw_delay_line_group_delay_ns = 0
saw_delay_line_m3_dB_bandwidth_MHz = 0
saw_delay_line_flicker_phase_noise_coefficient_radSq = 0
saw_delay_line_white_phase_noise_coefficient_radSqPerHz = 0
saw_delay_line_noise_figure_dB = 0
saw_delay_line_noise_factor = 0

saw_resonator_center_frequency_MHz = 0
saw_resonator_insertion_loss_dB = 0
saw_resonator_quality_factor = 0
saw_resonator_m3_dB_bandwidth_MHz = 0
saw_resonator_flicker_phase_noise_coefficient_radSq = 0
saw_resonator_flicker_phase_noise_coefficient_radSqPerHz = 0
saw_resonator_noise_figure_dB = 0
saw_resonator_noise_factor = 0

amplifier_output_power_level_dBm = 0
amplifier_output_power_level_W = 0
amplifier_flicker_phase_noise_coefficient_radSq = 0
amplifier_flicker_phase_noise_coefficient_radSqPerHz = 0
amplifier_noise_figure_dB = 0
amplifier_noise_factor = 0

power_splitter_flicker_phase_noise_coefficient_radSq = 0
power_splitter_flicker_phase_noise_coefficient_radSqPerHz = 0
power_splitter_noise_figure_dB = 0
power_splitter_noise_factor = 0

oscillator_output_power_level_dBm = 0
oscillator_output_power_level_W = 0
saw_input_power_level_dBm = 0
saw_input_power_level_W = 0
saw_output_power_level_dBm = 0
saw_output_power_level_W = 0

kB = 1.38064852e-23
N_J = 0
N_dBmPerHz = 0

bm1_saw = 0
b0_saw = 0
Sphi_saw = 0

bm1_amplifier = 0
b0_amplifier = 0
Sphi_amplifier = 0

bm1_power_splitter = 0
b0_power_splitter = 0
Sphi_power_splitter = 0

Sphi_CL_saw = 0
Sphi_CL_amplifier = 0
Sphi_CL_power_splitter = 0
Sphi_CL = 0

Sf_CL_saw = 0
Sf_CL_amplifier = 0
Sf_CL_power_splitter = 0
Sf_CL = 0

# function that gets the absolute path to a resource (works in IDE and for PyInstaller)
def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

# function for analyzing and storing all given input parameters
def checkAndGetAllInputsParameters(values):

    # definitions of global variables
    global f
    global start_frequency_Hz
    global stop_frequency_Hz
    global number_of_points

    global generic_saw_type
    global generic_temperature_K

    global saw_delay_line_synchronous_frequency_MHz
    global saw_delay_line_insertion_loss_dB
    global saw_delay_line_group_delay_ns
    global saw_delay_line_m3_dB_bandwidth_MHz
    global saw_delay_line_flicker_phase_noise_coefficient_radSq

    global saw_resonator_center_frequency_MHz
    global saw_resonator_insertion_loss_dB
    global saw_resonator_quality_factor
    global saw_resonator_flicker_phase_noise_coefficient_radSq

    global amplifier_output_power_level_dBm
    global amplifier_flicker_phase_noise_coefficient_radSq
    global amplifier_noise_figure_dB

    global power_splitter_flicker_phase_noise_coefficient_radSq
    global power_splitter_noise_figure_dB
    
    # read parameters from input fields
    start_frequency_Hz = float(values['key_start_frequency_Hz'])
    stop_frequency_Hz = float(values['key_stop_frequency_Hz'])
    number_of_points = int(values['key_number_of_points'])
    f = np.logspace(np.log10(start_frequency_Hz), np.log10(stop_frequency_Hz), number_of_points)

    if (bool(values['key_generic_saw_type_delay_line']) == True):
        generic_saw_type = SAW_TYPE_DELAY_LINE
    elif (bool(values['key_generic_saw_type_resonator']) == True):
        generic_saw_type = SAW_TYPE_RESONATOR
    
    generic_temperature_K = float(values['key_generic_temperature_K'])

    saw_delay_line_synchronous_frequency_MHz = float(values['key_saw_delay_line_synchronous_frequency_MHz'])
    saw_delay_line_insertion_loss_dB = float(values['key_saw_delay_line_insertion_loss_dB'])
    saw_delay_line_group_delay_ns = float(values['key_saw_delay_line_group_delay_ns'])
    saw_delay_line_m3_dB_bandwidth_MHz = float(values['key_saw_delay_line_m3_dB_bandwidth_MHz'])
    saw_delay_line_flicker_phase_noise_coefficient_radSq = float(values['key_saw_delay_line_flicker_phase_noise_coefficient_radSq'])
    
    saw_resonator_center_frequency_MHz = float(values['key_saw_resonator_center_frequency_MHz'])
    saw_resonator_insertion_loss_dB = float(values['key_saw_resonator_insertion_loss_dB'])
    saw_resonator_quality_factor = float(values['key_saw_resonator_quality_factor'])
    saw_resonator_flicker_phase_noise_coefficient_radSq = float(values['key_saw_resonator_flicker_phase_noise_coefficient_radSq'])

    amplifier_output_power_level_dBm = float(values['key_amplifier_output_power_level_dBm'])
    amplifier_flicker_phase_noise_coefficient_radSq = float(values['key_amplifier_flicker_phase_noise_coefficient_radSq'])
    amplifier_noise_figure_dB = float(values['key_amplifier_noise_figure_dB'])
    
    power_splitter_flicker_phase_noise_coefficient_radSq = float(values['key_power_splitter_flicker_phase_noise_coefficient_radSq'])
    power_splitter_noise_figure_dB = float(values['key_power_splitter_noise_figure_dB'])

# function for performing all the calculations
def performAllCalculations():

    # definitions of global variables
    global generic_saw_type
    global generic_temperature_K

    global saw_delay_line_synchronous_frequency_MHz
    global saw_delay_line_insertion_loss_dB
    global saw_delay_line_group_delay_ns
    global saw_delay_line_m3_dB_bandwidth_MHz
    global saw_delay_line_flicker_phase_noise_coefficient_radSq
    global saw_delay_line_white_phase_noise_coefficient_radSqPerHz
    global saw_delay_line_noise_figure_dB
    global saw_delay_line_noise_factor

    global saw_resonator_center_frequency_MHz
    global saw_resonator_insertion_loss_dB
    global saw_resonator_quality_factor
    global saw_resonator_m3_dB_bandwidth_MHz
    global saw_resonator_flicker_phase_noise_coefficient_radSq
    global saw_resonator_flicker_phase_noise_coefficient_radSqPerHz
    global saw_resonator_noise_figure_dB
    global saw_resonator_noise_factor

    global amplifier_output_power_level_dBm
    global amplifier_flicker_phase_noise_coefficient_radSq
    global amplifier_flicker_phase_noise_coefficient_radSqPerHz
    global amplifier_noise_figure_dB
    global amplifier_noise_factor

    global power_splitter_flicker_phase_noise_coefficient_radSq
    global power_splitter_flicker_phase_noise_coefficient_radSqPerHz
    global power_splitter_noise_figure_dB
    global power_splitter_noise_factor

    global oscillator_output_power_level_dBm
    global oscillator_output_power_level_W
    global saw_input_power_level_dBm
    global saw_input_power_level_W
    global saw_output_power_level_dBm
    global saw_output_power_level_W

    global N_J
    global N_dBmPerHz

    global bm1_saw
    global b0_saw
    global Sphi_saw

    global bm1_amplifier
    global b0_amplifier
    global Sphi_amplifier

    global bm1_power_splitter
    global b0_power_splitter
    global Sphi_power_splitter 

    global Sphi_CL_saw
    global Sphi_CL_amplifier
    global Sphi_CL_power_splitter
    global Sphi_CL

    global Sf_CL_saw
    global Sf_CL_amplifier
    global Sf_CL_power_splitter
    global Sf_CL

    # calculate power levels
    oscillator_output_power_level_dBm = amplifier_output_power_level_dBm - 3
    saw_input_power_level_dBm = amplifier_output_power_level_dBm - 3

    if (generic_saw_type == SAW_TYPE_DELAY_LINE):
        saw_output_power_level_dBm = saw_input_power_level_dBm - saw_delay_line_insertion_loss_dB
    elif (generic_saw_type == SAW_TYPE_RESONATOR):
        saw_output_power_level_dBm = saw_input_power_level_dBm - saw_resonator_insertion_loss_dB

    amplifier_output_power_level_W = 1e-3 * np.power(10, amplifier_output_power_level_dBm/10)
    oscillator_output_power_level_W = 1e-3 * np.power(10, oscillator_output_power_level_dBm/10)
    saw_input_power_level_W = 1e-3 * np.power(10, saw_input_power_level_dBm/10)
    saw_output_power_level_W = 1e-3 * np.power(10, saw_output_power_level_dBm/10)

    # calculate overall noise energy
    N_J = kB * generic_temperature_K
    N_dBmPerHz = 10 * np.log10(N_J/1e-3)

    # calculate further SAW delay line parameters and update GUI
    saw_delay_line_noise_figure_dB = saw_delay_line_insertion_loss_dB
    saw_delay_line_noise_factor = np.power(10, saw_delay_line_noise_figure_dB/10)
    saw_delay_line_white_phase_noise_coefficient_radSqPerHz = float(format(saw_delay_line_noise_factor * N_J / saw_input_power_level_W, '.2e'))
    window['key_saw_delay_line_noise_figure_dB'].update(str(saw_delay_line_noise_figure_dB))
    window['key_saw_delay_line_white_phase_noise_coefficient_radSqPerHz'].update(str(saw_delay_line_white_phase_noise_coefficient_radSqPerHz))

    # calculate further SAW resonator parameters and update GUI
    saw_resonator_noise_figure_dB = saw_resonator_insertion_loss_dB
    saw_resonator_noise_factor = np.power(10, saw_resonator_noise_figure_dB/10)
    saw_resonator_white_phase_noise_coefficient_radSqPerHz = float(format(saw_resonator_noise_factor * N_J / saw_input_power_level_W, '.2e'))
    saw_resonator_m3_dB_bandwidth_MHz = saw_resonator_center_frequency_MHz / saw_resonator_quality_factor
    window['key_saw_resonator_noise_figure_dB'].update(str(saw_resonator_noise_figure_dB))
    window['key_saw_resonator_white_phase_noise_coefficient_radSqPerHz'].update(str(saw_resonator_white_phase_noise_coefficient_radSqPerHz))
    window['key_saw_resonator_m3_dB_bandwidth_MHz'].update(str(saw_resonator_m3_dB_bandwidth_MHz))

    # calculate further amplifier parameters and update GUI
    amplifier_noise_factor = np.power(10, amplifier_noise_figure_dB/10)
    amplifier_white_phase_noise_coefficient_radSqPerHz = float(format(amplifier_noise_factor * N_J / saw_output_power_level_W, '.2e'))
    window['key_amplifier_white_phase_noise_coefficient_radSqPerHz'].update(str(amplifier_white_phase_noise_coefficient_radSqPerHz))

    # calculate further power splitter parameters and update GUI
    power_splitter_noise_factor = np.power(10, power_splitter_noise_figure_dB/10)
    power_splitter_white_phase_noise_coefficient_radSqPerHz = float(format(power_splitter_noise_factor * N_J / amplifier_output_power_level_W, '.2e'))
    window['key_power_splitter_white_phase_noise_coefficient_radSqPerHz'].update(str(power_splitter_white_phase_noise_coefficient_radSqPerHz))

    # calculate SAW phase noise
    if (generic_saw_type == SAW_TYPE_DELAY_LINE):
        bm1_saw = saw_delay_line_flicker_phase_noise_coefficient_radSq
        b0_saw = saw_delay_line_white_phase_noise_coefficient_radSqPerHz
    elif (generic_saw_type == SAW_TYPE_RESONATOR == True):
        bm1_saw = saw_resonator_flicker_phase_noise_coefficient_radSq
        b0_saw = saw_resonator_white_phase_noise_coefficient_radSqPerHz
    Sphi_saw = bm1_saw * np.power(f,-1) + b0_saw * np.power(f,0)

    # calculate amplifier phase noise
    bm1_amplifier = amplifier_flicker_phase_noise_coefficient_radSq
    b0_amplifier = amplifier_white_phase_noise_coefficient_radSqPerHz
    Sphi_amplifier = bm1_amplifier * np.power(f,-1) + b0_amplifier * np.power(f,0)

    # calculate power splitter phase noise
    bm1_power_splitter = power_splitter_flicker_phase_noise_coefficient_radSq
    b0_power_splitter = power_splitter_white_phase_noise_coefficient_radSqPerHz
    Sphi_power_splitter = bm1_power_splitter * np.power(f,-1) + b0_power_splitter * np.power(f,0)

    # calculate oscillator phase noise
    if (generic_saw_type == SAW_TYPE_DELAY_LINE):        
        HCLMagSq_N = 1 + np.power(2*f/(saw_delay_line_m3_dB_bandwidth_MHz*1e6), 2)
        HCLMagSq_D = np.power(2*f/(saw_delay_line_m3_dB_bandwidth_MHz*1e6), 2) + 2 * (1 - np.sqrt(1 + np.power(2*f/(saw_delay_line_m3_dB_bandwidth_MHz*1e6), 2)) * np.cos(2*np.pi*f*saw_delay_line_group_delay_ns/1e9))
        HCLMagSq = HCLMagSq_N / HCLMagSq_D
    elif (generic_saw_type == SAW_TYPE_RESONATOR == True):
        HCLMagSq = 1 + np.power((saw_resonator_m3_dB_bandwidth_MHz*1e6)/(2*saw_resonator_quality_factor*f), 2)

    Sphi_CL_saw = HCLMagSq * Sphi_saw
    Sphi_CL_amplifier = HCLMagSq * Sphi_amplifier
    Sphi_CL_power_splitter = HCLMagSq * Sphi_power_splitter
    Sphi_CL = HCLMagSq * (Sphi_saw + Sphi_amplifier + Sphi_power_splitter)

    # calculate oscillator frequency noise
    Sf_CL_saw = np.power(f, 2) * Sphi_CL_saw
    Sf_CL_amplifier = np.power(f, 2) * Sphi_CL_amplifier
    Sf_CL_power_splitter = np.power(f, 2) * Sphi_CL_power_splitter
    Sf_CL = np.power(f, 2) * Sphi_CL

# function plotting the power levels
def plotAllPowerLevels():

    image = Image.open(resource_path('block_diagram.png'))
    ImageDraw.Draw(image).text((660, 230), text=str(oscillator_output_power_level_dBm) + ' dBm', font=ImageFont.truetype("arial.ttf", 16), fill='green')
    ImageDraw.Draw(image).text((400, 100), text=str(saw_input_power_level_dBm) + ' dBm', font=ImageFont.truetype("arial.ttf", 16), fill='green')
    ImageDraw.Draw(image).text((400, 315), text=str(amplifier_output_power_level_dBm) + ' dBm', font=ImageFont.truetype("arial.ttf", 16), fill='green')
    ImageDraw.Draw(image).text((175, 315), text=str(saw_output_power_level_dBm) + ' dBm', font=ImageFont.truetype("arial.ttf", 16), fill='green')
    window['key_image_block_diagram'].update(data=ImageTk.PhotoImage(image))

# function for plotting the components phase noise
def plotComponentsPhaseNoise(firstCall):

    if firstCall == True:
        
        # initialize empty figure on program start
        figure = plt.figure()        
        figure.patch.set_facecolor(SYSTEM_COLOR)
        canvas = window.FindElement('key_canvas_components_phase_noise')
        figureCanvas = FigureCanvasTkAgg(figure, master=canvas.TKCanvas)
        figureCanvas.draw()
        figureCanvas.get_tk_widget().pack()
        return figureCanvas

    else:

        figure = plt.figure()
        ax1 = figure.add_subplot(111)
        figure.patch.set_facecolor(SYSTEM_COLOR)
        canvas = window.FindElement('key_canvas_components_phase_noise')
        figureCanvas = FigureCanvasTkAgg(figure, master=canvas.TKCanvas)

        # SAW phase noise
        if (generic_saw_type == SAW_TYPE_DELAY_LINE):
            ax1.loglog(f, Sphi_saw, label='SAW Delay Line')
        elif (generic_saw_type == SAW_TYPE_RESONATOR == True):
            ax1.loglog(f, Sphi_saw, label='SAW Delay Resonator')

        # amplifier phase noise
        ax1.loglog(f, Sphi_amplifier, label='Amplifier')

        # power splitter phase noise
        ax1.loglog(f, Sphi_power_splitter, label='Power Splitter')

        # further plot settings
        plt.xlabel('Frequency [Hz]')
        plt.xlim(f[0], f[-1])

        # first y axis (rad^2/Hz)
        exp_bottom = np.floor(np.log10(ax1.get_ylim())[0])
        exp_top = np.ceil(np.log10(ax1.get_ylim())[1])
        ax1.set_ylim(np.power(10, exp_bottom), np.power(10, exp_top))
        ax1.set_yticks(np.logspace(exp_bottom, exp_top, int((exp_top-exp_bottom)+1)))

        ax1.grid()
        ax1.set_ylabel('One-Sided Phase Noise Power Spectral Density $S_{\\varphi}(f)$ [rad$^2$/Hz]')
        ax1.legend()

        # second y axis (dBc/Hz)
        ax2 = ax1.twinx()
        ax2.set_ylabel('Two-Sided Phase Noise Power Spectral Density $\\mathscr{L}(f)$ [dBc/Hz]')
        ax2.set_ylim(10 * np.log10(ax1.get_ylim()) - 3)        
        
        # create figure
        figureCanvas.draw()
        figureCanvas.get_tk_widget().pack()

        return figureCanvas

# function for plotting the oscillator phase noise
def plotOscillatorPhaseNoise(firstCall):

    if firstCall == True:
        
        # initialize empty figure on program start
        figure = plt.figure()        
        figure.patch.set_facecolor(SYSTEM_COLOR)
        canvas = window.FindElement('key_canvas_oscillator_phase_noise')
        figureCanvas = FigureCanvasTkAgg(figure, master=canvas.TKCanvas)
        figureCanvas.draw()
        figureCanvas.get_tk_widget().pack()
        return figureCanvas

    else:

        figure = plt.figure()
        ax1 = figure.add_subplot(111)
        figure.patch.set_facecolor(SYSTEM_COLOR)
        canvas = window.FindElement('key_canvas_oscillator_phase_noise')
        figureCanvas = FigureCanvasTkAgg(figure, master=canvas.TKCanvas)

        # SAW phase noise
        if (generic_saw_type == SAW_TYPE_DELAY_LINE):
            ax1.loglog(f, Sphi_CL_saw, label='SAW Delay Line')
        elif (generic_saw_type == SAW_TYPE_RESONATOR == True):
            ax1.loglog(f, Sphi_CL_saw, label='SAW Delay Resonator')

        # amplifier phase noise
        ax1.loglog(f, Sphi_CL_amplifier, label='Amplifier')

        # power splitter phase noise
        ax1.loglog(f, Sphi_CL_power_splitter, label='Power Splitter')

        # total phase noise
        ax1.loglog(f, Sphi_CL, 'k--', label='Total Phase Noise')

        # further plot settings
        plt.xlabel('Frequency [Hz]')
        plt.xlim(f[0], f[-1])

        # first y axis (rad^2/Hz)
        exp_bottom = np.floor(np.log10(ax1.get_ylim())[0])
        exp_top = np.ceil(np.log10(ax1.get_ylim())[1])
        ax1.set_ylim(np.power(10, exp_bottom), np.power(10, exp_top))
        ax1.set_yticks(np.logspace(exp_bottom, exp_top, int((exp_top-exp_bottom)+1)))

        ax1.grid()
        ax1.set_ylabel('One-Sided Phase Noise Power Spectral Density $S_{\\varphi}(f)$ [rad$^2$/Hz]')
        ax1.legend()

        # second y axis (dBc/Hz)
        ax2 = ax1.twinx()
        ax2.set_ylabel('Two-Sided Phase Noise Power Spectral Density $\\mathscr{L}(f)$ [dBc/Hz]')
        ax2.set_ylim(10 * np.log10(ax1.get_ylim()) - 3)        
        
        # create figure
        figureCanvas.draw()
        figureCanvas.get_tk_widget().pack()

        return figureCanvas

# function for plotting the oscillator frequency noise
def plotOscillatorFrequencyNoise(firstCall):

    if firstCall == True:
        
        # initialize empty figure on program start
        figure = plt.figure()        
        figure.patch.set_facecolor(SYSTEM_COLOR)
        canvas = window.FindElement('key_canvas_oscillator_frequency_noise')
        figureCanvas = FigureCanvasTkAgg(figure, master=canvas.TKCanvas)
        figureCanvas.draw()
        figureCanvas.get_tk_widget().pack()
        return figureCanvas

    else:

        figure = plt.figure()
        ax1 = figure.add_subplot(111)
        figure.patch.set_facecolor(SYSTEM_COLOR)
        canvas = window.FindElement('key_canvas_oscillator_frequency_noise')
        figureCanvas = FigureCanvasTkAgg(figure, master=canvas.TKCanvas)

        # SAW frequency noise
        if (generic_saw_type == SAW_TYPE_DELAY_LINE):
            ax1.loglog(f, Sf_CL_saw, label='SAW Delay Line')
        elif (generic_saw_type == SAW_TYPE_RESONATOR == True):
            ax1.loglog(f, Sf_CL_saw, label='SAW Delay Resonator')

        # amplifier frequency noise
        ax1.loglog(f, Sf_CL_amplifier, label='Amplifier')

        # power splitter frequency noise
        ax1.loglog(f, Sf_CL_power_splitter, label='Power Splitter')

        # total frequency noise
        ax1.loglog(f, Sf_CL, 'k--', label='Total Frequency Noise')

        # further plot settings
        plt.xlabel('Frequency [Hz]')
        plt.xlim(f[0], f[-1])

        # first y axis (rad^2/Hz)
        exp_bottom = np.floor(np.log10(ax1.get_ylim())[0])
        exp_top = np.ceil(np.log10(ax1.get_ylim())[1])
        ax1.set_ylim(np.power(10, exp_bottom), np.power(10, exp_top))
        ax1.set_yticks(np.logspace(exp_bottom, exp_top, int((exp_top-exp_bottom)+1)))

        ax1.grid()
        ax1.set_ylabel('One-Sided Frequency Noise Power Spectral Density $S_{\\mathrm{f}}(f)$ [Hz$^2$/Hz]')
        ax1.legend()      
        
        # create figure
        figureCanvas.draw()
        figureCanvas.get_tk_widget().pack()

        return figureCanvas
        
# set standard color theme
sg.theme('SystemDefault')

# load block diagram
block_diagram = sg.Image(filename=resource_path('block_diagram.png'), key='key_image_block_diagram')

# layout of first column
column_first = [  
    [ sg.Text('Generic', font=('Arial', 10, 'bold')) ],
    [ sg.Text('SAW Type:', size=(LABEL_WIDTH,1)), sg.Radio('Delay Line', size=(RADIO_WIDTH,1), group_id=GROUP_ID_SAW_RADIO, default=True, key='key_generic_saw_type_delay_line') ],
    [ sg.Text('', size=(LABEL_WIDTH,1)), sg.Radio('Resonator', size=(RADIO_WIDTH,1), group_id=GROUP_ID_SAW_RADIO, default=False, key='key_generic_saw_type_resonator') ],
    [ sg.Text('Temperature [K]:', size=(LABEL_WIDTH,1)), sg.Input('290', size=(INPUT_WIDTH,1), key='key_generic_temperature_K') ],
    [ sg.Text('') ],
    [ sg.Text('SAW Delay Line', font=('Arial', 10, 'bold')) ],
    [ sg.Text('Synchronous Frequency [MHz]:', size=(LABEL_WIDTH,1)), sg.Input('100', size=(INPUT_WIDTH,1), key='key_saw_delay_line_synchronous_frequency_MHz') ],
    [ sg.Text('Insertion Loss [dB]:', size=(LABEL_WIDTH,1)), sg.Input('20', size=(INPUT_WIDTH,1), key='key_saw_delay_line_insertion_loss_dB') ],
    [ sg.Text('Group Delay [ns]:', size=(LABEL_WIDTH,1)), sg.Input('1500', size=(INPUT_WIDTH,1), key='key_saw_delay_line_group_delay_ns') ],
    [ sg.Text('-3 dB Bandwidth [MHz]:', size=(LABEL_WIDTH,1)), sg.Input('5', size=(INPUT_WIDTH,1), key='key_saw_delay_line_m3_dB_bandwidth_MHz') ],
    [ sg.Text('Flicker Phase Noise Coefficient [rad^2]:', size=(LABEL_WIDTH,1)), sg.Input('1e-13', size=(INPUT_WIDTH,1), key='key_saw_delay_line_flicker_phase_noise_coefficient_radSq') ],
    [ sg.Text('White Phase Noise Coefficient [rad^2/Hz]:', size=(LABEL_WIDTH,1)), sg.Input('', size=(INPUT_WIDTH,1), key='key_saw_delay_line_white_phase_noise_coefficient_radSqPerHz', readonly=True) ],
    [ sg.Text('Noise Figure [dB]:', size=(LABEL_WIDTH,1)), sg.Input('', size=(INPUT_WIDTH,1), key='key_saw_delay_line_noise_figure_dB', readonly=True) ],
    [ sg.Text('') ],
    [ sg.Text('SAW Resonator', font=('Arial', 10, 'bold')) ],
    [ sg.Text('Center Frequency [MHz]:', size=(LABEL_WIDTH,1)), sg.Input('100', size=(INPUT_WIDTH,1), key='key_saw_resonator_center_frequency_MHz') ],
    [ sg.Text('Insertion Loss [dB]:', size=(LABEL_WIDTH,1)), sg.Input('2', size=(INPUT_WIDTH,1), key='key_saw_resonator_insertion_loss_dB') ],
    [ sg.Text('Quality Factor [1]:', size=(LABEL_WIDTH,1)), sg.Input('500', size=(INPUT_WIDTH,1), key='key_saw_resonator_quality_factor') ],
    [ sg.Text('-3 dB Bandwidth [MHz]:', size=(LABEL_WIDTH,1)), sg.Input('', size=(INPUT_WIDTH,1), key='key_saw_resonator_m3_dB_bandwidth_MHz', readonly=True) ],
    [ sg.Text('Flicker Phase Noise Coefficient [rad^2]:', size=(LABEL_WIDTH,1)), sg.Input('1e-14', size=(INPUT_WIDTH,1), key='key_saw_resonator_flicker_phase_noise_coefficient_radSq') ],
    [ sg.Text('White Phase Noise Coefficient [rad^2/Hz]:', size=(LABEL_WIDTH,1)), sg.Input('', size=(INPUT_WIDTH,1), key='key_saw_resonator_white_phase_noise_coefficient_radSqPerHz', readonly=True) ],
    [ sg.Text('Noise Figure [dB]:', size=(LABEL_WIDTH,1)), sg.Input('', size=(INPUT_WIDTH,1), key='key_saw_resonator_noise_figure_dB', readonly=True) ]
]

# layout of third column
column_third = [    
    [ sg.Text('Amplifier', font=('Arial', 10, 'bold')) ],
    [ sg.Text('Output Power Level [dBm]:', size=(LABEL_WIDTH,1)), sg.Input('10', size=(INPUT_WIDTH,1), key='key_amplifier_output_power_level_dBm') ],
    [ sg.Text('Flicker Phase Noise Coefficient [rad^2]:', size=(LABEL_WIDTH,1)), sg.Input('1e-12', size=(INPUT_WIDTH,1), key='key_amplifier_flicker_phase_noise_coefficient_radSq') ],
    [ sg.Text('White Phase Noise Coefficient [rad^2/Hz]:', size=(LABEL_WIDTH,1)), sg.Input('', size=(INPUT_WIDTH,1), key='key_amplifier_white_phase_noise_coefficient_radSqPerHz', readonly=True) ],
    [ sg.Text('Noise Figure [dB]:', size=(LABEL_WIDTH,1)), sg.Input('2.5', size=(INPUT_WIDTH,1), key='key_amplifier_noise_figure_dB') ],
    [ sg.Text('') ],
    [ sg.Text('Power Splitter', font=('Arial', 10, 'bold')) ],
    [ sg.Text('Flicker Phase Noise Coefficient [rad^2]:', size=(LABEL_WIDTH,1)), sg.Input('1e-17', size=(INPUT_WIDTH,1), key='key_power_splitter_flicker_phase_noise_coefficient_radSq') ],
    [ sg.Text('White Phase Noise Coefficient [rad^2/Hz]:', size=(LABEL_WIDTH,1)), sg.Input('', size=(INPUT_WIDTH,1), key='key_power_splitter_white_phase_noise_coefficient_radSqPerHz', readonly=True) ],
    [ sg.Text('Noise Figure [dB]:', size=(LABEL_WIDTH,1)), sg.Input('3', size=(INPUT_WIDTH,1), key='key_power_splitter_noise_figure_dB', readonly=True) ],
    [ sg.Text('') ],
    [ sg.Text('Plot Settings', font=('Arial', 10, 'bold')) ],
    [ sg.Text('Start Frequency [Hz]:', size=(LABEL_WIDTH,1)), sg.Input('1', size=(INPUT_WIDTH,1), key='key_start_frequency_Hz') ],
    [ sg.Text('Stop Frequency [Hz]:', size=(LABEL_WIDTH,1)), sg.Input('1e6', size=(INPUT_WIDTH,1), key='key_stop_frequency_Hz') ],
    [ sg.Text('Number of Points [1]:', size=(LABEL_WIDTH,1)), sg.Input('6001', size=(INPUT_WIDTH,1), key='key_number_of_points') ],
    [ sg.Text('') ],
    [ sg.Button('Calculate and Plot'), sg.Button('Exit') ],
    [ sg.Text('') ],
    [ sg.Text('') ],
    [ sg.Text('Contact:') ],
    [ sg.Text('Phillip Durdaut') ],
    [ sg.Text('pdurdaut@googlemail.com') ],
    [ sg.Text('https://github.com/yildi1337') ]
]

# tab: block diagram
tab_block_diagram = [ 
    [ block_diagram ]
]

# tab: components phase noise
tab_components_phase_noise = [ 
    [ sg.Canvas(key='key_canvas_components_phase_noise') ]
]

# tab: oscillator phase noise
tab_oscillator_phase_noise = [ 
    [ sg.Canvas(key='key_canvas_oscillator_phase_noise') ]
]

# tab: oscillator frequency noise
tab_oscillator_frequency_noise = [ 
    [ sg.Canvas(key='key_canvas_oscillator_frequency_noise') ]
]

tabgroup = [
    [ 
        sg.TabGroup(
            [
                [
                    sg.Tab('Block Diagram', tab_block_diagram, element_justification='center'),
                    sg.Tab('Components Phase Noise', tab_components_phase_noise, element_justification='center'),
                    sg.Tab('Oscillator Phase Noise', tab_oscillator_phase_noise, element_justification='center'),
                    sg.Tab('Oscillator Frequency Noise', tab_oscillator_frequency_noise, element_justification='center')
                ]
            ], border_width=0, tab_background_color=SYSTEM_COLOR, selected_background_color='White', background_color=SYSTEM_COLOR
        )
    ]
] 

# overall layout
layout = [
    [
        sg.Column(column_first, element_justification='center'),
        sg.VSeperator(),
        sg.Column(tabgroup),
        sg.VSeperator(),
        sg.Column(column_third, element_justification='center')
    ]
]

# create new window
window = sg.Window("SAW Oscillator Phase Noise Calculator (2021)", layout, finalize=True, margins=(25, 25), location=(100,150))

# create initial plots
figureCanvasComponentsPhaseNoise = plotComponentsPhaseNoise(True)
figureCanvasOscillatorPhaseNoise = plotOscillatorPhaseNoise(True)
figureCanvasOscillatorFrequencyNoise = plotOscillatorFrequencyNoise(True)

# run
while True:

    event, values = window.read()

    if event == "Calculate and Plot":

        checkAndGetAllInputsParameters(values)
        performAllCalculations()
        plotAllPowerLevels()
        
        figureCanvasComponentsPhaseNoise.get_tk_widget().pack_forget()
        figureCanvasComponentsPhaseNoise = plotComponentsPhaseNoise(False)

        figureCanvasOscillatorPhaseNoise.get_tk_widget().pack_forget()
        figureCanvasOscillatorPhaseNoise = plotOscillatorPhaseNoise(False)

        figureCanvasOscillatorFrequencyNoise.get_tk_widget().pack_forget()
        figureCanvasOscillatorFrequencyNoise = plotOscillatorFrequencyNoise(False)

    if event == "Exit" or event == sg.WIN_CLOSED:
        break

window.close()
