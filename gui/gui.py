#!/usr/bin/env python
#

import os
import subprocess
import sys
from pathlib import Path
import pandas as pd
import PySimpleGUI as sg

#
#
#
# SETTINGS WINDOWS
def settings_window(settings):

    parameters_column = [
        [sg.Text("Run identifier:")],
        [sg.Input(settings["Docker"]["run_ID"], s=25, key="RUNID")],

        [sg.Text("Container:")],
        [sg.Input(settings["Docker"]["image"], s=25, key="PULL")],

        [sg.Text("Assembly options:"),
        sg.Combo(settings["Options"]["assemble"].split("|"), 
        default_value=settings["Parameters"]["assemble"], 
        s=10, key="ASSEM")],

        [sg.Text("Reordering options:"),
        sg.Combo(settings["Options"]["reorder"].split("|"), 
        default_value=settings["Parameters"]["reorder"], 
        s=10, key="REORD")],

        [sg.Text("Network options:"),
        sg.Combo(settings["Options"]["network"].split("|"), 
        default_value=settings["Parameters"]["network"], 
        s=10, key="NETW")],

    ]

    exit_column = [
        [sg.Button("Save settings", s=20)],
        [sg.Button("Exit", s=20)]
    ]

    layout = [
        [
        sg.Column(parameters_column),
        sg.VSeperator(),
        sg.Column(exit_column)
        ]
    ]

    window = sg.Window(title, layout, modal=True)

    while True:
        event, values = window.read()
        if event == sg.WINDOW_CLOSED:
            break
        if event == "Exit":
            break
        if event == "Save settings":
            # Writing changes to config file
            settings["Docker"]["run_ID"] = values["RUNID"]
            settings["Docker"]["image"] = values["PULL"]
            settings["Parameters"]["assemble"] = values["ASSEM"]
            settings["Parameters"]["reorder"] = values["REORD"]
            settings["Parameters"]["network"] = values["NETW"]

            # Displaying success message
            sg.popup_no_titlebar("Settings saved")

    window.close()

def check_window(settings):
    sg.popup_scrolled(settings, title="Current settings")

#
#
#
#
# MAIN WINDOW
def main_window():
    # Defining functions
    def validate(filepath):
            if filepath and Path(filepath).exists():
                return True
            sg.popup_error("Filepaths incorrect you phantastic nonce")
            return False

    # Producing layout columns
    files_column = [
        [sg.Text('Input folder'), sg.Input(key="IN", s=(20,2)), sg.FolderBrowse()],
        [sg.Text('Output folder'), sg.Input(key="OUT", s=(19,2)), sg.FolderBrowse()],
        [sg.Button('Change settings', s=(8,2)), sg.Button('Check settings', s=(8,2))]
    ]

    tools_column = [
        [sg.Button('Assemble reads', s=20)],
        [sg.Button('Reorder genomes', s=20)],
        [sg.Button('Produce network', s=20)],
        [sg.Button('CLOSE', s=20)]
    ]

    layout = [
        [
        sg.Column(tools_column),
        sg.VSeperator(),
        sg.Column(files_column)
        ]
    ]

    # Create window
    window = sg.Window(title, layout)
    
    # EVENTS LOOP
    while True:
        event, values = window.read()
        if event == 'CLOSE':
            window.close()
            break
        if event == sg.WIN_CLOSED:
            os.system("notify-send 'Operation cancelled'")
            window.close()
            break
        if event == 'Check settings':
            check_window(settings)
        if event == 'Change settings':
            settings_window(settings)

        # BUTTONS
        if event == 'Assemble reads':

            os.system("notify-send 'RUNNING ASSEMBLY'")
            os.system("echo 'RUNNING ASSEMBLY'")

            # STANDARD BLOCK
            Input = str(values["IN"])
            Output = str(values["OUT"])
            Image = str(settings["Docker"]["image"])
            Command = str(settings["Commands"]["assemble"])
            Parameter = str(settings["Parameters"]["assemble"])

            if validate(values["IN"]) and validate(values["OUT"]):
                subprocess.Popen(["docker run -v %s:/phanatic/IN \
                    -v %s:/phanatic/OUT \
                    %s %s %s" %(Input,Output,Image,Command,Parameter)], shell=True)

        if event == 'Reorder genomes':

            os.system("notify-send 'RUNNING ANALYSIS'")
            os.system("echo 'RUNNING ANALYSIS'")

            # STANDARD BLOCK
            Input = str(values["IN"])
            Output = str(values["OUT"])
            Database = str(Path.cwd())
            Image = str(settings["Docker"]["image"])
            Command = str(settings["Commands"]["reorder"])
            Parameter = str(settings["Parameters"]["reorder"])

            if validate(values["IN"]) and validate(values["OUT"]):
                subprocess.Popen(["docker run -v %s:/phanatic/IN \
                    -v %s:/phanatic/OUT \
                    -v %s:/phanatic/DATA \
                    %s %s %s" %(Input,Output,Database,Image,Command,Parameter)], shell=True)

        if event == 'Produce network':

            os.system("notify-send 'PRODUCING NETWORK'")
            os.system("echo 'PRODUCING NETWORK'")

            # STANDARD BLOCK
            Input = str(values["IN"])
            Output = str(values["OUT"])
            Image = str(settings["Docker"]["image"])
            Command = str(settings["Commands"]["network"])
            Parameter = str(settings["Parameters"]["network"])

            if validate(values["IN"]) and validate(values["OUT"]):
                subprocess.Popen(["docker run -v %s:/phanatic/IN \
                    -v %s:/phanatic/OUT \
                    %s %s %s" %(Input,Output,Image,Command,Parameter)], shell=True)

    window.close()

#
#
#
# Boilerplate code (?)
if __name__ == "__main__":
    # Setting config file
    config_path = Path.cwd()
    settings = sg.UserSettings(path=config_path, filename="config.ini", \
        use_config_file=True, convert_bools_and_none=True)
    # Linking config file
    title = settings["GUI"]["title"]
    theme = settings["GUI"]["theme"]
    font_family = settings["GUI"]["font_family"]
    font_size = int(settings["GUI"]["font_size"])
    # Using config file
    sg.theme(theme)
    sg.set_options(font=(font_family, font_size))
    main_window()
