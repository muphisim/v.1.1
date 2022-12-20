import tkinter as tk
from tkinter import *
from tkinter import ttk
from tkinter import messagebox as MessageBox
import Pmw, sys
#import webbrowser

# ===================================================================
# STYLES FOR MuPhiSim USER INTERFACE
# ===================================================================
#002147
def MuPhiSimHeader(place, string, rownum, colnum):
	tk.Label(place, text=string, bg="#E8E8E8", fg="black", height=2, width=10, font = ("Open Sans",14, 'bold'), anchor=W, justify=LEFT).grid(row=rownum, column=colnum, sticky="nsew", columnspan=2)

def MuPhiSimSubHeader(place, string, rownum, colnum):
	tk.Label(place, text=string, fg="black", bg = "#E8E8E8", height=1, font = ("Open Sans",12, 'bold'), anchor=W, justify=LEFT).grid(row=rownum, column=colnum, sticky="nsew", pady = 5)

def MuPhiSimBodyLabel(place, string, rownum, colnum):
	tk.Label(place, text=string, fg="black", bg = "#E8E8E8", font = ("Open Sans",10, 'bold'), justify= CENTER, wraplength=450).grid(row=rownum, column=colnum, sticky="nsew", pady = 5)

def MuPhiSimButton1(place, text, command):
	Button(place, text=text,command=command, bg="#002147", font = ("Open Sans",10))
