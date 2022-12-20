#Parser to convert Abaqus files into MuPhiSim input files
#GUI file
#Authors: Khariton Gorbunov, Thomas Kavanagh

# ===================================================================
# 1.0 IMPORTS
# ===================================================================
from tkinter import ttk
from tkinter import *
from MuPhiSimstyles import *
from tkinter.constants import DISABLED
from tkinter.filedialog import askopenfilename, asksaveasfilename
from parser import *
# ===================================================================
# 2.0 USER INTERFACE
# ===================================================================
# Create Instance
root = tk.Tk()
root.resizable(1, 1)

# Add a title
root.title('MuPhiSim© Input')

# Set GUI Icon
imgicon = tk.PhotoImage(file='content/MuPhiSimicon.gif')
root.tk.call('wm', 'iconphoto', root._w, imgicon)
# TODO Add Mac compatability to icon - still displays py icon when ran on mac.

# Create Header
Header = MuPhiSimHeader(root, " MuPhiSim© Input", 0, 0)
#logo = tk.PhotoImage(file='content/MuPhiSim-1.gif')
#w = tk.Label(Header, image=logo)
#w.photo = logo
w = tk.Label(root, text="Hello Tkinter!")

# Set Background colour
root.configure(background="#E8E8E8")

# ===================================================================
# 2.1 Add Menu Bar
# ===================================================================

# Menu Bar Functions
def ExitApplication():
	ExitMsgBox = MessageBox.askquestion('Exit Application',
									 'Are you sure you want to exit the application?')
	if ExitMsgBox == 'yes':
		root.destroy()
	else:
		return
# TODO Add message box upon closing window via the 'X'.
#  Add flashing box to prevent 'X' being clicked straight away.

def OpenUserMan():
	webbrowser.open_new(r'content/MuPhiSimManualV2.pdf')

menubar = Menu(root)

def OpenEULA():
	EULAWindow = tk.Toplevel(root)
	EULAWindow.resizable(0, 0)
	EULAWindow.configure(background='#E8E8E8')
	EULAWindow.tk.call('wm', 'iconphoto', EULAWindow._w, imgicon)
	# TODO Add Macintosh compatability to icon - still displays py icon when ran on mac.

	MuPhiSimHeader(EULAWindow, " MuPhiSim© Input", 0, 0)
	MuPhiSimSubHeader(EULAWindow, "    End User License Agreement", 1, 0)
	MuPhiSimBodyLabel(EULAWindow,
					"To use this software you must accept the terms of the End User License Agreement (EULA):",
					3, 0)
	EULA = Pmw.ScrolledText(EULAWindow,
		borderframe=5,
		vscrollmode='dynamic',
		hscrollmode='dynamic',
		text_width=70,
		text_height=20,
		text_wrap=WORD)
	EULA.grid(row=5, column=0, padx=30, pady=10)
	EULA.insert('end', open('/home/thomas/Thomas/MuPhiSiminput/LICENSE.md', 'r').read())
	EULA.configure(text_state='disabled')
	EULAQuit = tk.Button(EULAWindow, text='Cancel', command=EULAWindow.destroy, bg="#002147",
						 fg="white",
						font=("Open Sans", 10,))
	EULAQuit.grid(pady=(0, 10))

def OpenReadME():
	ReadMeWindow = tk.Toplevel(root)
	ReadMeWindow.configure(background='#E8E8E8')
	ReadMeWindow.tk.call('wm', 'iconphoto', ReadMeWindow._w, imgicon)

	MuPhiSimHeader(ReadMeWindow, " MuPhiSim© Input", 0, 0)
	MuPhiSimSubHeader(ReadMeWindow, "    User QuickHelp (README)", 1, 0)
	MuPhiSimBodyLabel(ReadMeWindow,
					"To use this software, please read the following:",
					3, 0)
	ReadMe = Pmw.ScrolledText(ReadMeWindow,
		borderframe=5,
		vscrollmode='dynamic',
		hscrollmode='dynamic',
		text_width=70,
		text_height=20,
		text_wrap=WORD)
	ReadMe.grid(row=5, column=0, padx=30, pady=10)
	ReadMe.insert('end', open('content/README.md', 'r').read())
	ReadMe.configure(text_state='disabled')
	ReadMeQuit = tk.Button(ReadMeWindow, text='Cancel', command=ReadMeWindow.destroy, bg="#002147",
						 fg="white",
						font=("Open Sans", 10,))
	ReadMeQuit.grid(pady=(0, 10))

# File menu
filemenu = Menu(menubar, tearoff=0)
menubar.add_cascade(label="File", menu=filemenu)
filemenu.add_command(label="Exit", command=ExitApplication)

# Help menu
helpmenu = Menu(menubar, tearoff=0)
menubar.add_cascade(label="Help", menu=helpmenu)
helpmenu.add_command(label="License", command=OpenEULA)
helpmenu.add_command(label="User Guide", command=OpenUserMan)
helpmenu.add_command(label="Read Me", command=OpenReadME)

root.config(menu=menubar)

# ===================================================================
# 2.2 File I/O (IO)
# ===================================================================
s = ttk.Style()
s.configure('TLabelframe', background="#E8E8E8")
IOFrame = tk.LabelFrame(root, text="1. File I/O",fg="black", bg = "#E8E8E8", font =
					("Open Sans",10, 'bold'), height=100)
IOFrame.grid(column=0, row=1, padx=(7,5), pady=10)

# Labels
inputLabel = ttk.Label(IOFrame, text="Open File:",foreground="black", background = "#E8E8E8",
					   font = ("Open Sans",9)).grid(column=0, row=0,padx=(10,5),pady=(5, 2.5))
outputLabel = ttk.Label(IOFrame, text="Save As:",foreground="black", background  = "#E8E8E8",
						font = ("Open Sans",9)).grid(column=0, row=1,padx=(10,5),pady=(0,5))


Nsets=["All"]
Elsets=["All"]
def getSetNames(fileName):
	global Nsets
	global Elsets
	SetsDict=getSets(fileName)
	for curSet in SetsDict:
		if(curSet["Label"]=="Elset"):
			Elsets.append(curSet["Name"])
		if(curSet["Label"]=="Nset"):
			Nsets.append(curSet["Name"])


def OpenFile():
	name = askopenfilename(initialdir="MuPhiSiminputv2/MuPhiSiminput/content",filetypes =
						(("Abaqus File", "*.inp"),("All Files","*.*")), title = "Choose a file...")
	inputText.set(name)
	try:
		with open(name,'r') as UseFile:UseFile.read()
	except:
		print("No file exists")
	getSetNames(name)
	mainCode()

def SaveFile():
	name = asksaveasfilename(initialdir="MuPhiSiminputv2/MuPhiSiminput/content",
							filetypes =(("MuPhiSim File", "*.inp"),("All Files","*.*")),
							title = "Save file as...")
	outputText.set(name)

# Text-box Entry
inputText = tk.StringVar()
inputTextEntered = tk.Entry(IOFrame, width=14, textvariable=inputText, font = ("Open Sans",9))
inputTextEntered.grid(column = 1, row = 0,padx=(0,5),pady=(5,2.5))

outputText = tk.StringVar()
outputTextEntered = tk.Entry(IOFrame, width=14, textvariable=outputText,font = ("Open Sans",9))
outputTextEntered.grid(column = 1, row = 1,padx=(0,5),pady=(0,5))

# Buttons
def browseFile():
	browseButton.configure(text='Browsed')

def saveAs():
	saveButton.configure(text='Saved')

browseButton = tk.Button(IOFrame, text="Browse...", command=OpenFile,
						 foreground="black", background = "#d0d4d9", font = ("Open Sans",9))
browseButton.grid(column=2, row=0,padx=(0,10),pady=(5, 2.5))

saveButton = tk.Button(IOFrame, text="Save as...", command=SaveFile,
					   foreground="black", background = "#d0d4d9", font = ("Open Sans",9))
saveButton.grid(column=2, row=1,padx=(0,10),pady=(0,5))


def mainCode():
	global Nsets
	global Elsets
	# ===================================================================
	# 2.3 Region Set (RS)
	# ===================================================================
	# Build frame
	RSFrame = tk.LabelFrame(root, text="2. Region Set",fg="black", bg = "#E8E8E8", font = ("Open Sans",10, 'bold'), height=100)
	RSFrame.grid(column=1, row=1, padx=5, pady=10)

	# Labels
	SurfaceOrder = tk.Label(RSFrame, text="Surface element order:",foreground="black", background = "#E8E8E8",
						font = ("Open Sans",9)).grid(column=0, row=0,padx=(12,5),pady=(9, 10))
	FEMRSNodes = tk.Label(RSFrame, text="FEM Nodes Set:",foreground="black", background = "#E8E8E8",
						font = ("Open Sans",9)).grid(column=0, row=1,padx=(12,5),pady=(9, 10))
	FEMRSElements = tk.Label(RSFrame, text="FEM Elements Set:",foreground="black", background = "#E8E8E8",
						font = ("Open Sans",9)).grid(column=0, row=2,padx=(12,5),pady=(9, 10))

	MMRSNodes= tk.Label(RSFrame, text="MM Nodes Set:",foreground="black", background  = "#E8E8E8",
						font = ("Open Sans",9)).grid(column=0, row=3,padx=(12,5),pady=(3,11))


	MMRSElements = tk.Label(RSFrame, text="MM Elements Set:",foreground="black", background  = "#E8E8E8",
						font = ("Open Sans",9)).grid(column=0, row=4,padx=(12,5),pady=(3,11))
	NRCPU = tk.Label(RSFrame, text="Number of CPUs:",foreground="black", background  = "#E8E8E8",
						font = ("Open Sans",9)).grid(column=0, row=5,padx=(12,5),pady=(3,11))
	# Textbox Entry
	#self.Variables.append(None)
	#self.Variables_TxtEnt.append(ttk.Combobox(self.frame, width=11, state='readonly',values=Sols_Type, justify='center'))
	#self.Variables_TxtEnt[0].grid(column=0, row=self.position+1)
	NrofCPU=tk.StringVar()
	SurfaceOrder_TextEntered = ttk.Combobox(RSFrame, width=11, state='readonly',values=["Linear", "Quadratic"], justify='center')
	SurfaceOrder_TextEntered.grid(column = 1, row = 0,padx=(0,12),pady=(7,10))

	FEM_Nodes_TextEntered = ttk.Combobox(RSFrame, width=11, state='readonly',values=Nsets, justify='center')
	FEM_Nodes_TextEntered.grid(column = 1, row = 1,padx=(0,12),pady=(7,10))

	FEM_Elements_Text = tk.StringVar()
	FEM_Elements_TextEntered  = ttk.Combobox(RSFrame, width=11, state='readonly',values=Elsets, justify='center')
	FEM_Elements_TextEntered.grid(column = 1, row = 2,padx=(0,12),pady=(7,10))

	MM_Nodes_TextEntered  = ttk.Combobox(RSFrame, width=11, state='readonly',values=Nsets, justify='center')
	MM_Nodes_TextEntered.grid(column = 1, row = 3,padx=(0,12),pady=(0,7))

	MM_Elements_TextEntered  = ttk.Combobox(RSFrame, width=11, state='readonly',values=Elsets, justify='center')
	MM_Elements_TextEntered.grid(column = 1, row = 4,padx=(0,12),pady=(0,7))

	CPU_Elements_TextEntered  = tk.Entry(RSFrame, width=12, textvariable=NrofCPU,font = ("Open Sans",9))
	CPU_Elements_TextEntered.grid(column = 1, row = 5,padx=(0,12),pady=(0,7))



	###################################################################################
	CPUObj=None
	CPUCanvFrame=None
	CPUCanvas=None
	#SlvrScrollY=None
	CPUFrameInternal=None
	def CPUBut():
		nonlocal CPUCanvFrame
		if(CPUCanvFrame!=None):
			CPUCanvFrame.grid_forget()
		nonlocal CPUCanvas
		if(CPUCanvas!=None):
			CPUCanvas.grid_forget()
		nonlocal CPUFrameInternal
		if(CPUFrameInternal!=None):
			CPUFrameInternal.grid_forget()

		# Build frame
		CPUCanvFrame = Frame(RSFrame, borderwidth=0, bg="#E8E8E8")
		CPUCanvFrame.grid(column=0, row=6, columnspan=1)
		CPUCanvas = tk.Canvas(CPUCanvFrame)#, width=root.winfo_width(), height=60)
		#SlvrScrollY = tk.Scrollbar(SlvrCanvFrame, orient="vertical", command=SlvrCanvas.yview)
		CPUFrameInternal = tk.Frame(CPUCanvas)#,width=root.winfo_width())#, height = 60)
		NrofCPU = CPU_Elements_TextEntered.get()


		nonlocal CPUObj
		CPUObj=None
		if (NrofCPU.isdigit()):
			CPUObj=cpuBox(0, CPUFrameInternal, int(NrofCPU))
			#for i in range (0,int(NrofCPU)):
				# Widgets
				#CPUArray.append(choiceBoxSol(2*i,CPUFrameInternal))
			#SlvrScrollY.set(0.0, 0.0)
		else:
			IntegerErrorMssg = MessageBox.showerror("Error", "Please enter a whole number.")



		# put the frame in the canvas
		CPUCanvas.create_window(0, 0, window=CPUFrameInternal)
		# make sure everything is displayed before configuring the scrollregion
		CPUCanvas.update_idletasks()
		#SlvrCanvas.configure(scrollregion=SlvrCanvas.bbox('all'),
		#                 yscrollcommand=SlvrScrollY.set)
		CPUFrameInternal.pack(side="top", fill="both", expand=True)
		CPUCanvas.pack(fill='both', expand=True, side='left')
		root.update_idletasks()
		#SlvrScrollY.pack(fill='y', side='right')
	# Add Button
	NumCPUButton = Button(RSFrame, text="Confirm",foreground="black", command= CPUBut,
						background = "#d0d4d9", font = ("Open Sans",9))
	NumCPUButton.grid(column=2, row=5,padx=(0,10),pady=(0,5))
	###################################################################################

	# ===================================================================
	# 2.4 Solvers (Slvr)
	# ===================================================================
	# Build frame
	SlvrFrame = tk.LabelFrame(root, text="3. Solver Options",fg="black", bg = "#E8E8E8", font = ("Open Sans",10, 'bold') )#, width = root.winfo_width())
	SlvrFrame.grid(column=0, row=3, padx=20, pady=(0, 10), columnspan=2)

	# Build Button Frame
	SlvrButtonFrame = Frame(SlvrFrame, bg = "#E8E8E8")
	SlvrButtonFrame.grid(column=0, row=0, padx=(5,5), pady=(7, 8))

	# Labels
	SlvrLabel = tk.Label(SlvrButtonFrame, text="No. of Solvers:",foreground="black", background  = "#E8E8E8",
						font = ("Open Sans",9)).grid(column=0, row=0,padx=(12,5))


	# No of Solvers Entry
	NumSol = tk.StringVar()
	NumSlvrEntry = tk.Entry(SlvrButtonFrame, width=12, textvariable=NumSol,font = ("Open Sans",9))
	NumSlvrEntry.grid(column = 1, row = 0,padx=(0,5))
	#SlvrFrame.pack()
	root.update_idletasks()


	class cpuBox:
		def __init__(self, position, frame, cpunum):
			self.position = position
			self.frame=frame
			self.Variables =[]
			self.Labels= []
			self.Variables_TxtEnt=[]
			self.cpunum = cpunum
			if cpunum!=0:
				self.Labels.append(ttk.Label(self.frame, text="Name", justify='center'))
				self.Labels[0].grid(column=0, row=self.position)
				self.Labels.append(ttk.Label(self.frame, text="Nodes set", justify='center'))
				self.Labels[1].grid(column=0, row=self.position+1)
				self.Labels.append(ttk.Label(self.frame, text="Elements set", justify='center'))
				self.Labels[2].grid(column=0, row=self.position+2)
			for i in range(0,cpunum):
				self.Variables.append(tk.StringVar())
				self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=11, textvariable=self.Variables[-1],justify='center'))
				self.Variables_TxtEnt[-1].grid(column=i+1, row= self.position)

				self.Variables.append(ttk.Combobox(self.frame, width=11, state='readonly',values=Nsets, justify='center'))
				self.Variables[-1].grid(column=i+1, row= self.position+1)
				#self.Variables.append(tk.StringVar())
				self.Variables.append(ttk.Combobox(self.frame, width=11, state='readonly',values=Elsets, justify='center'))
				self.Variables[-1].grid(column=i+1, row= self.position+2)
				#self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=11, textvariable=self.Variables[i],justify='center'))
				#self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)



	class choiceBoxSol:

		def __init__(self, positon, frame):
			Sols_Type = ["IMPLICIT STATIC", "EXPLICIT", "IMPLICIT"]
			self.position=positon
			self.frame = frame
			self.Variables =[]
			self.Labels= []
			self.Variables_TxtEnt=[]

			self.Labels.append(ttk.Label(self.frame, text="Solver\nType", justify='center'))
			self.Labels.append(ttk.Label(self.frame, text="Time Scale\nFactor", justify='center'))
			self.Labels.append(ttk.Label(self.frame, text="Output File\nNumber", justify='center'))
			self.Labels.append(ttk.Label(self.frame, text="Internal\nVariables", justify='center'))
			self.Labels.append(ttk.Label(self.frame, text="Time", justify='center'))
			self.Labels.append(ttk.Label(self.frame, text="Displacement\nBC Sets", justify='center'))
			self.Labels.append(ttk.Label(self.frame, text="BC DOFs", justify='center'))
			self.Labels.append(ttk.Label(self.frame, text="Force output\nSets", justify='center'))
			self.Labels.append(ttk.Label(self.frame, text="Output\nfrequency", justify='center'))
			self.Labels.append(ttk.Label(self.frame, text="Force output\nDOFs", justify='center'))


			self.Variables.append(None)
			self.Variables_TxtEnt.append(ttk.Combobox(self.frame, width=11, state='readonly',values=Sols_Type, justify='center'))
			self.Variables_TxtEnt[0].grid(column=0, row=self.position+1)
			#self.Variables_TxtEnt[0].bind("<<ComboboxSelected>>",self.callback)
			for i,label in enumerate(self.Labels):
				label.grid(column=i, row=self.position)
				self.Variables.append(tk.StringVar())
				self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=11, textvariable=self.Variables[i],justify='center'))
				self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)



	SolversArray=[]
	SlvrCanvFrame=None
	SlvrCanvas=None
	#SlvrScrollY=None
	SlvrFrameInternal=None
	def SolverBut():
		nonlocal SlvrCanvFrame
		if(SlvrCanvFrame!=None):
			SlvrCanvFrame.grid_forget()
		nonlocal SlvrCanvas
		if(SlvrCanvas!=None):
			SlvrCanvas.grid_forget()
		nonlocal SlvrFrameInternal
		if(SlvrFrameInternal!=None):
			SlvrFrameInternal.grid_forget()

		# Build frame
		SlvrCanvFrame = Frame(SlvrFrame, borderwidth=0, bg="#E8E8E8")
		SlvrCanvFrame.grid(column=0, row=1, columnspan=1)
		SlvrCanvas = tk.Canvas(SlvrCanvFrame)#, width=root.winfo_width(), height=60)
		#SlvrScrollY = tk.Scrollbar(SlvrCanvFrame, orient="vertical", command=SlvrCanvas.yview)
		SlvrFrameInternal = tk.Frame(SlvrCanvas)#,width=root.winfo_width())#, height = 60)
		NumSol = NumSlvrEntry.get()


		nonlocal SolversArray
		SolversArray=[]
		if (NumSol.isdigit()):
			for i in range (0,int(NumSol)):
				# Widgets
				SolversArray.append(choiceBoxSol(2*i,SlvrFrameInternal))
			#SlvrScrollY.set(0.0, 0.0)
		else:
			IntegerErrorMssg = MessageBox.showerror("Error", "Please enter a whole number.")



		# put the frame in the canvas
		SlvrCanvas.create_window(0, 0, window=SlvrFrameInternal)
		# make sure everything is displayed before configuring the scrollregion
		SlvrCanvas.update_idletasks()
		#SlvrCanvas.configure(scrollregion=SlvrCanvas.bbox('all'),
		#                 yscrollcommand=SlvrScrollY.set)
		SlvrFrameInternal.pack(side="top", fill="both", expand=True)
		SlvrCanvas.pack(fill='both', expand=True, side='left')
		#SlvrScrollY.pack(fill='y', side='right')

	# Add Button
	NumSlvrButton = Button(SlvrButtonFrame, text="Confirm",foreground="black", command= SolverBut,
						background = "#d0d4d9", font = ("Open Sans",9))
	NumSlvrButton.grid(column=2, row=0,padx=(0,10),pady=(0,5))


	# ===================================================================
	# 2.5 Material Laws (ML)
	# ===================================================================
	# Build frame
	MLFrame = tk.LabelFrame(root, text="4. Material Laws",fg="black", bg = "#E8E8E8", font = ("Open Sans",10, 'bold'), width = 800)
	MLFrame.grid(column=0, row=4, padx=20, pady=(0, 10),columnspan=2)

	# Build Button Frame
	MLButtonFrame = Frame(MLFrame, bg = "#E8E8E8")
	MLButtonFrame.grid(column=0, row=0, padx=(5,5), pady=(7, 8))
	# Labels
	MLLabel = tk.Label(MLButtonFrame, text="No. of Material Laws:",foreground="black", background  = "#E8E8E8",
						font = ("Open Sans",9)).grid(column=0, row=0,padx=(12,5))
	# No of Material Laws Entry
	NumML = tk.StringVar()
	NumMLEntry = tk.Entry(MLButtonFrame, width=12, textvariable=NumML,font = ("Open Sans",9))
	NumMLEntry.grid(column = 1, row = 0,padx=(0,5))
	root.update_idletasks()


	class choiceBoxMat:

		def __init__(self, positon, frame):
			MLTbl_Type = ["HyperElastic, St-Venant-Kirchhoff", "HyperElastic, Neo-Hookean", "ViscoElastic", "FiberGrowth",
									"AreaGrowth", "ViscoGrowth", "ViscoElastic",
									"ViscoMechAxialGrowth", "FHNelec", "ContractAxialGrowth", 										"Temperature", "Linear-Elastic","HyperElastic, St-Venant-Kirchhoff-Thermomechanics", "3DprintingTemperature","3DPrinting-St-Venant-Kirchhoff-Thermoelasticity" ] ###To be modified if constitutive model needs to be added
			self.position=positon
			self.frame = frame
			self.Variables =[]
			self.Labels= []
			self.Variables_TxtEnt=[]

			self.Labels.append(ttk.Label(self.frame, text="Material\nName", justify='center'))
			self.Labels.append(ttk.Label(self.frame, text="Material\nSet", justify='center'))
			self.Labels.append(ttk.Label(self.frame, text="Material\nDensity", justify='center'))
			self.Labels.append(ttk.Label(self.frame, text="Material\nModel", justify='center'))
			for i,label in enumerate(self.Labels):
				label.grid(column=i, row=self.position)


			self.Variables.append(tk.StringVar())
			self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=11, textvariable=self.Variables[0],justify='center'))
			self.Variables_TxtEnt[0].grid(column=0, row=self.position+1)

			self.Variables.append(None)#consistency in length
			self.Variables_TxtEnt.append(ttk.Combobox(self.frame, width=11, state='readonly',values=Elsets, justify='center'))
			self.Variables_TxtEnt[1].grid(column=1, row=self.position+1)

			self.Variables.append(tk.StringVar())
			self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=11, textvariable=self.Variables[2],justify='center'))
			self.Variables_TxtEnt[2].grid(column=2, row=self.position+1)


			self.Variables.append(None)#consistency in length
			self.Variables_TxtEnt.append(ttk.Combobox(self.frame, width=11, state='readonly',values=MLTbl_Type, justify='center'))
			self.Variables_TxtEnt[3].grid(column=3, row=self.position+1)
			self.Variables_TxtEnt[3].bind("<<ComboboxSelected>>",self.callback)




		def callback(self, event=None):
			for i in range(4, len(self.Variables_TxtEnt)):#remove and empty all variables after 4
				self.Variables_TxtEnt[4].destroy()
				self.Variables_TxtEnt.pop(4)
				self.Variables.pop(4)
				self.Labels[4].destroy()
				self.Labels.pop(4)
			columntrack = 4
			if(self.Variables_TxtEnt[3].get() == "HyperElastic, Neo-Hookean"):
				self.Labels.append(ttk.Label(self.frame, text="Young's Modulus, Poisson Ratio", justify='center'))
				for i in range(columntrack, columntrack+1):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=26, textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
			if(self.Variables_TxtEnt[3].get() == "HyperElastic, St-Venant-Kirchhoff"):
				self.Labels.append(ttk.Label(self.frame, text="Young's Modulus, Poisson Ratio, Stress State (optional)", justify='center'))
				for i in range(columntrack, columntrack+1):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=36, textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
			if( self.Variables_TxtEnt[3].get() == "AreaGrowth" or self.Variables_TxtEnt[3].get() == "FiberGrowth" or self.Variables_TxtEnt[3].get() == "Iso Morpho Growth"):
				self.Labels.append(ttk.Label(self.frame, text="Young's Modulus, Poisson Ratio, Growth Multiplier", justify='center'))
				for i in range(columntrack, columntrack+1):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=36, textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
			if(self.Variables_TxtEnt[3].get() == "FHNelec"): 
				self.Labels.append(ttk.Label(self.frame, text="Nx, Ny, Nz, Isotropic Component of D, Anisotropic Component of D, Activation Parameter a, Activation Parameter b, Time Scale Difference, Equilibrium Potential, Stimulus Current", justify='center'))
				for i in range(columntrack, columntrack+1):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=140, textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
			if(self.Variables_TxtEnt[3].get() == "Temperature"): 
				self.Labels.append(ttk.Label(self.frame, text="Specific heat (C), C_tempDepend, Conductivity (k), k_tempDepend, Temp Initial, Temp Ref1, Temp Ref2", justify='center'))
				for i in range(columntrack, columntrack+1):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=98, textvariable=self.Variables[i],justify='center'))		
				self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)					
			if(self.Variables_TxtEnt[3].get() == "ViscoElastic"):
				self.Labels.append(ttk.Label(self.frame, text="Bulk Modulus", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Number Of\nBranches", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Shear Modulii\nμ∞, μ1,..", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Viscosities\nη1, η2,..", justify='center'))
				for i in range(columntrack, columntrack+4):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=11, textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
			if(self.Variables_TxtEnt[3].get() == "ViscoGrowth"):
				self.Labels.append(ttk.Label(self.frame, text="Bulk Modulus", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Number Of\nBranches", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Shear Modulii\nμ∞, μ1,..", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Viscosities\nη1, η2,..", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Growth Multiplier", justify='center'))
				for i in range(columntrack, columntrack+5):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=11, textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
			if(self.Variables_TxtEnt[3].get() == "ViscoMechAxialGrowth"):
				self.Labels.append(ttk.Label(self.frame, text="Bulk Modulus", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Number Of\nBranches", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Shear Modulii\nμ∞, μ1,..", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Viscosities\nη1, η2,..", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Growth Multiplier, Polymerisation Rate, Depolymerisation Rate, Cytoskeletal Density, Growth Direction", justify='center'))
				lenght=lambda x: 80 if x==columntrack+5-1 else 11
				for i in range(columntrack, columntrack+5):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=lenght(i), textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
			if(self.Variables_TxtEnt[3].get() == "ContractAxialGrowth"):
				self.Labels.append(ttk.Label(self.frame, text="Bulk Modulus", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Number Of\nBranches", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Shear Modulii\nμ∞, μ1,..", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Viscosities\nη1, η2,..", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Growth Multiplier, Polymerisation Rate, Depolymerisation Rate, Cytoskeletal Density, Contractility Rate, Contractility Stress Threshold, Contractility Growth Direction", justify='center'))
				lenght=lambda x: 120 if x==columntrack+5-1 else 11
				for i in range(columntrack, columntrack+5):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=lenght(i), textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
			if(self.Variables_TxtEnt[3].get() == "Linear-Elastic"):
				self.Labels.append(ttk.Label(self.frame, text="Young's Modulus, Poisson Ratio, Stress State (optional)", justify='center'))
				for i in range(columntrack, columntrack+1):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=36, textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
			if(self.Variables_TxtEnt[3].get() == "HyperElastic, St-Venant-Kirchhoff-Thermomechanics"): 
				self.Labels.append(ttk.Label(self.frame, text="Young's Modulus, Poisson Ratio, Thermal expansion, Initial Temp, Stress State (optional)", justify='center'))
				for i in range(columntrack, columntrack+1):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=98, textvariable=self.Variables[i],justify='center'))		
				self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
			if(self.Variables_TxtEnt[3].get() == "3DprintingTemperature"): 
				self.Labels.append(ttk.Label(self.frame, text="C Pow-Sol , C L, C1 Pow , C1 Sol , L1 Pow-L , L2 Pow-L , L3 Pow-L , L1 Sol-L , L2 Sol-L ,L3 Sol-L , k1 Pow, k2 Pow , k L, k Sol , k1 Pow-L , k2 Pow-L , k1 Sol-L , Temp ref , Temp Pow, Temp Sol, Temp L", justify='center'))
				for i in range(columntrack, columntrack+1):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=140, textvariable=self.Variables[i],justify='center'))		
				self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
			if(self.Variables_TxtEnt[3].get() == "3DPrinting-St-Venant-Kirchhoff-Thermoelasticity"): 
				self.Labels.append(ttk.Label(self.frame, text="\u03C1 Liq, \u03C1 Pow-L, E Pow , E L , E Sol, E1 Pow-L , E1 Sol-L , E2 Sol-L , \u03C5 L, \u03C5 Sol, \u03C5 Pow-L, \u03C5 Sol-L, Log growth rate Sol-L, Temp midpoint, \u03B1 exp Sol-L, \u03B1 exp Sol, Temp ref , Temp ini, Temp L, Temp Sol-L", justify='center'))
				for i in range(columntrack, columntrack+1):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=140, textvariable=self.Variables[i],justify='center'))		
				self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)

	MaterialsArray=[]

	MLCanvFrame=None
	MLCanvas=None
	MLFrameInternal=None

	def MLBut():
		nonlocal MLCanvFrame
		if(MLCanvFrame!=None):
			MLCanvFrame.grid_forget()
		nonlocal MLCanvas
		if(MLCanvas!=None):
			MLCanvas.grid_forget()
		nonlocal MLFrameInternal
		if(MLFrameInternal!=None):
			MLFrameInternal.grid_forget()

		# Build frame
		MLCanvFrame = Frame(MLFrame, borderwidth=0, bg="#E8E8E8")
		MLCanvFrame.grid(column=0, row=1, columnspan=1)
		MLCanvas = tk.Canvas(MLCanvFrame, width=root.winfo_width(), height=60)
		#MLScrollY = tk.Scrollbar(MLCanvFrame, orient="vertical", command=MLCanvas.yview)
		MLFrameInternal = tk.Frame(MLCanvas,width=root.winfo_width())#, height = 500)

		NumML = NumMLEntry.get()
		nonlocal MaterialsArray
		MaterialsArray=[]
		if (NumML.isdigit()):



			for i in range (0,int(NumML)):
				# Widgets
				MaterialsArray.append(choiceBoxMat(2*i,MLFrameInternal))
			#MLScrollY.set(0.0, 0.0)
		else:
			IntegerErrorMssg = MessageBox.showerror("Error", "Please enter a whole number.")

		# put the frame in the canvas
		MLCanvas.create_window(0, 0, window=MLFrameInternal)

		# make sure everything is displayed before configuring the scroll region
		MLCanvas.update_idletasks()
		#MLCanvas.configure(scrollregion=MLCanvas.bbox('all'),
		#                 yscrollcommand=MLScrollY.set)
		MLFrameInternal.pack(side="top", fill="both", expand=True)
		MLCanvas.pack(fill='both', expand=True, side='left')
		#MLScrollY.pack(fill='both', side='right')

	# Add Button
	NumMLButton = Button(MLButtonFrame, text="Confirm",foreground="black", command= MLBut,
						background = "#d0d4d9", font = ("Open Sans",9))
	NumMLButton.grid(column=2, row=0,padx=(0,10),pady=(0,5))

	# ===================================================================
	# 2.5 NBC
	# ===================================================================
	# Build frame
	BCFrame = tk.LabelFrame(root, text="5. Boundary Conditions",fg="black", bg = "#E8E8E8", font = ("Open Sans",10, 'bold'), width = 800)
	BCFrame.grid(column=0, row=5, padx=20, pady=(0, 10), columnspan=2)

	# Build Force & Pressure Frame
	FBCFrame = Frame(BCFrame,bg = "#E8E8E8", width = 395)
	FBCFrame.grid(column=0,row=0)


	# Build Button Frame
	FBCButtonFrame = Frame(FBCFrame, bg = "#E8E8E8")
	FBCButtonFrame.grid(column=0, row=0, padx=(5,5), pady=(7, 8))

	# Labels
	FBCLabel = tk.Label(FBCButtonFrame, text="No. of Force\nBoundary Conditions:",foreground="black", background  = "#E8E8E8",
						font = ("Open Sans",9)).grid(column=0, row=0,padx=(12,5))

	# No of Material Laws Entry
	NumFBC = tk.StringVar()
	NumFBCEntry = tk.Entry(FBCButtonFrame, width=12, textvariable=NumFBC,font = ("Open Sans",9))
	NumFBCEntry.grid(column = 1, row = 0,padx=(0,5))

	class choiceBoxBC:

		def __init__(self, positon, frame):
			MLTbl_Type = ["Pressure ramp", "Pressure inst", "Force", "Current inst", "Heat inst", "Volumetric heat flux"] ###Place where BC have to be added
			self.position=positon
			self.frame = frame
			self.Variables =[]
			self.Labels= []
			self.Variables_TxtEnt=[]

			self.Labels.append(ttk.Label(self.frame, text="BC Type", justify='center'))
			for i,label in enumerate(self.Labels):
				label.grid(column=i, row=self.position)

			self.Variables.append(None)#consistency in length
			self.Variables_TxtEnt.append(ttk.Combobox(self.frame, width=11, state='readonly',values=MLTbl_Type, justify='center'))
			self.Variables_TxtEnt[0].grid(column=0, row=self.position+1)
			self.Variables_TxtEnt[0].bind("<<ComboboxSelected>>",self.callback)




		def callback(self, event=None):
			for i in range(1, len(self.Variables_TxtEnt)):#remove and empty all variables after 4
				self.Variables_TxtEnt[1].destroy()
				self.Variables_TxtEnt.pop(1)
				self.Variables.pop(1)
				self.Labels[1].destroy()
				self.Labels.pop(1)
			columntrack = 1
			if(self.Variables_TxtEnt[0].get() == "Pressure inst" or self.Variables_TxtEnt[0].get() == "Pressure ramp"):
				self.Labels.append(ttk.Label(self.frame, text="Elements\nSet", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Start\nTime", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="End\nTime", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Pressure\nValue", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Surface", justify='center'))
				self.Variables_TxtEnt.append(ttk.Combobox(self.frame, width=11, state='readonly',values=Elsets, justify='center'))
				self.Variables_TxtEnt[columntrack].grid(column=columntrack, row=self.position+1)
				self.Variables.append(None)
				self.Labels[columntrack].grid(column=columntrack, row=self.position)
				columntrack += 1
				for i in range(columntrack, columntrack+4):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=11, textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)


			if( self.Variables_TxtEnt[0].get() == "Force"):
				self.Labels.append(ttk.Label(self.frame, text="Nodes\nSet", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Start\nTime", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="End\nTime", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Force\nVector", justify='center'))
				self.Variables_TxtEnt.append(ttk.Combobox(self.frame, width=11, state='readonly',values=Nsets, justify='center'))
				self.Variables_TxtEnt[columntrack].grid(column=columntrack, row=self.position+1)
				self.Variables.append(None)
				self.Labels[columntrack].grid(column=columntrack, row=self.position)
				columntrack += 1
				for i in range(columntrack, columntrack+3):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=11, textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
					
					
			if( self.Variables_TxtEnt[0].get() == "Current inst"):  
				self.Labels.append(ttk.Label(self.frame, text="Elements\nSet", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Start\nTime", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="End\nTime", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Current\nValue", justify='center'))
				self.Variables_TxtEnt.append(ttk.Combobox(self.frame, width=11, state='readonly',values=Elsets, justify='center'))
				self.Variables_TxtEnt[columntrack].grid(column=columntrack, row=self.position+1)
				self.Variables.append(None)
				self.Labels[columntrack].grid(column=columntrack, row=self.position)
				columntrack += 1
				for i in range(columntrack, columntrack+3):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=11, textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
					
			if(self.Variables_TxtEnt[0].get() == "Heat inst"): 	
				self.Labels.append(ttk.Label(self.frame, text="Elements\nSet", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Start\nTime", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="End\nTime", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Surface heat\nflux Vector", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Surface", justify='center'))
				self.Variables_TxtEnt.append(ttk.Combobox(self.frame, width=11, state='readonly',values=Elsets, justify='center'))
				self.Variables_TxtEnt[columntrack].grid(column=columntrack, row=self.position+1)
				self.Variables.append(None)
				self.Labels[columntrack].grid(column=columntrack, row=self.position)
				columntrack += 1
				for i in range(columntrack, columntrack+4):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=11, textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)
                    
			if(self.Variables_TxtEnt[0].get() == "Volumetric heat flux"): 	

				self.Labels.append(ttk.Label(self.frame, text="Elements\nSet", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Start\nTime", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="End\nTime", justify='center'))
				self.Labels.append(ttk.Label(self.frame, text="Heat flux\nvalue", justify='center'))
				self.Variables_TxtEnt.append(ttk.Combobox(self.frame, width=11, state='readonly',values=Elsets, justify='center'))
				self.Variables_TxtEnt[columntrack].grid(column=columntrack, row=self.position+1)
				self.Variables.append(None)
				self.Labels[columntrack].grid(column=columntrack, row=self.position)
				columntrack += 1
				for i in range(columntrack, columntrack+3):
					self.Labels[i].grid(column=i, row=self.position)
					self.Variables.append(tk.StringVar())
					self.Variables_TxtEnt.append(ttk.Entry(self.frame, width=11, textvariable=self.Variables[i],justify='center'))
					self.Variables_TxtEnt[i].grid(column=i, row= self.position+1)



	FBCCanvFrame = None
	FBCCanvas = None
	FBCFrameInternal = None
	BCsArray=[]
	def FBCBut():
		nonlocal FBCCanvFrame
		if(FBCCanvFrame!=None):
			FBCCanvFrame.grid_forget()
		nonlocal FBCCanvas
		if(FBCCanvas!=None):
			FBCCanvas.grid_forget()
		nonlocal FBCFrameInternal
		if(FBCFrameInternal!=None):
			FBCFrameInternal.grid_forget()


		NumFBC = NumFBCEntry.get()
		if (NumFBC.isdigit()):
			# Build frame
			FBCCanvFrame = Frame(FBCFrame, borderwidth=0, bg="#E8E8E8")
			FBCCanvFrame.grid(column=0, row=1, columnspan=1)
			FBCCanvas = tk.Canvas(FBCCanvFrame, width=800, height=60)
			#FBCScrollY = tk.Scrollbar(FBCCanvFrame, orient="vertical", command=FBCCanvas.yview)
			FBCFrameInternal = tk.Frame(FBCCanvas,width=800)#, height = 500)

			nonlocal BCsArray
			BCsArray=[]

			for i in range (0,int(NumFBC)):

				BCsArray.append(choiceBoxBC(2*i,FBCFrameInternal))

		else:
			IntegerErrorMssg = MessageBox.showerror("Error", "Please enter a whole number.")


		# put the frame in the canvas
		FBCCanvas.create_window(0, 0, window=FBCFrameInternal)
		# make sure everything is displayed before configuring the scrollregion
		FBCCanvas.update_idletasks()
		#FBCCanvas.configure(scrollregion=FBCCanvas.bbox('all'),
		#                 yscrollcommand=FBCScrollY.set)

		FBCFrameInternal.pack(fill='both', expand=True, side='left')
		FBCCanvas.pack(fill='both', expand=True, side='left')
		#FBCScrollY.pack(fill='both', side='right')


	# Add Button
	NumFBCButton = Button(FBCButtonFrame, text="Confirm",foreground="black", command= FBCBut,
						background = "#d0d4d9", font = ("Open Sans",9))
	NumFBCButton.grid(column=2, row=0,padx=(0,10),pady=(0,5))

	# ===================================================================
	# 2.6 Execute Fn
	# ===================================================================
	def Execute():

		Var_inputText = inputText.get()
		if not inputText.get():
			EmptyInputErrorMssg = MessageBox.showerror("Error", "Please choose a file to open.")

		Var_outputText = outputText.get()
		if not outputText.get():
			EmptyInputErrorMssg = MessageBox.showerror("Error", "Please choose name to save as.")


		Var_SurfaceOrder_TextEntered = SurfaceOrder_TextEntered.get()
		Var_FEM_Nodes_TextEntered = FEM_Nodes_TextEntered.get()
		Var_FEM_Elements_TextEntered = FEM_Elements_TextEntered.get()
		Var_MM_Nodes_TextEntered = MM_Nodes_TextEntered.get()
		Var_MM_Elements_TextEntered = MM_Elements_TextEntered.get()


		Var_SlvrTbl_SlvrType_TxtEnt=[]
		Var_SlvrTbl_TimeScaleFctr_TxtEnt=[]
		Var_SlvrTbl_OutFileNo_TxtEnt=[]
		Var_SlvrTbl_IntVar_TxtEnt=[]
		Var_SlvrTbl_Time_TxtEnt=[]
		Var_SlvrTbl_Disp_TxtEnt=[]
		Var_SlvrTbl_DOF_TxtEnt=[]
		Var_SlvrTbl_Forces_TxtEnt=[]
		Var_SlvrTbl_DOFForces_TxtEnt=[]
		Var_SlvrTbl_FREQ_TxtEnt=[]


		# Solver Variables
		if not NumSlvrEntry.get():
			EmptyInputErrorMssg = MessageBox.showerror("Error", "Please insert a value for 'Number of Solvers'")
		for solver in SolversArray:
			Var_SlvrTbl_SlvrType_TxtEnt.append(solver.Variables_TxtEnt[0].get())
			Var_SlvrTbl_TimeScaleFctr_TxtEnt.append(solver.Variables_TxtEnt[1].get())
			Var_SlvrTbl_OutFileNo_TxtEnt.append(solver.Variables_TxtEnt[2].get())
			Var_SlvrTbl_IntVar_TxtEnt.append(solver.Variables_TxtEnt[3].get())
			Var_SlvrTbl_Time_TxtEnt.append(solver.Variables_TxtEnt[4].get())
			Var_SlvrTbl_Disp_TxtEnt.append(solver.Variables_TxtEnt[5].get())
			Var_SlvrTbl_DOF_TxtEnt.append(solver.Variables_TxtEnt[6].get())
			Var_SlvrTbl_Forces_TxtEnt.append(solver.Variables_TxtEnt[7].get())
			Var_SlvrTbl_FREQ_TxtEnt.append(solver.Variables_TxtEnt[8].get())
			Var_SlvrTbl_DOFForces_TxtEnt.append(solver.Variables_TxtEnt[9].get())



		# Variables_ Material Laws
		if not NumMLEntry.get():
			EmptyInputErrorMssg = MessageBox.showerror("Error", "Please insert a value for 'Number of Material Laws'")
		Var_MLTbl_Name_TxtEnt=[]
		Var_MLTbl_Set_TxtEnt=[]
		Var_MLTbl_Density_TxtEnt=[]
		Var_MLTbl_Model_TxtEnt=[]
		Var_MLTbl_Parameters_TxtEnt=[]
		for material in MaterialsArray:
			Var_MLTbl_Name_TxtEnt.append(material.Variables_TxtEnt[0].get())
			Var_MLTbl_Set_TxtEnt.append(material.Variables_TxtEnt[1].get())
			Var_MLTbl_Density_TxtEnt.append(material.Variables_TxtEnt[2].get())
			Var_MLTbl_Model_TxtEnt.append(material.Variables_TxtEnt[3].get())
			tempArrPar=[]
			for i in range(4, len(material.Variables_TxtEnt)):
				tempArrPar.append(material.Variables_TxtEnt[i].get())
			Var_MLTbl_Parameters_TxtEnt.append(tempArrPar)


		# Boundary Conditions
		Var_FBCTbl_ForceType_TxtEnt =[]
		Var_FBCTbl_ForceSet_TxtEnt =[]
		Var_FBCTbl_TimeL_TxtEnt = []
		Var_FBCTbl_TimeU_TxtEnt = []
		Var_FBCTbl_Value_TxtEnt = []
		Var_FBCTbl_Surface_TxtEnt = []
		for BC in BCsArray:
			Var_FBCTbl_ForceType_TxtEnt.append(BC.Variables_TxtEnt[0].get())
			Var_FBCTbl_ForceSet_TxtEnt.append(BC.Variables_TxtEnt[1].get())
			Var_FBCTbl_TimeL_TxtEnt.append(BC.Variables_TxtEnt[2].get())
			Var_FBCTbl_TimeU_TxtEnt.append(BC.Variables_TxtEnt[3].get())
			Var_FBCTbl_Value_TxtEnt.append(BC.Variables_TxtEnt[4].get())
			if(BC.Variables_TxtEnt[0].get()!="Force" and BC.Variables_TxtEnt[0].get()!="Current inst" and BC.Variables_TxtEnt[0].get()!="Volumetric heat flux" ): 
				Var_FBCTbl_Surface_TxtEnt.append(BC.Variables_TxtEnt[5].get())
			#print(Var_FBCTbl_ForceSet_TxtEnt[0])
		Var_CPUs_TxtEnt=[]
		if CPUObj:
			for CPU in CPUObj.Variables:
				Var_CPUs_TxtEnt.append(CPU.get())

		writeFileOut(Var_inputText,Var_outputText,Var_FEM_Nodes_TextEntered,Var_FEM_Elements_TextEntered,
		                Var_MM_Nodes_TextEntered,Var_MM_Elements_TextEntered,NumSol.get(),Var_SlvrTbl_SlvrType_TxtEnt,
		                Var_SlvrTbl_TimeScaleFctr_TxtEnt,Var_SlvrTbl_OutFileNo_TxtEnt,Var_SlvrTbl_IntVar_TxtEnt,
		                Var_SlvrTbl_Time_TxtEnt,Var_SlvrTbl_Disp_TxtEnt,Var_SlvrTbl_DOF_TxtEnt,
						NumFBCEntry.get(),Var_FBCTbl_ForceType_TxtEnt,Var_FBCTbl_ForceSet_TxtEnt,Var_FBCTbl_TimeL_TxtEnt, Var_FBCTbl_TimeU_TxtEnt,
						Var_FBCTbl_Value_TxtEnt,Var_FBCTbl_Surface_TxtEnt,NumMLEntry.get(),Var_MLTbl_Name_TxtEnt,Var_MLTbl_Set_TxtEnt,
						Var_MLTbl_Density_TxtEnt,Var_MLTbl_Model_TxtEnt,Var_MLTbl_Parameters_TxtEnt,Var_SlvrTbl_Forces_TxtEnt,Var_SlvrTbl_DOFForces_TxtEnt,Var_SlvrTbl_FREQ_TxtEnt,
						Var_CPUs_TxtEnt,Var_SurfaceOrder_TextEntered)


	# ===================================================================
	# 2.6 Execute/Exit Buttons (CE)
	# ===================================================================
	# Build frame
	CEFrame = Frame(root, borderwidth=0, bg="#E8E8E8")
	CEFrame.grid(column=0, row=1, padx=20, pady=(0, 20), sticky=tk.E, columnspan=2) 

	# Create Buttons
	ExecuteButton = ttk.Button(CEFrame, text="Execute", command=Execute)
	ExecuteButton.grid(column=0, row=0, sticky=tk.E)

	ExitButton = ttk.Button(CEFrame, text="Exit", command=ExitApplication)
	ExitButton.grid(column=1, row=0, sticky=tk.E)

	# ===================================================================
	# Status
	# ===================================================================



	# ===================================================================
	# Initiate GUI
	# ===================================================================
root.mainloop()
