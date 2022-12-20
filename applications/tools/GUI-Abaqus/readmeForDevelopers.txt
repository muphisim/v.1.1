Guidance for developers by Khariton Gorbunov

The GUI code consists of two scripts, MuPhiSiminput.py and parser.py. MuPhiSiminput.py handles the graphics and the user input, while parser will take data from MuPhiSiminput and write it out accordingly.

In MuPhiSiminput, every box is an object (class), and will have optional parameters, whose structure should be quite obvious. When the user presses execute, the parser script is called, taking all user defined variables. Every parameter for every class will be added to an array (for instance if a user wishes two solvers, and specifies type "IMPLICIT" and "EXPLICIT" respectivly, then there will be an array called Slvr_TBL_type=["IMPLICIT","EXPLICIT"]), to be used for processing in a loop that goes through number of items of each class (in this example 2 solvers).

In parser, the code loops over number of solvers, then over Neumann boundary conditions and then materials. Inside of these, it loops over the relevant parameters, printing the relevant string outputs. Unless otherwise specified, every item in the array will be printed out as a string on a separate line. However, some will have special functions, such as "writeBoundary" which handles the Boundary output. The sets are stored as dictionaries, and are either labled Nset or Elset (nodes set or elements set).

