#class to define configure file
#Author(s): Van Dung NGUYEN, 2022
#


"""In the following, the following determination is considered a group of nodes or of elements:
    - "All": means all elements available in problem
    - an integer: a particular node or a particular element depending on the context
    - a tuple of two integers (dim, phy): correspond to the dimension and the identification numbers of the groups. This tuple must be unique for each group
"""

class Configure :
    def getFEMConfigure(self):
        """ FEM options are provided within a dictionary as
        self.FEM = {
                       #support All or [(dim,phys), ...] or nothing (all are MM elements)
                       "Support": tuple of location or "All",
                    }
       
        """
        raise NotImplementedError('The method getFEMConfigure() must be defined in Configure!')
        
    def getMMConfigure(self):
        """ MM options are provided within a dictionary as
        self.MM = {
                       #support All or [(dim,phys), ...] or nothing (all are FEM elements)
                       "Support": tuple of location or "All",
                    }
       
        """
        raise NotImplementedError('The method getMMConfigure() must be defined in Configure!')
    
    def getPartionsByElementGroups(self):
        """ partition options are provided within a dictionary as
        self.partitions = {
                      0: [list of element, list of nodes]
                      ...
                    }
       
        """
        return None
    
    def getSolversConfigure(self):
        """get solver options need to be provided with a tuple of dictionaries, each dictionary corresponds to a solver, 
        The structure is defined as follows
        
        self.solvers = ({
                            #solver type
                            "Solver": "IMPLICIT STATIC",
                            #scale factor to estimate time sep for solver
                            "Scale Factor": 1e10,
                            #number of steps, "Scale Factor" is not used if Number Of Steps is set
                            "Number Of Steps": 70,
                            #output specification, each component of the list corresponds a line in the oxfemm input file
                            "Outputs": ((1000,), ("All" , "max"),),
                            # start and end time
                            "Time": (0.,10.),
                            #FSI with locations it is considered
                            "FSI": ((1,2),), 
                            "DispBC":(
                                        # each bc contains (type, location(dim, physical) list of prescried components and prescribed values
                                        ("DISPLACEMENT RAMP",(2,43), ((0,0.),(1,0.),(2,0.))), 
                                        ("DISPLACEMENT RAMP",(2,51), ((0,0.1),(1,0.15),(2,0.2))),
                                        #("DISPLACEMENT RAMP",(2,43), ((3,273.),)),
                                        #("DISPLACEMENT RAMP",(2,51), ((3,500.),)),
                                     ), 
                             "FSI": ((dim, phy), ...), # location
                             "Extraction":(
                                            # each extraction contains 
                                            # (1) type=NODE or ELEMENT, 
                                            # (2) an user-define string as the first part of the file saved,
                                            # (3) the quantity to save, for NODE, Unknow, ExternalForce or InternalForce should be used
                                            #                       for ELEMENT, FXX, ..., FZZ, PXX, ...,PZZ, DEFO_ENERGY, EXTERNAL_ENERGY, KIN_ENERGY 
                                            #                       should be used
                                            # (4) operation = one of the following operations: Mean, Max, Min, Sum
                                            # (5) location
                                            # if type=NODE, the field (6) must be filled to define the component to save
                                            ("NODE", "MyExtract", "Unknown",Mean,  5, 1), 
                                            # group of nodes defined by (1, 3) and component
                                            ("NODE", "MyExtract", "Unknown",Mean, (1, 3), 1),
                                            # extraction for a group (1,3)
                                            "ELEMENT", "MyExtract", "DEFO_ENERGY",Mean, (1, 3)),
                                            # extract over the whole domain 
                                            "ELEMENT", "MyExtract", "DEFO_ENERGY",Mean, "All"),
                                         )
                            },
                        },
                        {
                            next solver
                        }
                        )
          
        """
        raise NotImplementedError('The method getSolversConfigure() must be defined in Configure!')
        
    def getNeumannBCsConfigure(self):
        """Neumann BCs are defined by a dictionary as
        self.NeumannBCs = {
                Type: (
                        (BC, location, time interval, values),
                      )
            }
            
            Type : "PressureBC" for pressure BC, "TractionBC" for traction BC, "Flux ExtraDof" for extra Dof, #"GRAVITY" for gravity
            BC: a BC keyword depending on Type as "PRESSURE RAMP", "TRACTION RAMP",  "RADIATION", ...
            location: whre the BC is applied by a pair (dim, phy)
            time interval: a pair specify start and end time
            value: depending on case, pressure = a float, traction = a tuple of 3 values, 
          It is noted that for CONVECTION, and RADIATION,
            # type, location, time interval, h and sink temperature
            #("CONVECTION", (2,52), (0., 0.5), (800.,0)),
            # type, location, time interval, A and sink temperature
            #("RADIATION", (2,52), (0., 0.5), (500.,0)),
        """
        raise NotImplementedError('The method getNeumannBCsConfigure() must be defined in Configure!')
        
        
    def getMaterialsConfigure(self):
        """Material parameters must be provided by a tuple of dictionaries as
         self.materials = ({
                #material name
                "Name": "Mat1",
                #density
                "Density": 1.,
                #material type
                "Type": "Linear-Elastic",
                # parameter of the law, each component of this list corresponds a line in the input file
                "Parameters": [(1000,0.3),],
                # support All or [(dim,phys), ...]
                #"ExtraDof": ("Temperature",(1.e-10, 0, 5, 0, 0, 0, 0)),
                "List of Elements": "All"
                #"List of Elements": [(3,1),]
                },
                )
        """
        raise NotImplementedError('The method getMaterialsConfigure() must be defined in Configure!')
        
        
    def getInitialBCsConfigure(self):
        """When extra-Dof is used, must be defined by a dictionary
            self.initialBCs = {
                    #"INITIAL CONDITIONS EXTRADOF": (273.,)
                    }
        """
     
