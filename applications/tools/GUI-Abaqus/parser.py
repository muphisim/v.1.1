#Parser to convert Abaqus files into MuPhiSim input files
#Logic script
#Author: Khariton Gorbunov
import math

def sameExpression(line,nodestr):#strings to be compared
    match=True
    for char1,char2 in zip(line, nodestr):
        if(char1!=char2):
            match=False
            break
    if (match):
        return True
    else:
        return False


def isTag(currentLine, tags):
    matchingTag=""
    for tag in tags:
        match=1
        for char1,char2 in zip(tag, currentLine):
            if(char1!=char2):
                match=0
                break
        if(match):
            return(1)

    return(0)

def addFromUntil(linesArray, outputArray, startkey, stopkey):
    outputArray[:]=[]
    started=False
    newArray=[]
    done=False
    for index, line in enumerate(linesArray):
        if(sameExpression(line, startkey) and done==False):
            started=True
        if(started==False):
            newArray.append(line)
        if(started):
            outputArray.append(line)
            if(sameExpression(linesArray[index+1], stopkey)):
                done=True
                started=False
    linesArray[:]=newArray
    return
def processSet(curset, label):
    nameind = curset[0].find("=")+1
    nameind
    name=""
    for i in range(nameind, len(curset[0])):
        if(curset[0][i]==","):
            break
        name+=curset[0][i]

    generate=False
    if("generate" in curset[0]):
        generate = True;
    else:
        generate = False


    parametersArray=[]
    numArr=[]
    if(generate==False):
        for i in range(1,len(curset)):
            curset[i]=curset[i].replace(" ", "")
            numRowArr=[]
            currNum=""
            for j in range (0, len(curset[i])):
                char = curset[i][j]
                if(char==","):
                    numRowArr.append(currNum)
                    currNum=""
                else:
                    currNum+=char
                if(j==len(curset[i])-1):
                    numRowArr.append(currNum)
                    currNum=""
            numArr+=numRowArr
    else:
        curset[1]=curset[1].replace(" ", "")
        numRowArr=[]
        currNum=""
        for j in range (0, len(curset[1])):
            char = curset[1][j]
            if(char==","):
                numRowArr.append(currNum)
                currNum=""
            else:
                currNum+=char
            if(j==len(curset[1])-1):
                numRowArr.append(currNum)
                currNum=""
        parametersArray+=numRowArr
        startNode=int(parametersArray[0])
        lastNode=int(parametersArray[1])
        difference=int(parametersArray[2])
        currentNumber = startNode
        while(currentNumber<=lastNode):
            numArr.append(str(currentNumber))
            currentNumber+=difference
    #sometimes people set comma and then no variable, remove
    numArrF=[val for val in numArr if val!=""]
    numArr[:]=numArrF
    setDict={"Name":name, "Values":numArr, "Label" : label}
    return(setDict)

def writeSet(outputFile, array):
    for x in array:
        outputFile.write(x+"\n")

def writeSetAndDOF(outputFile, array,DOFs, sharenodes, numberset, includeDisp=True):
    for item in array:
        for i in range(0, len(DOFs)):
            if (DOFs[i] != "?"):
                if numberset==0: ###Cases with FEM and MM nodes: only BC defined for the first Set will apllied on shared nodes
                	if(includeDisp):
                    		outputFile.write(item+","+str(i)+","+DOFs[i]+"\n")
                	else:
                    		outputFile.write(item+","+str(i)+"\n")
                elif numberset>0 and item not in sharenodes:  		               
                	if(includeDisp):
                    		outputFile.write(item+","+str(i)+","+DOFs[i]+"\n")
                	else:
                    		outputFile.write(item+","+str(i)+"\n")
                    		
                    		
def writeSetAndVec(outputFile, array,Vec):
    for item in array:
        outstring=""
        outstring+=item
        outstring+=","
        for i,component in enumerate(Vec):
            outstring+=str(component)
            if(i!=len(Vec)-1):
                outstring+=","
        outputFile.write(outstring+"\n")

def writeOutSet(namevar, label, outputFile, Sets):
    found=False
    if namevar == "All" or namevar == "all":
        outputFile.write("All\n")
    else:
        for setdict in Sets:
            if (setdict["Label"]==label):
                if(setdict["Name"]==namevar):
                    writeSet(outputFile, setdict["Values"])
                    found=True
                    break
        if(found==False):
            raise Exception("The following set of nodes does not exist")

def writeBoundary(namevar, label, outputFile, Sets, DOFs, sharednodes, numberset,includeDisp=True):
    found=False
    for setdict in Sets:
        if (setdict["Label"]==label):
            if(setdict["Name"]==namevar):
                writeSetAndDOF(outputFile, setdict["Values"],DOFs,sharednodes,numberset,includeDisp)
                found=True
                break
    if(found==False):
        raise Exception("The following set of nodes does not exist")
           

def writeForce(namevar, label, outputFile, Sets,Vec):
    found=False
    for setdict in Sets:
        if (setdict["Label"]==label):
            if(setdict["Name"]==namevar):
                writeSetAndVec(outputFile, setdict["Values"],Vec)
                found=True
                break
    if(found==False):
        raise Exception("The following set of nodes does not exist")

def separateCom(inputstring):
    inputstring=inputstring.replace(" ", "")
    numRowArr=[]
    currNum=""
    for j in range (0, len(inputstring)):
        char = inputstring[j]
        if(char==","):
            numRowArr.append(currNum)
            currNum=""
        else:
            currNum+=char
        if(j==len(inputstring)-1):
            numRowArr.append(currNum)
            currNum=""
    numRowArrF=[val for val in numRowArr if val!=""]
    numRowArr[:]=numRowArrF
    return numRowArr

def getSets(Var_inputText):
    inputName=Var_inputText


    inputFile = open(inputName,"r")
    outputArray=[]
    addingToOutput=0
    linesArrayp = inputFile.readlines()
    #remove all \n
    linesArray = [x.replace('\n', '') for x in linesArrayp]

    nsets=[]
    elsets=[]
    currentSet=[""]
    while(len(currentSet)!=0):
        currentSet=[]
        addFromUntil(linesArray, currentSet, "*Nset", "*")
        if(len(currentSet)!=0):
            nsets.append(currentSet)

    currentSet=[""]
    while(len(currentSet)!=0):
        currentSet=[]
        addFromUntil(linesArray, currentSet, "*Elset", "*")
        if(len(currentSet)!=0):
            elsets.append(currentSet)


    Sets=[]
    for i in range (0, len(nsets)):
        Sets.append(processSet(nsets[i],"Nset"))
    for i in range (0, len(elsets)):
        Sets.append(processSet(elsets[i],"Elset"))

    return Sets

def writeFileOut(Var_inputText,Var_outputText,Var_FEM_Nodes_TextEntered,Var_FEM_Elements_TextEntered,
                Var_MM_Nodes_TextEntered,Var_MM_Elements_TextEntered,NumSol,Var_SlvrTbl_SlvrType_TxtEnt,
                Var_SlvrTbl_TimeScaleFctr_TxtEnt,Var_SlvrTbl_OutFileNo_TxtEnt,Var_SlvrTbl_IntVar_TxtEnt,
                Var_SlvrTbl_Time_TxtEnt,Var_SlvrTbl_Disp_TxtEnt,Var_SlvrTbl_DOF_TxtEnt,
                NumFBCEntry,Var_FBCTbl_ForceType_TxtEnt,Var_FBCTbl_ForceSet_TxtEnt,Var_FBCTbl_TimeL_TxtEnt, Var_FBCTbl_TimeU_TxtEnt,
                Var_FBCTbl_Value_TxtEnt,Var_FBCTbl_Surface_TxtEnt,NumMLEntry,Var_MLTbl_Name_TxtEnt,Var_MLTbl_Set_TxtEnt,
                Var_MLTbl_Density_TxtEnt,Var_MLTbl_Model_TxtEnt,Var_MLTbl_Parameters_TxtEnt,Var_SlvrTbl_Forces_TxtEnt,Var_SlvrTbl_DOFForces_TxtEnt,Var_SlvrTbl_FREQ_TxtEnt,
                Var_CPUs_TxtEnt, Var_SurfaceOrder_TextEntered):
    #outputName = "output2.inp"
    #inputName = "Job-1.inp"
    outputName=Var_outputText
    inputName=Var_inputText


    outputFile= open(outputName,"w+")
    inputFile = open(inputName,"r")
    tags = ["*Node", "*Element, type="]
    outputArray=[]
    addingToOutput=0
    linesArrayp = inputFile.readlines()
    #remove all \n
    linesArray = [x.replace('\n', '') for x in linesArrayp]

    outputArrayTemp=[]
    #addFromUntil(linesArray, outputArrayTemp, "*Node", "*End Step")
    #linesArray[:]=outputArrayTemp
    finalArray=[]
    #outputArrayTemp=[]
    addFromUntil(linesArray, outputArrayTemp, "*Node", "*")
    finalArray=finalArray+outputArrayTemp
    #outputArrayTemp=[]

    addFromUntil(linesArray, outputArrayTemp, "*Element, type=", "*")
    if(Var_SurfaceOrder_TextEntered == "Linear" and ("CPS6" in outputArrayTemp[0] or "C3D10" in outputArrayTemp[0]) ):
        outputArrayTemp[0]+="_lin_surf"
    finalArray=finalArray+outputArrayTemp
    #outputArrayTemp=[]
    nsets=[]
    elsets=[]
    currentSet=[""]
    while(len(currentSet)!=0):
        currentSet=[]
        addFromUntil(linesArray, currentSet, "*Nset", "*")
        if(len(currentSet)!=0):
            nsets.append(currentSet)

    currentSet=[""]
    while(len(currentSet)!=0):
        currentSet=[]
        addFromUntil(linesArray, currentSet, "*Elset", "*")
        if(len(currentSet)!=0):
            elsets.append(currentSet)


    Sets=[]
    for i in range (0, len(nsets)):
        Sets.append(processSet(nsets[i],"Nset"))
    for i in range (0, len(elsets)):
        Sets.append(processSet(elsets[i],"Elset"))

    for line in finalArray:
        line+="\n"
        outputFile.write(line)

    FEMNodes = Var_FEM_Nodes_TextEntered#"all"
    FEMElements =Var_FEM_Elements_TextEntered #"all"#user input, set name

    MMNodes =Var_MM_Nodes_TextEntered #""
    MMElements = Var_MM_Elements_TextEntered#""
    CPUArr = Var_CPUs_TxtEnt

    a=len(FEMNodes)==0
    b=len(FEMElements)==0
    c=len(MMNodes)==0
    d=len(MMElements)==0

    conditional = ((not c) and d)or(c and (not d))or(b and d)or(a and (not b))or((not a) and b)

    if( conditional ):
        raise Exception("Please define sets for MM/FEM nodes and elements")


    if(len(FEMNodes)!=0):
        #after user sets the variable
        outputFile.write("*FEM Nodes\n")
        writeOutSet(FEMNodes, "Nset", outputFile, Sets)

    if(len(FEMElements)!=0):
        #after user sets the variable
        outputFile.write("*FEM Elements\n")
        writeOutSet(FEMElements, "Elset", outputFile, Sets)

    if(len(MMNodes)!=0):
        #after user sets the variable
        outputFile.write("*MM Nodes\n")
        writeOutSet(MMNodes, "Nset", outputFile, Sets)

    if(len(MMElements)!=0):
        #after user sets the variable
        outputFile.write("*MM Elements\n")
        writeOutSet(MMElements, "Elset", outputFile, Sets)

    if (len(FEMNodes)!=0 and len(MMNodes)!=0):
        sharedNodes=[]
        found1=False
        found2=False
        setFEM=set([])
        setMM=set([])
        for setdict in Sets:
            ###if (setdict["Label"]=="NSet"): I remove it because in other case it doesn't work
                if(setdict["Name"]==FEMNodes):
                    setFEM=set(setdict["Values"])
                    found1=True
                if(setdict["Name"]==MMNodes):
                    setMM=set(setdict["Values"])
                    found2=True
                if(found1 and found2):
                    break
        sharedNodes=list(setFEM.intersection(setMM))
    else:
         sharedNodes=[]

    if(len(CPUArr)>0):
        cpunr=1
        for i in range(0,len(CPUArr),3):
            outputFile.write("*NSET=" + CPUArr[i]+"\n")
            outputFile.write("*CPU=" +str(cpunr)+"\n")
            writeOutSet(CPUArr[i+1], "Nset", outputFile, Sets)
            cpunr+=1
        cpunr=1
        for i in range(0,len(CPUArr),3):
            outputFile.write("*ELSET=" + CPUArr[i]+"\n")
            outputFile.write("*CPU=" +str(cpunr)+"\n")
            writeOutSet(CPUArr[i+2], "Elset", outputFile, Sets)
            cpunr+=1

    numberOfSolvers = NumSol #"2" #user input
    it=0#track iterations
    #for each solver
    SolverTypeA=Var_SlvrTbl_SlvrType_TxtEnt#["IMPLICIT STATIC", "IMPLICIT DYNAMIC"]
    TimeScaleFactorA=Var_SlvrTbl_TimeScaleFctr_TxtEnt#["10","1"]
    OutputFileNumberA=Var_SlvrTbl_OutFileNo_TxtEnt#["100","1"]
    InternalVariablesA=Var_SlvrTbl_IntVar_TxtEnt#["All, max","60,23,34"]#what? not implemented.
    TimeA=Var_SlvrTbl_Time_TxtEnt#["1","2"]

    DisplacementSetsIn=Var_SlvrTbl_Disp_TxtEnt#["Boundary1, Boundary1", "Boundary1, Boundary1, Boundary1"]
    ForcesSetsIn=Var_SlvrTbl_Forces_TxtEnt
    outFrequencies=Var_SlvrTbl_FREQ_TxtEnt
    DisplacementBCsetsA=[]
    ForcessetsA=[]
    NumberOfDisplacementBCsA=[]
    NumberOfForcesA=[]
    for s in DisplacementSetsIn:
        s=s.replace(" ", "")
        xs=s.split(",")
        DisplacementBCsetsA.append(xs)#fix optionality
        NumberOfDisplacementBCsA.append(str(len(xs)))
    for s in ForcesSetsIn:
        s=s.replace(" ", "")
        xs=s.split(",")
        ForcessetsA.append(xs)
        if(len(xs)==1 and xs[0]==""):
            NumberOfForcesA.append("0")
        else:
            NumberOfForcesA.append(str(len(xs)))
    DisplacementCoordinatesIn=Var_SlvrTbl_DOF_TxtEnt#["(0, -, -), (0,-,0)","(1,1,1),(2,2,2),(3,3,3)"]
    ForcesCoordinatesIn=Var_SlvrTbl_DOFForces_TxtEnt
    coordinatesArrayFull=[]
    coordinatesArrayFullForce=[]
    for it in range(0,len(DisplacementCoordinatesIn)):
        s=DisplacementCoordinatesIn[it]
        s=s.replace(" ", "")
        startedScan=False
        coordinatesArr=[]
        currentCoor=""
        for i in range(0,len(s)):
            if(startedScan and s[i]!=")"):
                currentCoor+=s[i]

            if(s[i]=="("):
                startedScan=True
            if(startedScan and s[i]==")"):
                startedScan=False
                coordinatesArr.append(currentCoor)
                currentCoor=""
        coordinatesNArr=[]
        for item in coordinatesArr:
            tempArrItem=item.split(",")
            for it2 in range (0, len(tempArrItem)):
                if(tempArrItem[it2]=="?"):
                    tempArrItem[it2]="?"
            coordinatesNArr.append(tempArrItem)

        coordinatesArrayFull.append(coordinatesNArr)



    for it in range(0,len(ForcesCoordinatesIn)):
        s=ForcesCoordinatesIn[it]
        s=s.replace(" ", "")
        startedScan=False
        coordinatesArr=[]
        currentCoor=""
        for i in range(0,len(s)):
            if(startedScan and s[i]!=")"):
                currentCoor+=s[i]

            if(s[i]=="("):
                startedScan=True
            if(startedScan and s[i]==")"):
                startedScan=False
                coordinatesArr.append(currentCoor)
                currentCoor=""
        coordinatesNArr=[]
        for item in coordinatesArr:
            tempArrItem=item.split(",")
            for it2 in range (0, len(tempArrItem)):
                if(tempArrItem[it2]=="?"):
                    tempArrItem[it2]="?"
            coordinatesNArr.append(tempArrItem)

        coordinatesArrayFullForce.append(coordinatesNArr)



    ResetAllNeumannIn=[]#Var_SlvrTbl_Reset_TxtEnt#["1 ,3 ,4", "2, 1"]
    ResetAllNeumannTime=[]
    for s in ResetAllNeumannIn:
        s=s.replace(" ", "")
        xs=s.split(",")
        ResetAllNeumannTime.append(xs)
    NrNeumannBC=NumFBCEntry#"1"
    NeumannBCType=Var_FBCTbl_ForceType_TxtEnt#["Pressure ramp"]##all
    NeumannLowerTime=Var_FBCTbl_TimeL_TxtEnt#["0"]##these
    NeumannUpperTime=Var_FBCTbl_TimeU_TxtEnt#["1"]##array
    NeumannValue=Var_FBCTbl_Value_TxtEnt#["100000000"]##must match in length to NrNeumannBC
    NeumannSets=Var_FBCTbl_ForceSet_TxtEnt#["_Surf-1_S3"]
    NeumannSurface=Var_FBCTbl_Surface_TxtEnt#["3"]
    NrOfMaterials=NumMLEntry#"1"
    MaterialName=Var_MLTbl_Name_TxtEnt#["khariksMaterial"]
    Density=Var_MLTbl_Density_TxtEnt#["1000"]
    ConstitutiveModel=Var_MLTbl_Model_TxtEnt#["Elastic"]
    MaterialParameters=Var_MLTbl_Parameters_TxtEnt#[["1E9, 0.3"]]
    MaterialSets=Var_MLTbl_Set_TxtEnt#["all"]



    Autoscan = "no"#only working for displacments, do not implement yet


    for j in range(0,int(numberOfSolvers)):
        it=0#track iterations
        #for each solver
        SolverType=SolverTypeA[j]
        TimeScaleFactor=TimeScaleFactorA[j]
        OutputFileNumber=OutputFileNumberA[j]
        InternalVariables=InternalVariablesA[j]#what? not implemented.
        Time=TimeA[j]
        NumberOfDisplacementBCs=NumberOfDisplacementBCsA[j]


        if(len(SolverType)==0 or len(TimeScaleFactor)==0 or len(OutputFileNumber)==0 or len(Time)==0 ):
            raise Exception("Please fill in the mandatory fields")


        outputFile.write("*SOLVER, "+SolverType+"\n")
        outputFile.write("*Scale Factor\n")
        outputFile.write(TimeScaleFactor+"\n")
        outputFile.write("*Outputs\n")
        outputFile.write(OutputFileNumber+"\n")
        if(len(InternalVariables)!=0):
            outputFile.write(InternalVariables+"\n")
        if(j==0):
            initialTime="0"
            finalTime=Time
        else:
            initialTime=finalTime
            finalTime=str( float(finalTime)+float(Time) )#convert time so it atches between solvers
        outputFile.write("*Time\n")
        outputFile.write(initialTime+","+finalTime+"\n")
        outputFile.write("*BOUNDARY, TYPE=DISPLACEMENT RAMP\n")

        if(Autoscan=="no"):
            for jj in range(0, int(NumberOfDisplacementBCsA[j])):
                DisplacementBCsets=DisplacementBCsetsA[j][jj]#for each numberofdisplacementbcs
                #DisplacementX = coordinatesArrayFull[j][jj][0]
                #DisplacementY= coordinatesArrayFull[j][jj][1]
                #DisplacementZ = coordinatesArrayFull[j][jj][2]
                writeBoundary(DisplacementBCsets, "Nset", outputFile, Sets,coordinatesArrayFull[j][jj],sharedNodes, jj)
            if(int(NumberOfForcesA[j])>0):
                outputFile.write("*Forces\n")
                outputFile.write(outFrequencies[j]+"\n")
                for jj in range(0, int(NumberOfForcesA[j])):
                    Forcessets=ForcessetsA[j][jj]#for each numberofdisplacementbcs
                    #DisplacementX = coordinatesArrayFull[j][jj][0]
                    #DisplacementY= coordinatesArrayFull[j][jj][1]
                    #DisplacementZ = coordinatesArrayFull[j][jj][2]
                    writeBoundary(Forcessets, "Nset", outputFile, Sets,coordinatesArrayFullForce[j][jj],sharedNodes, jj, False)
        elif(Autoscan=="yes"):#not implemented, will error
            if(int(numberOfSolvers)!=1):
                raise Exception("Autoscan feature works only for 1 solver")
            else:
                boundaryDefinitionTemp=["dummy, will be cleared out"]
                boundaryDefinitionFin=[]
                while(len(boundaryDefinitionTemp) !=0):
                    addFromUntil(linesArray, boundaryDefinitionTemp, "*Boundary", "*")
                    if(len(boundaryDefinitionTemp)>0):
                        boundaryDefinitionTemp.pop(0)
                    for k in boundaryDefinitionTemp:
                        k=separateCom(k)
                        boundaryDefinitionFin.append(k)
                for jj in range(0,len(boundaryDefinitionFin)):
                    dispVal="0"
                    xyz=["","",""]
                    if (len(boundaryDefinitionFin[jj])==4):
                        dispVal=boundaryDefinitionFin[jj][3]
                    if(boundaryDefinitionFin[jj][1] != boundaryDefinitionFin[jj][2]):
                        raise Exception("Something went wrong. Please input  the boundary conditions manually")
                    xyz[ int(boundaryDefinitionFin[jj][1]) - 1]=dispVal
                    writeBoundary(boundaryDefinitionFin[jj][0], "Nset", outputFile, Sets,xyz[0],xyz[1],xyz[2],sharedNodes)
        outputFile.write("*END SOLVER\n")
    for j in range(0, int(NrNeumannBC) ):
        if(len(NeumannBCType[j]))>0:
            pressure=False
            forcevec=[]
            if(NeumannBCType[j].upper() == "pressure inst".upper() or NeumannBCType[j].upper()=="pressure ramp".upper() or NeumannBCType[j].upper()=="Heat inst".upper()):
                pressure=True
            if(not pressure):
                forcevecstr=NeumannValue[j].replace(" ", "")
                forcevecstr=forcevecstr.replace("(","")
                forcevecstr=forcevecstr.replace(")","")
                forcevec=forcevecstr.split(",")
                forcemag=0
                for component in forcevec:
                    forcemag+=float(component)**2
                forcemag=math.sqrt(forcemag)
                NeumannValue[j]=str(forcemag)
                for i in range (0, len(forcevec) ):
                    forcevec[i]=float(forcevec[i])/forcemag
            if(NeumannBCType[j].upper()=="Heat inst".upper()): 
            	 NeumannValue[j]=NeumannValue[j].replace("(","") 
            	 NeumannValue[j]=NeumannValue[j].replace(")","")
            if(NeumannBCType[j].upper() == "Current inst".upper() or NeumannBCType[j].upper()=="Heat inst".upper() or NeumannBCType[j].upper() == "Volumetric heat flux".upper()): 
            	outputFile.write("*FLUX EXTRADOF, TYPE="+NeumannBCType[j].upper() +"\n")
            else:
            	outputFile.write("*BOUNDARY, TYPE="+NeumannBCType[j].upper() +"\n")
            parametersLine=NeumannLowerTime[j] +", "+NeumannUpperTime[j]+", "+NeumannValue[j]
            if(pressure):
                parametersLine=parametersLine+", "+NeumannSurface[0]+"\n"
                NeumannSurface.pop(0)#keep track of pressures vs nodes
            else:
                parametersLine+="\n"
            outputFile.write(parametersLine)
            if(pressure):
                writeOutSet(NeumannSets[j], "Elset", outputFile, Sets) 
            elif (NeumannBCType[j].upper() == "Current inst".upper() or NeumannBCType[j].upper() == "Volumetric heat flux".upper()): 
                writeOutSet(NeumannSets[j], "Elset", outputFile, Sets) 
            else:
                #writeOutSet(NeumannSets[j], "Nset", outputFile, Sets)
                writeForce(NeumannSets[j], "Nset", outputFile, Sets,forcevec)


    for j in range (0, int(NrOfMaterials) ):
    	if(ConstitutiveModel[j] == "FHNelec" or ConstitutiveModel[j] == "Temperature" or ConstitutiveModel[j] == "3DprintingTemperature"): 
        	outputFile.write("*ExtraDof\n")
        	outputFile.write("*"+ConstitutiveModel[j]+"\n")#.replace(" ", "")
        	for k in range(0, len(MaterialParameters[j])):
            		outputFile.write(MaterialParameters[j][k]+"\n")
    	else:
        	outputFile.write("*MATERIAL, name="+MaterialName[j].replace(" ", "")+"\n")
        	outputFile.write("*Density\n")
        	outputFile.write(Density[j]+"\n")
        	outputFile.write("*"+ConstitutiveModel[j]+"\n")#.replace(" ", "")
        	for k in range(0, len(MaterialParameters[j])):
            		outputFile.write(MaterialParameters[j][k]+"\n")
    	outputFile.write("*List of Elements\n")
    	if(MaterialSets[j].upper()=="ALL"):
    		outputFile.write("All\n")
    	else:
    		writeOutSet(MaterialSets[j], "Elset", outputFile, Sets)
    outputFile.write("*END\n")
