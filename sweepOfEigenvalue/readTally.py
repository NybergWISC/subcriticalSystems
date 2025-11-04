import sys
import numpy as np
import matplotlib.pyplot as plt



class ReadTallyOut:

   def __init__(self,filename):

      with open(filename,"r") as f:
         data = f.readlines()

      self.statepoints=[]
      first_state = True
      sp = None
      sp = dict()

      for j,line in enumerate(data):

         if "heating" in line:
            lineArray = line.split()
            if len(lineArray) != 4:
               print("Error: incorrect line read")
            else:
               sp["heating_val"] = lineArray[1]
               sp["heating_error"] = lineArray[3]
         if " Fission Rate" in line:
            lineArray = line.split()
            if len(lineArray) != 5:
               print("Error: incorrect line read")
            else:
               sp["fission_rate"] = lineArray[2]
               sp["fission_error"] = lineArray[4]
         if "Nu-Fission Rate" in line:
            lineArray = line.split()
            if len(lineArray) != 5:
               print("Error: incorrect line read")
            else:
               sp["nu_fission_rate"] = lineArray[2]
               sp["nu_fission_error"] = lineArray[4]
         if "Xt" in line:
            lineArray = line.split()
            if len(lineArray) != 4:
               print("Error: incorrect line read")
               print("lineArray[1] = " + str(lineArray[1]))
               print("lineArray[3] = " + str(lineArray[3]))
            else:
               sp["trit_val"] = lineArray[1]
               sp["trit_error"] = lineArray[3]
      self.statepoints.append(sp)

def funcReadTally(parameter_label):
    path='./tallies.out'
    heatDict = ReadTallyOut(path)
    print(heatDict.statepoints[0]["heating_val"])
    with open('output_verbose.txt', 'a') as f:
        f.write(" Param Value: "+ str(parameter_label))
        f.write(" Heating: "+ str(heatDict.statepoints[0]["heating_val"]))
        f.write(" Error: "+ str(heatDict.statepoints[0]["heating_error"]))
        f.write(" Fission Rate: "+ str(heatDict.statepoints[0]["fission_rate"]))
        f.write(" Error: "+ str(heatDict.statepoints[0]["fission_error"]))
        f.write(" Nu-Fission Rate: "+ str(heatDict.statepoints[0]["nu_fission_rate"]))
        f.write(" Error: "+ str(heatDict.statepoints[0]["nu_fission_error"]))        
        f.write(" Tritium Breeding Rate: "+ str(heatDict.statepoints[0]["trit_val"]))
        f.write(" Error: "+ str(heatDict.statepoints[0]["trit_error"] + "\n"))
    with open('output.txt', 'a') as g:
        g.write(str(parameter_label) + ", ")
        g.write(str(heatDict.statepoints[0]["heating_val"]) + ", ")
        g.write(str(heatDict.statepoints[0]["heating_error"]) + ", ")
        g.write(str(heatDict.statepoints[0]["fission_rate"]) + ", ")
        g.write(str(heatDict.statepoints[0]["fission_error"]) + ", ")
        g.write(str(heatDict.statepoints[0]["nu_fission_rate"]) + ", ")
        g.write(str(heatDict.statepoints[0]["nu_fission_error"]) + ", ")        
        g.write(str(heatDict.statepoints[0]["trit_val"]) + ", ")
        g.write(str(heatDict.statepoints[0]["trit_error"]) + "\n")

funcReadTally(sys.argv[1])
