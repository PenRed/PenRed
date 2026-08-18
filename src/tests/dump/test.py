import pyPenred
import numpy as np
import os
import shutil
from pathlib import Path

currentFolder = Path(".")
folder1 = currentFolder / "res1"
folder2 = currentFolder / "res2"
folder3 = currentFolder / "res3"

folder1.mkdir(exist_ok=True)
folder2.mkdir(exist_ok=True)
folder3.mkdir(exist_ok=True)

sim1 = pyPenred.runFromFile("config1.in")

for f in currentFolder.glob("*.dat"):
    shutil.move(str(f), str(folder1 / f.name))
    
#sim2 = pyPenred.runFromFile("config2.in")
sim2 = pyPenred.simulation.create()
sim2.configFromFile("config2.in")

sim2.addDumps(["th0test.dump", "th1test.dump", "th3test.dump"])
for f in currentFolder.glob("*.dat"):
    shutil.move(str(f), str(folder2 / f.name))

    
sim2.addDumps(["th0test.dump", "th1test.dump", "th2test.dump", "th3test.dump"])
for f in currentFolder.glob("*.dat"):
    shutil.move(str(f), str(folder3 / f.name))
