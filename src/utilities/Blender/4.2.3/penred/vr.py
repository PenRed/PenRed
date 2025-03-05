
def createInteractionForcing(f, name, partType, interactionIndex, factorType, factor, weightWindow, bodies):
    f.write(f"# Interaction Forcing configuration for '{name}'\n")
    f.write(f"VR/IForcing/{name}/particle \"{partType}\"\n")
    f.write(f"VR/IForcing/{name}/interaction {interactionIndex}\n")

    if factorType == "AVERAGE":
        f.write(f"VR/IForcing/{name}/factor -{factor}\n")
    else:
        f.write(f"VR/IForcing/{name}/factor  {factor}\n")

    f.write(f"VR/IForcing/{name}/min-weight {weightWindow[0]:.4e}\n")
    f.write(f"VR/IForcing/{name}/max-weight {weightWindow[1]:.4e}\n")

    for body in bodies:
        f.write(f"VR/IForcing/{name}/bodies/{body} true\n")
    
    f.write("\n")

def createBremssSplitting(f, name, factor, bodies):
    f.write(f"# Bremsstrahlung splitting configuration for '{name}'\n")
    f.write(f"VR/bremss/{name}/splitting {factor}\n")
    for body in bodies:
        f.write(f"VR/bremss/{name}/bodies/{body} true\n")    
    f.write("\n")

def createXRaySplitting(f, name, bodies):
    f.write(f"# X-Ray splitting configuration for '{name}'\n")
    f.write(f"VR/photon/{name}/type \"XRAY_SPLITTING\"\n")
    for body in bodies:
        f.write(f"VR/photon/{name}/bodies/{body[0]}/splitting {body[1]}\n")    
    f.write("\n")
